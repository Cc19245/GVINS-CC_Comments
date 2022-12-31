#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <gvins/LocalSensorExternalTrigger.h>
#include <sensor_msgs/NavSatFix.h>

#include "estimator.h"
#include "parameters.h"
#include "utility/visualization.h"

using namespace gnss_comm;

//! 疑问：这是啥？
//; 解答：这个是在寻找GNSS和VIO对齐的时间的时候，在靠前0.05秒找，因为不能把离当前太近的都删掉了
#define MAX_GNSS_CAMERA_DELAY 0.05  

std::unique_ptr<Estimator> estimator_ptr;

std::condition_variable con;
double current_time = -1;
queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<std::vector<ObsPtr>> gnss_meas_buf;  //; buf中每个位置存的都是一个vector的观测，就是每个位置都会有多个卫星的观测
queue<sensor_msgs::PointCloudConstPtr> relo_buf;
int sum_of_wait = 0;

std::mutex m_buf;
std::mutex m_state;
std::mutex i_buf;
std::mutex m_estimator;

double latest_time;
Eigen::Vector3d tmp_P;
Eigen::Quaterniond tmp_Q;
Eigen::Vector3d tmp_V;
Eigen::Vector3d tmp_Ba;
Eigen::Vector3d tmp_Bg;
Eigen::Vector3d acc_0;
Eigen::Vector3d gyr_0;
bool init_feature = 0;
bool init_imu = 1;
double last_imu_t = -1;

std::mutex m_time;
double next_pulse_time;      // PPS触发时间
bool next_pulse_time_valid;  // 如果进入PPS触发的回调函数，那么这个就会设置成true
double time_diff_gnss_local; // PPS触发时间和VI传感器实际被触发的时间之间的差值
bool time_diff_valid;   // 如果这个是false，则对于收到的gnss观测数据直接不会存储
double latest_gnss_time;
double tmp_last_feature_time;
uint64_t feature_msg_counter;
int skip_parameter;

void predict(const sensor_msgs::ImuConstPtr &imu_msg)
{
    double t = imu_msg->header.stamp.toSec();
    if (init_imu)
    {
        latest_time = t;
        init_imu = 0;
        return;
    }
    double dt = t - latest_time;
    latest_time = t;

    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    Eigen::Vector3d linear_acceleration{dx, dy, dz};

    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Eigen::Vector3d angular_velocity{rx, ry, rz};

    Eigen::Vector3d un_acc_0 = tmp_Q * (acc_0 - tmp_Ba) - estimator_ptr->g;

    Eigen::Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - tmp_Bg;
    tmp_Q = tmp_Q * Utility::deltaQ(un_gyr * dt);

    Eigen::Vector3d un_acc_1 = tmp_Q * (linear_acceleration - tmp_Ba) - estimator_ptr->g;

    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);

    tmp_P = tmp_P + dt * tmp_V + 0.5 * dt * dt * un_acc;
    tmp_V = tmp_V + dt * un_acc;

    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

void update()
{
    TicToc t_predict;
    latest_time = current_time;
    tmp_P = estimator_ptr->Ps[WINDOW_SIZE];
    tmp_Q = estimator_ptr->Rs[WINDOW_SIZE];
    tmp_V = estimator_ptr->Vs[WINDOW_SIZE];
    tmp_Ba = estimator_ptr->Bas[WINDOW_SIZE];
    tmp_Bg = estimator_ptr->Bgs[WINDOW_SIZE];
    acc_0 = estimator_ptr->acc_0;
    gyr_0 = estimator_ptr->gyr_0;

    queue<sensor_msgs::ImuConstPtr> tmp_imu_buf = imu_buf;
    for (sensor_msgs::ImuConstPtr tmp_imu_msg; !tmp_imu_buf.empty(); tmp_imu_buf.pop())
        predict(tmp_imu_buf.front());

}

// 同步一帧图像和多个IMU、GNSS观测的数据
bool getMeasurements(std::vector<sensor_msgs::ImuConstPtr> &imu_msg, 
    sensor_msgs::PointCloudConstPtr &img_msg, std::vector<ObsPtr> &gnss_msg)
{
    //; 注意这个地方很有意思，是按照顺序进行或的，也就是如果imu不是空，这里就能过去。
    //; 所以如果gnss一直收不到，那么也不影响单纯的VIO运行
    if (imu_buf.empty() || feature_buf.empty() || (GNSS_ENABLE && gnss_meas_buf.empty()))
        return false;
    
    double front_feature_ts = feature_buf.front()->header.stamp.toSec();

    // 最新的IMU时间比图像时间还早，说明imu还没到
    if (!(imu_buf.back()->header.stamp.toSec() > front_feature_ts))
    {
        //ROS_WARN("wait for imu, only should happen at the beginning");
        sum_of_wait++;
        return false;
    }
    double front_imu_ts = imu_buf.front()->header.stamp.toSec();
    // 最老的图像时间比最老的IMU时间还老，那么只能把图像丢掉
    while (!feature_buf.empty() && front_imu_ts > front_feature_ts)
    {
        ROS_WARN("throw img, only should happen at the beginning");
        feature_buf.pop();
        front_feature_ts = feature_buf.front()->header.stamp.toSec();
    }
    //; ----------- 至此，就找到了和IMU能够对齐的图像时间 


    if (GNSS_ENABLE)
    {
        front_feature_ts += time_diff_gnss_local;  // 补偿图像时间，和GNSS时间对齐
        double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
        // 把太老的GNSS数据全部丢掉
        while (!gnss_meas_buf.empty() && front_gnss_ts < front_feature_ts-MAX_GNSS_CAMERA_DELAY)
        {
            ROS_WARN("throw gnss, only should happen at the beginning");
            gnss_meas_buf.pop();
            if (gnss_meas_buf.empty()) return false;
            front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
        }
        //! 疑问：如果是在室内，此时GNSS全部失效，那这里返回false，岂不是VIO都不能运行了？
        if (gnss_meas_buf.empty())
        {
            ROS_WARN("wait for gnss...");
            return false;
        }
        //; 如果在时间容忍范围内，则这个gnss观测数据就可以使用
        //! 疑问：但是为什么还是用了一个时间容忍来寻找呢?
        else if (abs(front_gnss_ts-front_feature_ts) < MAX_GNSS_CAMERA_DELAY)
        {
            gnss_msg = gnss_meas_buf.front();
            gnss_meas_buf.pop();
        }
    }

    img_msg = feature_buf.front();
    feature_buf.pop();

    // 最后，把所有可用的IMU序列找出来
    while (imu_buf.front()->header.stamp.toSec() < img_msg->header.stamp.toSec() + estimator_ptr->td)
    {
        imu_msg.emplace_back(imu_buf.front());
        imu_buf.pop();
    }
    imu_msg.emplace_back(imu_buf.front());
    if (imu_msg.empty())
        ROS_WARN("no imu between two image");
    return true;
}

/* imu消息存进buffer，同时按照imu频率预测predict位姿并发送(IMU状态递推并发布[P,Q,V,header])，这样就可以提高里程计频率*/
void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg)
{
    if (imu_msg->header.stamp.toSec() <= last_imu_t)
    {
        ROS_WARN("imu message in disorder!");
        return;
    }
    last_imu_t = imu_msg->header.stamp.toSec();

    m_buf.lock();
    imu_buf.push(imu_msg);
    m_buf.unlock();
    con.notify_one();

    last_imu_t = imu_msg->header.stamp.toSec();

    {
        std::lock_guard<std::mutex> lg(m_state);
        predict(imu_msg);
        std_msgs::Header header = imu_msg->header;
        header.frame_id = "world";
        if (estimator_ptr->solver_flag == Estimator::SolverFlag::NON_LINEAR)
            pubLatestOdometry(tmp_P, tmp_Q, tmp_V, header);
    }
}


/**
 * @brief 订阅GPS, Galileo, BeiDou 星历信息
 *   1.把ROS消息转成星历Ephem的数据结构
 *   2.把星历的数据结构存储到estimator的成员变量中
 * 
 * @param[in] ephem_msg 
 */
void gnss_ephem_callback(const GnssEphemMsgConstPtr &ephem_msg)
{
    EphemPtr ephem = msg2ephem(ephem_msg); //将ROS星历信息，转换成相应的Ephem数据
    estimator_ptr->inputEphem(ephem);
}

/* 同上，订阅GLONASS ephemeris，应该是与另外三导航系统的星历格式不同*/
void gnss_glo_ephem_callback(const GnssGloEphemMsgConstPtr &glo_ephem_msg)
{
    GloEphemPtr glo_ephem = msg2glo_ephem(glo_ephem_msg);
    estimator_ptr->inputEphem(glo_ephem);
}


/* 订阅GNSS broadcast ionospheric parameters，电离层参数*/
//卫星信号在传播的过程中会受到电离层和对流层的影响，且如果建模不正确或不考虑两者的影响，会导致定位结果变差，
// 因此，通常都会对两者进行建模处理；后面我们在选择卫星信号时，会考虑卫星的仰角，也是因为对于仰角小的卫星，
// 其信号在电离层和对流层中经过的时间较长，对定位影响大，这样的卫星我们就会排除
void gnss_iono_params_callback(const StampedFloat64ArrayConstPtr &iono_msg)
{
    double ts = iono_msg->header.stamp.toSec();
    std::vector<double> iono_params;
    std::copy(iono_msg->data.begin(), iono_msg->data.end(), std::back_inserter(iono_params));
    //; 注意：这里可以看出来电离层参数一定是8个，具体应该和内部的定义有关
    assert(iono_params.size() == 8);
    estimator_ptr->inputIonoParams(ts, iono_params); //更新电离层的参数
}

/**
 * @brief 订阅 GNSS 的原始 measurements，包含Code Pseudorange 和Doppler Measurement
 * 
 */
void gnss_meas_callback(const GnssMeasMsgConstPtr &meas_msg)
{
    //; 解析gnss观测信息，注意这个和上面的星历信息有些不同，上面星历信息每次就是1个，
    //! 这里的观测信息一次有多个，为什么？
    std::vector<ObsPtr> gnss_meas = msg2meas(meas_msg); // 从ros消息解析GNSS测量值

    latest_gnss_time = time2sec(gnss_meas[0]->time);    // 将时间转换为秒,such as x.y秒

    // cerr << "gnss ts is " << std::setprecision(20) << time2sec(gnss_meas[0]->time) << endl;
    if (!time_diff_valid)   
        return;

    m_buf.lock();
    //TODO： std::move! 还是要好好看看
    gnss_meas_buf.push(std::move(gnss_meas)); //将解析之后GNSS测量push到queue中，并且gnss_meas_buf是一个全局变量
    m_buf.unlock();
    con.notify_one();
}

/* feature回调函数，将feature_msg放入feature_buf*/
void feature_callback(const sensor_msgs::PointCloudConstPtr &feature_msg)
{
    ++ feature_msg_counter;

    if (skip_parameter < 0 && time_diff_valid)
    {
        const double this_feature_ts = feature_msg->header.stamp.toSec()+time_diff_gnss_local;
        if (latest_gnss_time > 0 && tmp_last_feature_time > 0)
        {
            if (abs(this_feature_ts - latest_gnss_time) > abs(tmp_last_feature_time - latest_gnss_time))
                skip_parameter = feature_msg_counter%2;       // skip this frame and afterwards
            else
                skip_parameter = 1 - (feature_msg_counter%2);   // skip next frame and afterwards
        }
        // cerr << "feature counter is " << feature_msg_counter << ", skip parameter is " << int(skip_parameter) << endl;
        tmp_last_feature_time = this_feature_ts;
    }

    if (skip_parameter >= 0 && int(feature_msg_counter%2) != skip_parameter)
    {
        m_buf.lock();
        feature_buf.push(feature_msg);
        m_buf.unlock();
        con.notify_one();
    }
}

/* 订阅external trigger info of the local sensor，VI传感器的外部触发信息*/
/** trigger_msg记录的是VI传感器被GNSS脉冲触发时的时间，也可以理解成图像的命名（以时间命名），
 * 这个和真正的gnss时间是有区别的，因为存在硬件延迟等，这也是后面为什么会校正local和world时间的原因*/
void local_trigger_info_callback(const gvins::LocalSensorExternalTriggerConstPtr &trigger_msg)
{
    std::lock_guard<std::mutex> lg(m_time);

    //; 如果之前记录过PPS触发时间
    if (next_pulse_time_valid)
    {
        // next_pulse_time记录了PPS触发时间
        // trigger_msg记录了VI传感器实际被触发的时间
        // time_diff_gnss_local表示PPS触发时间和VI传感器实际被触发的时间之间的差值，应该是由于硬件传输延迟等原因
        time_diff_gnss_local = next_pulse_time - trigger_msg->header.stamp.toSec();
        estimator_ptr->inputGNSSTimeDiff(time_diff_gnss_local);
        if (!time_diff_valid)       // just get calibrated
            std::cout << "time difference between GNSS and VI-Sensor got calibrated: "
                << std::setprecision(15) << time_diff_gnss_local << " s\n";
        time_diff_valid = true;  //; 标记PPS和本地VI时间插值已经被标定过
    }
}

/**
 * @brief GNSS接受机的PPS触发信号的回调函数，内部存储PPS的触发时间
 * 
 * @param[in] tp_msg 
 */
void gnss_tp_info_callback(const GnssTimePulseInfoMsgConstPtr &tp_msg)
{
    // 先把GPS时间转成gtime_t数据结构
    gtime_t tp_time = gpst2time(tp_msg->time.week, tp_msg->time.tow); //convert week and tow in gps time to gtime_t struct
    // 根据不同的卫星系统，将gps时间进一步处理
    if (tp_msg->utc_based || tp_msg->time_sys == SYS_GLO)
        tp_time = utc2gpst(tp_time);
    else if (tp_msg->time_sys == SYS_GAL)
        tp_time = gst2time(tp_msg->time.week, tp_msg->time.tow);
    else if (tp_msg->time_sys == SYS_BDS)
        tp_time = bdt2time(tp_msg->time.week, tp_msg->time.tow);
    else if (tp_msg->time_sys == SYS_NONE)
    {
        std::cerr << "Unknown time system in GNSSTimePulseInfoMsg.\n";
        return;
    }
    double gnss_ts = time2sec(tp_time);

    std::lock_guard<std::mutex> lg(m_time);
    next_pulse_time = gnss_ts; // 记录PPS触发时间
    next_pulse_time_valid = true;
}

void restart_callback(const std_msgs::BoolConstPtr &restart_msg)
{
    if (restart_msg->data == true)
    {
        ROS_WARN("restart the estimator!");
        m_buf.lock();
        while(!feature_buf.empty())
            feature_buf.pop();
        while(!imu_buf.empty())
            imu_buf.pop();
        m_buf.unlock();
        m_estimator.lock();
        estimator_ptr->clearState();
        estimator_ptr->setParameter();
        m_estimator.unlock();
        current_time = -1;
        last_imu_t = 0;
    }
    return;
}

/**
 * @brief 程序的主要入口，初始化，及后续优化均在这里调用
 * 
 */
void process()
{
    while (true) 
    {
        std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>> measurements;
        std::vector<sensor_msgs::ImuConstPtr> imu_msg;
        sensor_msgs::PointCloudConstPtr img_msg;
        std::vector<ObsPtr> gnss_msg;   // gnss的观测信息，对于一个图像帧，会有多个卫星的观测，因此是vector

        // Step 1. 同步IMU、图像和GNSS数据
        std::unique_lock<std::mutex> lk(m_buf);
        con.wait(lk, [&]
                 {
                    //; 这帧图像和上一帧图像之间包括：一帧图像特征点、多个IMU数据、多个gnss数据
                    return getMeasurements(imu_msg, img_msg, gnss_msg);
                 });
        lk.unlock();
        m_estimator.lock();
        double dx = 0, dy = 0, dz = 0, rx = 0, ry = 0, rz = 0;

        // Step 2. 执行IMU预积分
        for (auto &imu_data : imu_msg)
        {
            double t = imu_data->header.stamp.toSec();
            double img_t = img_msg->header.stamp.toSec() + estimator_ptr->td;
            if (t <= img_t)
            { 
                if (current_time < 0)
                    current_time = t;
                double dt = t - current_time;
                ROS_ASSERT(dt >= 0);
                current_time = t;
                dx = imu_data->linear_acceleration.x;
                dy = imu_data->linear_acceleration.y;
                dz = imu_data->linear_acceleration.z;
                rx = imu_data->angular_velocity.x;
                ry = imu_data->angular_velocity.y;
                rz = imu_data->angular_velocity.z;
                /*处理IMU数据，包括更新预积分量，和通过中值积分得到当前PQV，为后端优化提供优化状态量的初始值*/
                estimator_ptr->processIMU(dt, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz)); /* processIMU <<<<<<<<<<------------------------*/
                //printf("imu: dt:%f a: %f %f %f w: %f %f %f\n",dt, dx, dy, dz, rx, ry, rz);
            }
            else //针对最后一个imu数据，做一个简单的线性插值
            {
                double dt_1 = img_t - current_time;
                double dt_2 = t - img_t;
                current_time = img_t;
                ROS_ASSERT(dt_1 >= 0);
                ROS_ASSERT(dt_2 >= 0);
                ROS_ASSERT(dt_1 + dt_2 > 0);
                double w1 = dt_2 / (dt_1 + dt_2);
                double w2 = dt_1 / (dt_1 + dt_2);
                dx = w1 * dx + w2 * imu_data->linear_acceleration.x;
                dy = w1 * dy + w2 * imu_data->linear_acceleration.y;
                dz = w1 * dz + w2 * imu_data->linear_acceleration.z;
                rx = w1 * rx + w2 * imu_data->angular_velocity.x;
                ry = w1 * ry + w2 * imu_data->angular_velocity.y;
                rz = w1 * rz + w2 * imu_data->angular_velocity.z;
                estimator_ptr->processIMU(dt_1, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
                //printf("dimu: dt:%f a: %f %f %f w: %f %f %f\n",dt_1, dx, dy, dz, rx, ry, rz);
            }
        }

        // Step 3. 处理GNSS观测和星历信息，放到estimator的类成员变量中
        //; 注意在处理过程中会对GNSS的观测和星历信息进行过滤，只留下满足要求、质量较好的信息
        if (GNSS_ENABLE && !gnss_msg.empty())
            estimator_ptr->processGNSS(gnss_msg);

        ROS_DEBUG("processing vision data with stamp %f \n", img_msg->header.stamp.toSec());

        // Step 4. 统计前端的特征点追踪信息
        TicToc t_s;
        map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> image;
        for (unsigned int i = 0; i < img_msg->points.size(); i++)
        {
            int v = img_msg->channels[0].values[i] + 0.5;
            int feature_id = v / NUM_OF_CAM;
            int camera_id = v % NUM_OF_CAM;
            double x = img_msg->points[i].x;
            double y = img_msg->points[i].y;
            double z = img_msg->points[i].z;
            double p_u = img_msg->channels[1].values[i];
            double p_v = img_msg->channels[2].values[i];
            double velocity_x = img_msg->channels[3].values[i];
            double velocity_y = img_msg->channels[4].values[i];
            ROS_ASSERT(z == 1);
            Eigen::Matrix<double, 7, 1> xyz_uv_velocity;
            xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
            image[feature_id].emplace_back(camera_id,  xyz_uv_velocity);
        }

        // Step 5. 重点：后端优化
        estimator_ptr->processImage(image, img_msg->header);

        // Step 6. 一次处理完成，进行一些统计信息计算
        double whole_t = t_s.toc();
        printStatistics(*estimator_ptr, whole_t);
        std_msgs::Header header = img_msg->header;
        header.frame_id = "world";

        pubOdometry(*estimator_ptr, header);
        pubKeyPoses(*estimator_ptr, header);
        pubCameraPose(*estimator_ptr, header);
        pubPointCloud(*estimator_ptr, header);
        pubTF(*estimator_ptr, header);
        pubKeyframe(*estimator_ptr);
        m_estimator.unlock();
        m_buf.lock();
        m_state.lock();
        if (estimator_ptr->solver_flag == Estimator::SolverFlag::NON_LINEAR)
            update();
        m_state.unlock();
        m_buf.unlock();
    }
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "gvins");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);
    readParameters(n);
    estimator_ptr.reset(new Estimator());
    estimator_ptr->setParameter();
#ifdef EIGEN_DONT_PARALLELIZE
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif

    registerPub(n);

    next_pulse_time_valid = false;
    time_diff_valid = false;
    latest_gnss_time = -1;
    tmp_last_feature_time = -1;
    feature_msg_counter = 0;

    if (GNSS_ENABLE)
        skip_parameter = -1;
    else
        skip_parameter = 0;

    ros::Subscriber sub_imu = n.subscribe(IMU_TOPIC, 2000, 
        imu_callback, ros::TransportHints().tcpNoDelay());
    ros::Subscriber sub_feature = n.subscribe("/gvins_feature_tracker/feature", 2000, 
        feature_callback);
    ros::Subscriber sub_restart = n.subscribe("/gvins_feature_tracker/restart", 2000, 
        restart_callback);

    //? add: gnss message
    ros::Subscriber sub_ephem, sub_glo_ephem, sub_gnss_meas, sub_gnss_iono_params;
    ros::Subscriber sub_gnss_time_pluse_info, sub_local_trigger_info;
    if (GNSS_ENABLE)
    {
        //; 1.订阅星历信息：卫星的位置、速度、时间偏差等信息
        sub_ephem = n.subscribe(GNSS_EPHEM_TOPIC, 100, gnss_ephem_callback); //GPS, Galileo, BeiDou ephemeris
        sub_glo_ephem = n.subscribe(GNSS_GLO_EPHEM_TOPIC, 100, gnss_glo_ephem_callback); //GLONASS ephemeris

        //; 2.订阅卫星的观测信息
        sub_gnss_meas = n.subscribe(GNSS_MEAS_TOPIC, 100, gnss_meas_callback); //GNSS raw measurement topic
        
        //; 3.订阅电离层延时相关信息
        sub_gnss_iono_params = n.subscribe(GNSS_IONO_PARAMS_TOPIC, 100, gnss_iono_params_callback); //GNSS broadcast ionospheric parameters，电离层参数

        // GNSS与VIO的时间是否在线同步
        //TODO: 看看下面两个时间同步的消息的频率是多少？
        if (GNSS_LOCAL_ONLINE_SYNC)    // 在线同步
        {
            // The PPS signal from the GNSS receiver is used to trigger the VI-Sensor to align the global time with the local time.
            // GNSS接收机的PPS信号被用来触发VI传感器
            sub_gnss_time_pluse_info = n.subscribe(GNSS_TP_INFO_TOPIC, 100, 
                gnss_tp_info_callback);  // PPS time info

            // external trigger info of the local sensor，VI传感器的外部触发信息
            sub_local_trigger_info = n.subscribe(LOCAL_TRIGGER_INFO_TOPIC, 100, 
                local_trigger_info_callback); 
        }
        else   // 如果不需要在线同步，那么直接设置成固定值即可
        {
            // GNSS_LOCAL_TIME_DIFF = fsSettings["gnss_local_time_diff"];
            time_diff_gnss_local = GNSS_LOCAL_TIME_DIFF;
            estimator_ptr->inputGNSSTimeDiff(time_diff_gnss_local);
            time_diff_valid = true;
        }
    }

    std::thread measurement_process{process};
    ros::spin();

    return 0;
}
