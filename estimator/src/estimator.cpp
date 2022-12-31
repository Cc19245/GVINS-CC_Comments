#include "estimator.h"

Estimator::Estimator(): f_manager{Rs}
{
    ROS_INFO("init begins");
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
        pre_integrations[i] = nullptr;
    clearState();
}

void Estimator::setParameter()
{
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = TIC[i];
        ric[i] = RIC[i];
    }
    f_manager.setRic(ric);
    ProjectionFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    ProjectionTdFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    td = TD;
}

void Estimator::clearState()
{
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();

        if (pre_integrations[i] != nullptr)
            delete pre_integrations[i];
        pre_integrations[i] = nullptr;
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d::Zero();
        ric[i] = Matrix3d::Identity();
    }

    for (auto &it : all_image_frame)
    {
        if (it.second.pre_integration != nullptr)
        {
            delete it.second.pre_integration;
            it.second.pre_integration = nullptr;
        }
    }

    solver_flag = INITIAL;
    first_imu = false,
    sum_of_back = 0;
    sum_of_front = 0;
    frame_count = 0;
    solver_flag = INITIAL;
    initial_timestamp = 0;
    all_image_frame.clear();
    td = TD;

    gnss_ready = false;
    anc_ecef.setZero();
    R_ecef_enu.setIdentity();
    para_yaw_enu_local[0] = 0;
    yaw_enu_local = 0;
    sat2ephem.clear();
    sat2time_index.clear();
    sat_track_status.clear();
    latest_gnss_iono_params.clear();
    std::copy(GNSS_IONO_DEFAULT_PARAMS.begin(), GNSS_IONO_DEFAULT_PARAMS.end(), 
        std::back_inserter(latest_gnss_iono_params));
    diff_t_gnss_local = 0;

    first_optimization = true;

    if (tmp_pre_integration != nullptr)
        delete tmp_pre_integration;
    if (last_marginalization_info != nullptr)
        delete last_marginalization_info;

    tmp_pre_integration = nullptr;
    last_marginalization_info = nullptr;
    last_marginalization_parameter_blocks.clear();

    f_manager.clearState();

    failure_occur = 0;
}

void Estimator::processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
    if (!first_imu)
    {
        first_imu = true;
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
    }

    if (!pre_integrations[frame_count])
    {
        pre_integrations[frame_count] = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};
    }
    if (frame_count != 0)
    {
        pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity);
        //if(solver_flag != NON_LINEAR)
            tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity);

        dt_buf[frame_count].push_back(dt);
        linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        angular_velocity_buf[frame_count].push_back(angular_velocity);

        int j = frame_count;         
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g;
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - g;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

void Estimator::processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, const std_msgs::Header &header)
{
    ROS_DEBUG("new image coming ------------------------------------------");
    ROS_DEBUG("Adding feature points %lu", image.size());
    // Step 1. 根据当前帧状态，判断边缘化老帧还是次新帧
    if (f_manager.addFeatureCheckParallax(frame_count, image, td))
        marginalization_flag = MARGIN_OLD;
    else
        marginalization_flag = MARGIN_SECOND_NEW;

    ROS_DEBUG("this frame is--------------------%s", marginalization_flag ? "reject" : "accept");
    ROS_DEBUG("%s", marginalization_flag ? "Non-keyframe" : "Keyframe");
    ROS_DEBUG("Solving %d", frame_count);
    ROS_DEBUG("number of feature: %d", f_manager.getFeatureCount());// 有效特征点数量
    Headers[frame_count] = header;

    ImageFrame imageframe(image, header.stamp.toSec());
    imageframe.pre_integration = tmp_pre_integration;
    all_image_frame.insert(make_pair(header.stamp.toSec(), imageframe));
    tmp_pre_integration = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};

    // Step 2. 初始化
    //; 如果需要在线标定VI的外参，那么进行标定外参
    if(ESTIMATE_EXTRINSIC == 2)
    {
        ROS_INFO("calibrating extrinsic param, rotation movement is needed");
        if (frame_count != 0)
        {
            vector<pair<Vector3d, Vector3d>> corres = f_manager.getCorresponding(frame_count - 1, frame_count);
            Matrix3d calib_ric;
            if (initial_ex_rotation.CalibrationExRotation(corres, pre_integrations[frame_count]->delta_q, calib_ric))
            {
                ROS_WARN("initial extrinsic rotation calib success");
                ROS_WARN_STREAM("initial extrinsic rotation: " << endl << calib_ric);
                ric[0] = calib_ric;
                RIC[0] = calib_ric;
                ESTIMATE_EXTRINSIC = 1;
            }
        }
    }

    if (solver_flag == INITIAL)
    {
        if (frame_count == WINDOW_SIZE)
        {
            bool result = false;
            //; 初始化
            if( ESTIMATE_EXTRINSIC != 2 && (header.stamp.toSec() - initial_timestamp) > 0.1)
            {
                result = initialStructure();
                initial_timestamp = header.stamp.toSec();
            }
            if(result)
            {
                solver_flag = NON_LINEAR;
                //非线性优化求解VIO
                solveOdometry();
                //滑动窗口
                slideWindow();
                // 移动无效地图点
                f_manager.removeFailures();
                ROS_INFO("Initialization finish!");
                // 保留信息
                last_R = Rs[WINDOW_SIZE];
                last_P = Ps[WINDOW_SIZE];
                last_R0 = Rs[0];
                last_P0 = Ps[0];
                
            }
            else
            {
                slideWindow();
            }
        }
        else
            frame_count++;
    }
    else
    {
        TicToc t_solve;
        solveOdometry();
        ROS_DEBUG("solver costs: %fms", t_solve.toc());

        if (failureDetection())
        {
            ROS_WARN("failure detection!");
            failure_occur = 1;
            clearState();
            setParameter();
            ROS_WARN("system reboot!");
            return;
        }

        TicToc t_margin;
        slideWindow();
        f_manager.removeFailures();
        ROS_DEBUG("marginalization costs: %fms", t_margin.toc());
        key_poses.clear();
        for (int i = 0; i <= WINDOW_SIZE; i++)
            key_poses.push_back(Ps[i]);

        last_R = Rs[WINDOW_SIZE];
        last_P = Ps[WINDOW_SIZE];
        last_R0 = Rs[0];
        last_P0 = Ps[0];
    }
}

// 存储星历信息到类成员变量中
void Estimator::inputEphem(EphemBasePtr ephem_ptr)
{
    double toe = time2sec(ephem_ptr->toe); // 将时间转换为秒,such as x.y秒
    // if a new ephemeris comes
    //; sat: satellite number，应该是卫星编号，如果没有这个卫星编号，那么说明观测到了一个新的卫星，就把它加到数组中
    //; toe: 卫星的时间，如果说之前已经观测过这个卫星，这次又观测到，那么时间肯定不同了，因此也属于一个新的观测，加入数组中
    if (sat2time_index.count(ephem_ptr->sat) == 0 || sat2time_index.at(ephem_ptr->sat).count(toe) == 0)
    {
        sat2ephem[ephem_ptr->sat].emplace_back(ephem_ptr);  
        sat2time_index[ephem_ptr->sat].emplace(toe, sat2ephem.at(ephem_ptr->sat).size()-1);
    }
}

// 存储电离层参数到类成员变量中
void Estimator::inputIonoParams(double ts, const std::vector<double> &iono_params)
{
    // 电离层参数必须是8个
    if (iono_params.size() != 8)    
        return;

    // update ionosphere parameters 更新电离层的参数
    latest_gnss_iono_params.clear();
    std::copy(iono_params.begin(), iono_params.end(), std::back_inserter(latest_gnss_iono_params));
}

void Estimator::inputGNSSTimeDiff(const double t_diff)
{
    diff_t_gnss_local = t_diff;
}


/**
 * @brief 输入当前帧图像匹配的gnss观测信息，进行处理后放到estimator的类成员变量中
 *    注意这里面会根据一些规则对接受到的卫星观测信息进行过滤，只会使用那种满足要求（比如比较稳定）的观测
 * 
 * @param[in] gnss_meas 
 */
void Estimator::processGNSS(const std::vector<ObsPtr> &gnss_meas)
{
    std::vector<ObsPtr> valid_meas;
    std::vector<EphemBasePtr> valid_ephems;

    // 遍历所有的卫星观测信息
    for (auto obs : gnss_meas)
    {
        // filter according to system
        //; 1.首先过滤观测的卫星系统，必须是四大卫星系统
        uint32_t sys = satsys(obs->sat, NULL);
        if (sys != SYS_GPS && sys != SYS_GLO && sys != SYS_GAL && sys != SYS_BDS)
            continue;

        // if not got cooresponding ephemeris yet
        //; 2.如果收到了观测，但是还没有收到星历中关于这个卫星的信息，那么这个观测也不能用
        if (sat2ephem.count(obs->sat) == 0)
            continue;

        // L1频段的
        //! 疑问：如果收到的卫星的频率vector是空，也不能用。这个没太明白
        if (obs->freqs.empty())    
            continue;       // no valid signal measurement
        //; 3.下面这个操作就是判断接收到的消息是否是L1频段的消息，只有L1频段的消息才使用
        int freq_idx = -1;
        L1_freq(obs, &freq_idx);
        if (freq_idx < 0)   
            continue;              // no L1 observation
        
        //; 4.星历信息要有效，否则不使用
        double obs_time = time2sec(obs->time);
        std::map<double, size_t> time2index = sat2time_index.at(obs->sat);
        double ephem_time = EPH_VALID_SECONDS;
        size_t ephem_index = -1;
        for (auto ti : time2index)
        {
            if (std::abs(ti.first - obs_time) < ephem_time)
            {
                ephem_time = std::abs(ti.first - obs_time);
                ephem_index = ti.second;
            }
        }
        if (ephem_time >= EPH_VALID_SECONDS)
        {
            cerr << "ephemeris not valid anymore\n";
            continue;
        }
        const EphemBasePtr &best_ephem = sat2ephem.at(obs->sat).at(ephem_index);

        // filter by tracking status
        //; 5.这个卫星的跟踪次数要足够，这样才比较稳定，否则不使用
        LOG_IF(FATAL, freq_idx < 0) << "No L1 observation found.\n";
        if (obs->psr_std[freq_idx]  > GNSS_PSR_STD_THRES ||
            obs->dopp_std[freq_idx] > GNSS_DOPP_STD_THRES)
        {
            sat_track_status[obs->sat] = 0;
            continue;
        }
        else
        {
            if (sat_track_status.count(obs->sat) == 0)
                sat_track_status[obs->sat] = 0;
            ++ sat_track_status[obs->sat];
        }
        if (sat_track_status[obs->sat] < GNSS_TRACK_NUM_THRES)
            continue;           // not being tracked for enough epochs

        // filter by elevation angle
        //; 6.卫星的仰角越小，说明是低空卫星，会受到更大的电离层延时、多路径效应的影响，也不使用
        if (gnss_ready)
        {
            Eigen::Vector3d sat_ecef;
            if (sys == SYS_GLO)
                sat_ecef = geph2pos(obs->time, std::dynamic_pointer_cast<GloEphem>(best_ephem), NULL);
            else
                sat_ecef = eph2pos(obs->time, std::dynamic_pointer_cast<Ephem>(best_ephem), NULL);
            double azel[2] = {0, M_PI/2.0};
            sat_azel(ecef_pos, sat_ecef, azel);
            if (azel[1] < GNSS_ELEVATION_THRES*M_PI/180.0)
                continue;
        }

        //; 经过上述的判断之后，把有效的观测和星历数据都存储起来
        valid_meas.push_back(obs);
        valid_ephems.push_back(best_ephem);
    }
    
    // 把处理后的观测信息和星历信息存到类成员变量中
    gnss_meas_buf[frame_count] = valid_meas;
    gnss_ephem_buf[frame_count] = valid_ephems;
}
/**
 * @brief 使用sfm算法处理图像，并将其与IMU对齐，与VI初始化相同
 * 
 * @return true 
 * @return false 
 */
bool Estimator::initialStructure()
{
    TicToc t_sfm;
    //check imu observibility，计算IMU预积分的加速度方差，以判断IMU是否活跃
    {
        // 计算加速度均值
        map<double, ImageFrame>::iterator frame_it;
        Vector3d sum_g;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            sum_g += tmp_g;
        }
        Vector3d aver_g;
        aver_g = sum_g * 1.0 / ((int)all_image_frame.size() - 1);
        // 计算加速度方差
        double var = 0;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
            //cout << "frame g " << tmp_g.transpose() << endl;
        }
        var = sqrt(var / ((int)all_image_frame.size() - 1));
        //ROS_WARN("IMU variation %f!", var);
        if(var < 0.25)
        {
            ROS_INFO("IMU excitation not enough!");
            //return false;
        }
    }
    // global sfm
    Quaterniond Q[frame_count + 1];
    Vector3d T[frame_count + 1];
    map<int, Vector3d> sfm_tracked_points;
    vector<SFMFeature> sfm_f;
    // 将特征对应的所有图像帧上的二维点存储在sfm_f中
    for (auto &it_per_id : f_manager.feature)
    {
        int imu_j = it_per_id.start_frame - 1; //start_frame代表观测到该特征的第一帧
        SFMFeature tmp_feature;
        tmp_feature.state = false;
        tmp_feature.id = it_per_id.feature_id;
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            Vector3d pts_j = it_per_frame.point;
            tmp_feature.observation.push_back(make_pair(imu_j, Eigen::Vector2d{pts_j.x(), pts_j.y()}));
        }
        sfm_f.push_back(tmp_feature);
    } 
    Matrix3d relative_R;
    Vector3d relative_T;
    int l;
    // 找出与最新帧相关性强且有足够视差的帧l
    if (!relativePose(relative_R, relative_T, l))
    {
        ROS_INFO("Not enough features or parallax; Move device around");
        return false;
    }
    GlobalSFM sfm;
    // sfm算法，三角化了特征点，并建立了以帧l为参考的位姿Q,T
    if(!sfm.construct(frame_count + 1, Q, T, l,
              relative_R, relative_T,
              sfm_f, sfm_tracked_points))
    {
        ROS_DEBUG("global SFM failed!");
        marginalization_flag = MARGIN_OLD;
        return false;
    }

    //solve pnp for all frame
    map<double, ImageFrame>::iterator frame_it;
    map<int, Vector3d>::iterator it;
    frame_it = all_image_frame.begin( );
    for (int i = 0; frame_it != all_image_frame.end( ); frame_it++)
    {
        // provide initial guess
        cv::Mat r, rvec, t, D, tmp_r;
         // 对于关键帧，直接赋值，不做PnP
        if((frame_it->first) == Headers[i].stamp.toSec())
        {
            frame_it->second.is_key_frame = true;
            frame_it->second.R = Q[i].toRotationMatrix() * RIC[0].transpose();
            frame_it->second.T = T[i];
            i++;
            continue;
        }
        if((frame_it->first) > Headers[i].stamp.toSec())
        {
            i++;
        }
        Matrix3d R_inital = (Q[i].inverse()).toRotationMatrix();
        Vector3d P_inital = - R_inital * T[i];
        cv::eigen2cv(R_inital, tmp_r);
        cv::Rodrigues(tmp_r, rvec); //罗德里格斯公式,将旋转矩阵转化为旋转向量，或将旋转向量转化为旋转矩阵
        cv::eigen2cv(P_inital, t);

        frame_it->second.is_key_frame = false;
        vector<cv::Point3f> pts_3_vector;
        vector<cv::Point2f> pts_2_vector;
        // 对该帧中所有的点建立对应的3d-2d匹配
        for (auto &id_pts : frame_it->second.points)
        {
            int feature_id = id_pts.first;
            for (auto &i_p : id_pts.second)
            {
                it = sfm_tracked_points.find(feature_id); // 该代码与最近的循环无关，提到上一层循环中应是等效的
                if(it != sfm_tracked_points.end())
                {
                    Vector3d world_pts = it->second;
                    cv::Point3f pts_3(world_pts(0), world_pts(1), world_pts(2));
                    pts_3_vector.push_back(pts_3);
                    Vector2d img_pts = i_p.second.head<2>();
                    cv::Point2f pts_2(img_pts(0), img_pts(1));
                    pts_2_vector.push_back(pts_2);
                }
            }
        }
        cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);     
        if(pts_3_vector.size() < 6)
        {
            cout << "pts_3_vector size " << pts_3_vector.size() << endl;
            ROS_DEBUG("Not enough points for solve pnp !");
            return false;
        }
        /* 使用PnP方法
        K为内参，因为用的二维点是相机归一化平面上的点，所以这里的K为单位阵
        D为畸变矩阵
        rvec, t分别是旋转向量和位移
        1 代表使用初始值迭代优化，Parameter used for SOLVEPNP_ITERATIVE. 
        采用的solvePnP方式，CV_P3P、CV_EPNP、CV_ITERATIVE，默认采用CV_ITERATIVE)
        */
        if (! cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1))
        {
            ROS_DEBUG("solve pnp fail!");
            return false;
        }
        
        cv::Rodrigues(rvec, r);
        MatrixXd R_pnp,tmp_R_pnp;
        cv::cv2eigen(r, tmp_R_pnp);
        R_pnp = tmp_R_pnp.transpose();
        MatrixXd T_pnp;
        cv::cv2eigen(t, T_pnp);
        T_pnp = R_pnp * (-T_pnp);
        frame_it->second.R = R_pnp * RIC[0].transpose();
        frame_it->second.T = T_pnp;
    }

    if (!visualInitialAlign())
    {
        ROS_WARN("misalign visual structure with IMU");
        return false;
    }
    return true;
}

/**
 * @brief 将视觉信息与IMU对齐，对齐到第0帧
 * 
 * @return true 
 * @return false 
 */
bool Estimator::visualInitialAlign()
{
    TicToc t_g;
    VectorXd x; // v1, ...,vn,g,s
    //solve scale
    bool result = VisualIMUAlignment(all_image_frame, Bgs, g, x);
    if(!result)
    {
        ROS_DEBUG("solve g failed!");
        return false;
    }

    // change state
    
    for (int i = 0; i <= frame_count; i++)
    {
        Matrix3d Ri = all_image_frame[Headers[i].stamp.toSec()].R;
        Vector3d Pi = all_image_frame[Headers[i].stamp.toSec()].T;
        Ps[i] = Pi;
        Rs[i] = Ri;
        all_image_frame[Headers[i].stamp.toSec()].is_key_frame = true;
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < dep.size(); i++)
        dep[i] = -1; //标记有效的地图点
    f_manager.clearDepth(dep); // 将所有特征点的深度置为-1

    //triangulate on cam pose , no tic
    Vector3d TIC_TMP[NUM_OF_CAM];
    for(int i = 0; i < NUM_OF_CAM; i++)
        TIC_TMP[i].setZero();
    ric[0] = RIC[0];
    f_manager.setRic(ric);
    // 在第0帧中的特征点
    f_manager.triangulate(Ps, &(TIC_TMP[0]), &(RIC[0])); 
    // 滑窗内，重新推算预积分
    double s = (x.tail<1>())(0);
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
    }
    // 开始将所有状态对齐到第0帧
    for (int i = frame_count; i >= 0; i--)
        // twi-tw0, 对齐到第0帧
        Ps[i] = s * Ps[i] - Rs[i] * TIC[0] - (s * Ps[0] - Rs[0] * TIC[0]);
    int kv = -1;
    map<double, ImageFrame>::iterator frame_i;
    for (frame_i = all_image_frame.begin(); frame_i != all_image_frame.end(); frame_i++)
    {
        if(frame_i->second.is_key_frame)
        {
            kv++;
            // 之前获得的速度是IMu坐标系，现在转到world系
            Vs[kv] = frame_i->second.R * x.segment<3>(kv * 3);
        }
    }
    // 为3D点恢复深度
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        it_per_id.estimated_depth *= s;
    }
    // 将所有P，V，Q全部对齐到第0帧，同时与重力方向对齐
    Matrix3d R0 = Utility::g2R(g); // 将枢纽帧l的重力方向转为旋转矩阵，yaw角置零，得到Rwj
    double yaw = Utility::R2ypr(R0 * Rs[0]).x(); // 
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0; //第0帧的yaw角置零
    g = R0 * g;
    //Matrix3d rot_diff = R0 * Rs[0].transpose();
    Matrix3d rot_diff = R0;
    for (int i = 0; i <= frame_count; i++)
    {
        Ps[i] = rot_diff * Ps[i]; //全部对齐到重力下，同时yaw角对齐到第0帧
        Rs[i] = rot_diff * Rs[i];
        Vs[i] = rot_diff * Vs[i];
    }

    ROS_DEBUG_STREAM("g0     " << g.transpose());
    ROS_DEBUG_STREAM("my R0  " << Utility::R2ypr(Rs[0]).transpose());

    return true;
}

/**
 * @brief GNSS和VIO之间的初始化
 * 
 * @return true 
 * @return false 
 */
bool Estimator::GNSSVIAlign()
{
    if (solver_flag == INITIAL)     // visual-inertial not initialized
        return false;
    
    if (gnss_ready)                 // GNSS-VI already initialized
        return true;
    
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
    {
        // 滑窗内gnss信息不能为空，
        //! 注意：这里又出现了对于gnss观测不能少于10个要求，不知道这个10个卫星的要求是这么来的？
        if (gnss_meas_buf[i].empty() || gnss_meas_buf[i].size() < 10) 
            return false;
    }

    // check horizontal velocity excitation，如果水平速度小于0.3m/s,GNSS对齐失败
    //; 这个地方和退化情况一样，速度太低的话多普勒频移都是噪声
    Eigen::Vector2d avg_hor_vel(0.0, 0.0);
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
        avg_hor_vel += Vs[i].head<2>().cwiseAbs();
    avg_hor_vel /= (WINDOW_SIZE+1);
    if (avg_hor_vel.norm() < 0.3)
    {
        std::cerr << "velocity excitation not enough for GNSS-VI alignment.\n";
        return false;
    }

    // 将当前GNSS信息记录起来
    std::vector<std::vector<ObsPtr>> curr_gnss_meas_buf;
    std::vector<std::vector<EphemBasePtr>> curr_gnss_ephem_buf;
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
    {
        curr_gnss_meas_buf.push_back(gnss_meas_buf[i]);
        curr_gnss_ephem_buf.push_back(gnss_ephem_buf[i]);
    }

    //; 传入观测信息、星历信息、最新的电离层参数等信息，构造一个GNSS-VIO的初始化器
    GNSSVIInitializer gnss_vi_initializer(curr_gnss_meas_buf, curr_gnss_ephem_buf, latest_gnss_iono_params);

    // Step 1. SPP单点定位，得到一个锚点的粗略的位置
    // 1. get a rough global location
    Eigen::Matrix<double, 7, 1> rough_xyzt; //; 返回值为锚点的位置xyz, 以及四个卫星系统的时间零偏t1-4
    rough_xyzt.setZero();
    if (!gnss_vi_initializer.coarse_localization(rough_xyzt))
    {
        std::cerr << "Fail to obtain a coarse location.\n";
        return false;
    }

    // Step 2. 使用多普勒频移进行yaw角的对齐
    // 2. perform yaw alignment
    std::vector<Eigen::Vector3d> local_vs;
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
        local_vs.push_back(Vs[i]);
    Eigen::Vector3d rough_anchor_ecef = rough_xyzt.head<3>();
    double aligned_yaw = 0;
    double aligned_rcv_ddt = 0; //时间漂移率
    if (!gnss_vi_initializer.yaw_alignment(local_vs, rough_anchor_ecef, aligned_yaw, aligned_rcv_ddt))
    {
        std::cerr << "Fail to align ENU and local frames.\n";
        return false;
    }
    // std::cout << "aligned_yaw is " << aligned_yaw*180.0/M_PI << '\n';

    // Step 3. 执行锚点的细化
    // 3. perform anchor refinement
    std::vector<Eigen::Vector3d> local_ps;
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
        local_ps.push_back(Ps[i]);
    Eigen::Matrix<double, 7, 1> refined_xyzt;
    refined_xyzt.setZero();
    if (!gnss_vi_initializer.anchor_refinement(local_ps, aligned_yaw, 
        aligned_rcv_ddt, rough_xyzt, refined_xyzt))
    {
        std::cerr << "Fail to refine anchor point.\n";
        return false;
    }
    // std::cout << "refined anchor point is " << std::setprecision(20) 
    //           << refined_xyzt.head<3>().transpose() << '\n';

    // restore GNSS states
    uint32_t one_observed_sys = static_cast<uint32_t>(-1);
    for (uint32_t k = 0; k < 4; ++k)
    {
        if (rough_xyzt(k+3) != 0)
        {
            one_observed_sys = k;
            break;
        }
    }
    // 保存时漂和时偏，对于没有观测到的卫星系统时偏，赋一个观测到的卫星系统时偏
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
    {
        para_rcv_ddt[i] = aligned_rcv_ddt;
        for (uint32_t k = 0; k < 4; ++k)
        {
           
            if (rough_xyzt(k+3) == 0)
                para_rcv_dt[i*4+k] = refined_xyzt(3+one_observed_sys) + aligned_rcv_ddt * i;
            else
                para_rcv_dt[i*4+k] = refined_xyzt(3+k) + aligned_rcv_ddt * i;
        }
    }
    anc_ecef = refined_xyzt.head<3>();
    R_ecef_enu = ecef2rotation(anc_ecef);
    yaw_enu_local = aligned_yaw;

    return true;
}


void Estimator::updateGNSSStatistics()
{
    R_enu_local = Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ());
    enu_pos = R_enu_local * Ps[WINDOW_SIZE];
    enu_vel = R_enu_local * Vs[WINDOW_SIZE];
    enu_ypr = Utility::R2ypr(R_enu_local*Rs[WINDOW_SIZE]);
    ecef_pos = anc_ecef + R_ecef_enu * enu_pos;
}

/**
 * @brief 找出一个参考的帧l，使其与最新帧有足够的相关性和视差，并求解两者的相对位姿
 * 
 */
bool Estimator::relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        vector<pair<Vector3d, Vector3d>> corres;
        corres = f_manager.getCorresponding(i, WINDOW_SIZE);
        // 与最后一帧有20个以上的匹配点
        if (corres.size() > 20)
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();
                sum_parallax = sum_parallax + parallax;

            }
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            if(average_parallax * 460 > 30 && m_estimator.solveRelativeRT(corres, relative_R, relative_T))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * 460, l);
                return true;
            }
        }
    }
    return false;
}

void Estimator::solveOdometry()
{
    if (frame_count < WINDOW_SIZE)
        return;
    if (solver_flag == NON_LINEAR)
    {
        TicToc t_tri;
        // 三角化
        f_manager.triangulate(Ps, tic, ric);
        ROS_DEBUG("triangulation costs %f", t_tri.toc());
        // 后端优化
        optimization();

        if (GNSS_ENABLE)
        {
            //; 如果还没有进行gnss的初始化，则对gnss进行初始化
            if (!gnss_ready)
            {
                gnss_ready = GNSSVIAlign();
            }
            if (gnss_ready)
            {
                updateGNSSStatistics();
            }
        }
    }
}

void Estimator::vector2double()
{
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        para_Pose[i][0] = Ps[i].x();
        para_Pose[i][1] = Ps[i].y();
        para_Pose[i][2] = Ps[i].z();
        Quaterniond q{Rs[i]};
        para_Pose[i][3] = q.x();
        para_Pose[i][4] = q.y();
        para_Pose[i][5] = q.z();
        para_Pose[i][6] = q.w();

        para_SpeedBias[i][0] = Vs[i].x();
        para_SpeedBias[i][1] = Vs[i].y();
        para_SpeedBias[i][2] = Vs[i].z();

        para_SpeedBias[i][3] = Bas[i].x();
        para_SpeedBias[i][4] = Bas[i].y();
        para_SpeedBias[i][5] = Bas[i].z();

        para_SpeedBias[i][6] = Bgs[i].x();
        para_SpeedBias[i][7] = Bgs[i].y();
        para_SpeedBias[i][8] = Bgs[i].z();
    }
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        para_Ex_Pose[i][0] = tic[i].x();
        para_Ex_Pose[i][1] = tic[i].y();
        para_Ex_Pose[i][2] = tic[i].z();
        Quaterniond q{ric[i]};
        para_Ex_Pose[i][3] = q.x();
        para_Ex_Pose[i][4] = q.y();
        para_Ex_Pose[i][5] = q.z();
        para_Ex_Pose[i][6] = q.w();
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        para_Feature[i][0] = dep(i);
    if (ESTIMATE_TD)
        para_Td[0][0] = td;
    
    para_yaw_enu_local[0] = yaw_enu_local;
    for (uint32_t k = 0; k < 3; ++k)
        para_anc_ecef[k] = anc_ecef(k);
}

void Estimator::double2vector()
{
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {

        Rs[i] = Quaterniond(para_Pose[i][6], para_Pose[i][3], 
                            para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
        
        Ps[i] = Vector3d(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);

        Vs[i] = Vector3d(para_SpeedBias[i][0], para_SpeedBias[i][1], para_SpeedBias[i][2]);

        Bas[i] = Vector3d(para_SpeedBias[i][3], para_SpeedBias[i][4], para_SpeedBias[i][5]);

        Bgs[i] = Vector3d(para_SpeedBias[i][6], para_SpeedBias[i][7], para_SpeedBias[i][8]);
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d(para_Ex_Pose[i][0], para_Ex_Pose[i][1], para_Ex_Pose[i][2]);
        ric[i] = Quaterniond(para_Ex_Pose[i][6], para_Ex_Pose[i][3],
                             para_Ex_Pose[i][4], para_Ex_Pose[i][5]).normalized().toRotationMatrix();
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        dep(i) = para_Feature[i][0];
    f_manager.setDepth(dep);
    if (ESTIMATE_TD)
        td = para_Td[0][0];
    
    if (gnss_ready)
    {
        yaw_enu_local = para_yaw_enu_local[0];
        for (uint32_t k = 0; k < 3; ++k)
            anc_ecef(k) = para_anc_ecef[k];
        R_ecef_enu = ecef2rotation(anc_ecef);
    }
}

bool Estimator::failureDetection()
{
    if (f_manager.last_track_num < 2)
    {
        ROS_INFO(" little feature %d", f_manager.last_track_num);
        //return true;
    }
    if (Bas[WINDOW_SIZE].norm() > 2.5)
    {
        ROS_INFO(" big IMU acc bias estimation %f", Bas[WINDOW_SIZE].norm());
        return true;
    }
    if (Bgs[WINDOW_SIZE].norm() > 1.0)
    {
        ROS_INFO(" big IMU gyr bias estimation %f", Bgs[WINDOW_SIZE].norm());
        return true;
    }
    /*
    if (tic(0) > 1)
    {
        ROS_INFO(" big extri param estimation %d", tic(0) > 1);
        return true;
    }
    */
    Vector3d tmp_P = Ps[WINDOW_SIZE];
    if ((tmp_P - last_P).norm() > 5)
    {
        ROS_INFO(" big translation");
        return true;
    }
    if (abs(tmp_P.z() - last_P.z()) > 1)
    {
        ROS_INFO(" big z translation");
        return true; 
    }
    Matrix3d tmp_R = Rs[WINDOW_SIZE];
    Matrix3d delta_R = tmp_R.transpose() * last_R;
    Quaterniond delta_Q(delta_R);
    double delta_angle;
    delta_angle = acos(delta_Q.w()) * 2.0 / 3.14 * 180.0;
    if (delta_angle > 50)
    {
        ROS_INFO(" big delta_angle ");
        //return true;
    }
    return false;
}

void Estimator::optimization()
{
    ceres::Problem problem;
    ceres::LossFunction *loss_function;
    //loss_function = new ceres::HuberLoss(1.0);
    loss_function = new ceres::CauchyLoss(1.0);

    // Step 1. 添加优化变量
    // Step 1.1. 添加滑窗中每一帧的位姿，作为优化变量
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_parameterization);
        problem.AddParameterBlock(para_SpeedBias[i], SIZE_SPEEDBIAS);
    }
    // Step 1.2. 添加VI的外参，作为优化变量
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Ex_Pose[i], SIZE_POSE, local_parameterization);
        // 如果不需要优化外参，设置constant
        if (!ESTIMATE_EXTRINSIC)
        {
            ROS_DEBUG("fix extinsic param");
            problem.SetParameterBlockConstant(para_Ex_Pose[i]);
        }
        else
            ROS_DEBUG("estimate extinsic param");
    }
    // Step 1.3. 添加VI的时间延时，作为优化变量
    if (ESTIMATE_TD)
    {
        problem.AddParameterBlock(para_Td[0], 1);
    }

    // Step 1.4. 如果GNSS已经初始化成功，则添加和GNSS有关的优化变量
    if (gnss_ready)
    {
        // Step 1.4.1. yaw角添加残差块
        problem.AddParameterBlock(para_yaw_enu_local, 1);
        Eigen::Vector2d avg_hor_vel(0.0, 0.0);
        for (uint32_t i = 0; i <= WINDOW_SIZE; ++i)
            avg_hor_vel += Vs[i].head<2>().cwiseAbs();   // cwiseAbs逐个元素取绝对值
        avg_hor_vel /= (WINDOW_SIZE+1);  // 统计滑窗内的速度平均值
        // cerr << "avg_hor_vel is " << avg_vel << endl;

        //; 论文退化情况1：平均速度<0.3, 此时速度低于多普勒频移的噪声水平，yaw角可能会被噪声破坏，因此固定yaw角不优化
        if (avg_hor_vel.norm() < 0.3)
        {
            // std::cerr << "velocity excitation not enough, fix yaw angle.\n";
            problem.SetParameterBlockConstant(para_yaw_enu_local);
        }
        //! 疑问：只要有一帧的卫星观测信息 < 10, 就设置yaw角固定？
        for (uint32_t i = 0; i <= WINDOW_SIZE; ++i)
        {
            if (gnss_meas_buf[i].size() < 10)
                problem.SetParameterBlockConstant(para_yaw_enu_local);
        }

        // Step 1.4.2. ECEF锚点平移添加残差块
        problem.AddParameterBlock(para_anc_ecef, 3);
        // problem.SetParameterBlockConstant(para_anc_ecef);

        // Step 1.4.3. 添加时钟误差和时钟漂移率残差块
        for (uint32_t i = 0; i <= WINDOW_SIZE; ++i)
        {
            //; 注意这个for循环是因为4大卫星系统都各自有一个时间的零偏
            for (uint32_t k = 0; k < 4; ++k)
                problem.AddParameterBlock(para_rcv_dt + i*4 + k, 1);
            problem.AddParameterBlock(para_rcv_ddt + i, 1);
        }
    }

    TicToc t_whole, t_prepare;
    //转成ceres的double类型
    vector2double();

    // Step 2 添加约束边，即各种因子
    // Step 2.1. 首次优化固定锚点变量，如果没有重启，则只会执行一次
    if (first_optimization)
    {
        std::vector<double> anchor_value;
        for (uint32_t k = 0; k < 7; ++k)
            anchor_value.push_back(para_Pose[0][k]);
        PoseAnchorFactor *pose_anchor_factor = new PoseAnchorFactor(anchor_value);
        problem.AddResidualBlock(pose_anchor_factor, NULL, para_Pose[0]);
        first_optimization = false;
    }

    // Step 2.2. 上次边缘化的先验信息
    if (last_marginalization_info)
    {
        // construct new marginlization_factor
        MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
        problem.AddResidualBlock(marginalization_factor, NULL,
                                 last_marginalization_parameter_blocks);
    }

    // Step 2.3. 添加预积分残差
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        int j = i + 1;
        if (pre_integrations[j]->sum_dt > 10.0)
            continue;
        IMUFactor* imu_factor = new IMUFactor(pre_integrations[j]);
        problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);
    }

    // Step 2.4. 添加gnss残差块
    if (gnss_ready)
    {   
        // Step 2.4.1. 添加伪距和多普勒因子
        //; 论文退化情况2：观测卫星数量少于4颗，实际无需处理，因此代码中也没有特殊的处理
        for(int i = 0; i <= WINDOW_SIZE; ++i)
        {
            // cerr << "size of gnss_meas_buf[" << i << "] is " << gnss_meas_buf[i].size() << endl;
            const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[i];
            const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[i];

            // 遍历这个关键帧的所有卫星观测信息
            for (uint32_t j = 0; j < curr_obs.size(); ++j)
            {   
                // sys表示这个观测属于哪个卫星系统的观测
                const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
                // sys_idx表示这个卫星的编号
                //! 疑问：这个编号的意思应该是1234，表示是4大卫星系统中的哪一个？
                const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);

                int lower_idx = -1;
                // 把gnss时间去掉和本地VI时间的偏执，尽量先保持一致
                const double obs_local_ts = time2sec(curr_obs[j]->time) - diff_t_gnss_local;
                // VIO的时间 > 卫星时间，则gnss时间在当前图像帧和上一个图像帧之间
                if (Headers[i].stamp.toSec() > obs_local_ts)
                    lower_idx = (i==0? 0 : i-1);
                // VIO的时间 < 卫星时间，则gnss时间在当前图像帧和下一个图像帧之间
                else
                    lower_idx = (i==WINDOW_SIZE? WINDOW_SIZE-1 : i);
                // 设置当前卫星时间对应于图像时间的前后帧
                const double lower_ts = Headers[lower_idx].stamp.toSec();
                const double upper_ts = Headers[lower_idx+1].stamp.toSec();
                // 时间比率，是指当前卫星时间占图像时间间隔的后半部分的时间
                const double ts_ratio = (upper_ts-obs_local_ts) / (upper_ts-lower_ts);
                //; 构造伪距和多普勒因子，传参（观测信息，星历信息，最新的电离层延迟，卫星时间占后面的时间比例）
                GnssPsrDoppFactor *gnss_factor = new GnssPsrDoppFactor(curr_obs[j], 
                    curr_ephem[j], latest_gnss_iono_params, ts_ratio);
                //; 把多普勒因子加到残差块中，其中变量分别为：
                // 伪距和多普勒因子的类，鲁棒核函数，
                // 前一帧图像的位姿，前一帧图像的速度零偏，后一帧图像的位姿，后一帧图像的速度零偏
                // 当前帧的接收机时间零偏，时间零偏变化率
                // yaw角的偏置，锚点的位置
                problem.AddResidualBlock(gnss_factor, NULL, 
                    para_Pose[lower_idx], para_SpeedBias[lower_idx], para_Pose[lower_idx+1], para_SpeedBias[lower_idx+1],
                    para_rcv_dt+i*4+sys_idx, para_rcv_ddt+i, 
                    para_yaw_enu_local, para_anc_ecef);
            }
        }

        // Step 2.4.2. 添加时间零偏和零偏变化率因子
        //; 论文退化情况3：不管有没有卫星观测，接收机时间零偏和零偏的变化率的约束都成立
        // build relationship between rcv_dt and rcv_ddt
        // 公式(23)，时间零偏靠零偏的变化率来约束
        // 注意这里for k < 4的遍历是因为4大卫星系统都有一个时间零偏
        for (size_t k = 0; k < 4; ++k)
        {
            for (uint32_t i = 0; i < WINDOW_SIZE; ++i)
            {
                // τ_k-1^k
                const double gnss_dt = Headers[i+1].stamp.toSec() - Headers[i].stamp.toSec();
                //; 时间零偏和零偏之间变化的因子
                DtDdtFactor *dt_ddt_factor = new DtDdtFactor(gnss_dt);
                //; 添加零偏和零偏变化率因子到残差块中，其中的参数：
                // 零偏和零偏变化率的类，鲁棒核函数
                // 前一帧的时间零偏，下一帧的时间零偏
                // 前一帧的时间零偏的变化率，下一帧的时间零偏的变化率
                problem.AddResidualBlock(dt_ddt_factor, NULL, 
                    para_rcv_dt+i*4+k, para_rcv_dt+(i+1)*4+k, 
                    para_rcv_ddt+i, para_rcv_ddt+i+1);
            }
        }

        // add rcv_ddt smooth factor
        // 公式(24)，时间零偏的变化率靠随机游走模型来约束
        for (int i = 0; i < WINDOW_SIZE; ++i)
        {
            //; 时间零偏变化率的随机游走因子
            DdtSmoothFactor *ddt_smooth_factor = new DdtSmoothFactor(GNSS_DDT_WEIGHT);
            //; 把时间零偏变化率的随机游走因子加到残差块中，其中的参数：
            // 时间零偏变化率随机游走因子，鲁棒核函数，前一阵的时间零偏变化率，后一帧的时间零偏变化率
            problem.AddResidualBlock(ddt_smooth_factor, NULL, para_rcv_ddt+i, para_rcv_ddt+i+1);
        }
    }

    int f_m_cnt = 0;
    int feature_index = -1;
    // Step 2.5. 添加视觉残差块
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;

        ++feature_index;

        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
        Vector3d pts_i = it_per_id.feature_per_frame[0].point;

        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            if (imu_i == imu_j)
            {
                continue;
            }
            Vector3d pts_j = it_per_frame.point;
            if (ESTIMATE_TD)
            {
                    ProjectionTdFactor *f_td = new ProjectionTdFactor(pts_i, pts_j,
                        it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                        it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                    problem.AddResidualBlock(f_td, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]);
            }
            else
            {
                ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                problem.AddResidualBlock(f, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]);
            }
            f_m_cnt++;
        }
    }

    ROS_DEBUG("visual measurement count: %d", f_m_cnt);
    ROS_DEBUG("prepare for ceres: %f", t_prepare.toc());

    // Step 3. 执行优化
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = NUM_ITERATIONS;
    //options.use_explicit_schur_complement = true;
    // options.minimizer_progress_to_stdout = true;
    options.use_nonmonotonic_steps = true;
    if (marginalization_flag == MARGIN_OLD)
        options.max_solver_time_in_seconds = SOLVER_TIME * 4.0 / 5.0;
    else
        options.max_solver_time_in_seconds = SOLVER_TIME;
    TicToc t_solver;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // cout << summary.BriefReport() << endl;
    // cout << summary.FullReport() << endl;
    ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
    ROS_DEBUG("solver costs: %f", t_solver.toc());

    while(para_yaw_enu_local[0] > M_PI)   para_yaw_enu_local[0] -= 2.0*M_PI;
    while(para_yaw_enu_local[0] < -M_PI)  para_yaw_enu_local[0] += 2.0*M_PI;
    // std::cout << "yaw is " << para_yaw_enu_local[0]*180/M_PI << std::endl;

    double2vector();

    // Step 4 优化完成，计算边缘化
    //仿照ceres的构造模式进行边缘化操作
    TicToc t_whole_marginalization;
    if (marginalization_flag == MARGIN_OLD)
    {
        MarginalizationInfo *marginalization_info = new MarginalizationInfo();
        //eigen 转成 double 类型
        vector2double();

        //添加先验约束
        if (last_marginalization_info)
        {
            vector<int> drop_set;
            for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
            {
                if (last_marginalization_parameter_blocks[i] == para_Pose[0] ||
                    last_marginalization_parameter_blocks[i] == para_SpeedBias[0])
                    drop_set.push_back(i);
            }
            // construct new marginlization_factor
            MarginalizationFactor *marginalization_factor = new MarginalizationFactor(
                last_marginalization_info);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(
                marginalization_factor, NULL, last_marginalization_parameter_blocks, drop_set);
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }
        else
        {
            std::vector<double> anchor_value;
            for (uint32_t k = 0; k < 7; ++k)
                anchor_value.push_back(para_Pose[0][k]);
            PoseAnchorFactor *pose_anchor_factor = new PoseAnchorFactor(anchor_value);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(pose_anchor_factor, 
                NULL, vector<double *>{para_Pose[0]}, vector<int>{0});
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }

        {
            //如果预积分时间过长 不可信
            if (pre_integrations[1]->sum_dt < 10.0)
            {
                IMUFactor* imu_factor = new IMUFactor(pre_integrations[1]);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(imu_factor, NULL,
                                                                           vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1]},
                                                                           vector<int>{0, 1});//（残差雅可比，核函数，涉及参数快，被边缘化索引）
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
        }

        //; 如果GNSS已经初始化过了，那么这里也要计算GNSS的边缘化结果
        if (gnss_ready)
        {
            for (uint32_t j = 0; j < gnss_meas_buf[0].size(); ++j)
            {
                const uint32_t sys = satsys(gnss_meas_buf[0][j]->sat, NULL);
                const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);

                const double obs_local_ts = time2sec(gnss_meas_buf[0][j]->time) - diff_t_gnss_local;
                const double lower_ts = Headers[0].stamp.toSec();
                const double upper_ts = Headers[1].stamp.toSec();
                const double ts_ratio = (upper_ts-obs_local_ts) / (upper_ts-lower_ts);

                GnssPsrDoppFactor *gnss_factor = new GnssPsrDoppFactor(gnss_meas_buf[0][j],
                    gnss_ephem_buf[0][j], latest_gnss_iono_params, ts_ratio);
                ResidualBlockInfo *psr_dopp_residual_block_info = new ResidualBlockInfo(gnss_factor, NULL,
                    vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1],
                        para_SpeedBias[1],para_rcv_dt+sys_idx, para_rcv_ddt,
                        para_yaw_enu_local, para_anc_ecef},
                    vector<int>{0, 1, 4, 5});
                marginalization_info->addResidualBlockInfo(psr_dopp_residual_block_info);
            }

            const double gnss_dt = Headers[1].stamp.toSec() - Headers[0].stamp.toSec();
            for (size_t k = 0; k < 4; ++k)
            {
                DtDdtFactor *dt_ddt_factor = new DtDdtFactor(gnss_dt);
                ResidualBlockInfo *dt_ddt_residual_block_info = new ResidualBlockInfo(dt_ddt_factor, NULL,
                    vector<double *>{para_rcv_dt+k, para_rcv_dt+4+k, para_rcv_ddt, para_rcv_ddt+1}, 
                    vector<int>{0, 2});
                marginalization_info->addResidualBlockInfo(dt_ddt_residual_block_info);
            }

            // margin rcv_ddt smooth factor
            DdtSmoothFactor *ddt_smooth_factor = new DdtSmoothFactor(GNSS_DDT_WEIGHT);
            ResidualBlockInfo *ddt_smooth_residual_block_info = new ResidualBlockInfo(ddt_smooth_factor, NULL,
                    vector<double *>{para_rcv_ddt, para_rcv_ddt+1}, vector<int>{0});
            marginalization_info->addResidualBlockInfo(ddt_smooth_residual_block_info);
        }

        //遍历视觉重投影约束
        {
            int feature_index = -1;
            //遍历特征点管理器
            for (auto &it_per_id : f_manager.feature)
            {
                it_per_id.used_num = it_per_id.feature_per_frame.size();
                //被两帧以上观察到，起始帧在倒数第二帧之前
                if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
                    continue;

                ++feature_index;
                int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
                //起始帧不等于0，跳过
                if (imu_i != 0)
                    continue;
                //该特征点在起始帧的3d点
                Vector3d pts_i = it_per_id.feature_per_frame[0].point;
                //遍历观测到该点的帧
                for (auto &it_per_frame : it_per_id.feature_per_frame)
                {
                    imu_j++;
                    if (imu_i == imu_j)
                        continue;

                    Vector3d pts_j = it_per_frame.point;
                    if (ESTIMATE_TD)
                    {
                        ProjectionTdFactor *f_td = new ProjectionTdFactor(pts_i, pts_j, 
                            it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                            it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f_td, 
                            loss_function, vector<double *>{para_Pose[imu_i], para_Pose[imu_j], 
                                para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]},
                            vector<int>{0, 3});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                    else
                    {
                        //视觉残差会有核函数
                        ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, 
                            loss_function, vector<double *>{para_Pose[imu_i], para_Pose[imu_j], 
                                para_Ex_Pose[0], para_Feature[feature_index]},
                            vector<int>{0, 3});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                }
            }
        }

        TicToc t_pre_margin;
        marginalization_info->preMarginalize();
        ROS_DEBUG("pre marginalization %f ms", t_pre_margin.toc());

        TicToc t_margin;
        marginalization_info->marginalize();
        ROS_DEBUG("marginalization %f ms", t_margin.toc());

        std::unordered_map<long, double *> addr_shift;
        for (int i = 1; i <= WINDOW_SIZE; i++)
        {
            addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
            addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
            for (uint32_t k = 0; k < 4; ++k)
                addr_shift[reinterpret_cast<long>(para_rcv_dt+i*4+k)] = para_rcv_dt+(i-1)*4+k;
            addr_shift[reinterpret_cast<long>(para_rcv_ddt+i)] = para_rcv_ddt+i-1;
        }
        for (int i = 0; i < NUM_OF_CAM; i++)
            addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
        if (ESTIMATE_TD)
        {
            addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
        }
        addr_shift[reinterpret_cast<long>(para_yaw_enu_local)] = para_yaw_enu_local;
        addr_shift[reinterpret_cast<long>(para_anc_ecef)] = para_anc_ecef;
        vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);

        if (last_marginalization_info)
            delete last_marginalization_info;
        last_marginalization_info = marginalization_info;
        last_marginalization_parameter_blocks = parameter_blocks;

    }
    else
    {   
        //边缘化次新帧
        if (last_marginalization_info &&
            std::count(std::begin(last_marginalization_parameter_blocks), std::end(last_marginalization_parameter_blocks), para_Pose[WINDOW_SIZE - 1]))
        {

            MarginalizationInfo *marginalization_info = new MarginalizationInfo();
            vector2double();
            if (last_marginalization_info)
            {
                vector<int> drop_set;
                for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
                {
                    //上一次边缘化的是关键帧，不会出现预积分约束
                    ROS_ASSERT(last_marginalization_parameter_blocks[i] != para_SpeedBias[WINDOW_SIZE - 1]);
                    //找到与倒数第二帧的位姿参数块
                    if (last_marginalization_parameter_blocks[i] == para_Pose[WINDOW_SIZE - 1])
                        drop_set.push_back(i);
                }
                // 构建新的残差块
                MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                               last_marginalization_parameter_blocks,
                                                                               drop_set);

                marginalization_info->addResidualBlockInfo(residual_block_info);
            }

            TicToc t_pre_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->preMarginalize();
            ROS_DEBUG("end pre marginalization, %f ms", t_pre_margin.toc());

            TicToc t_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->marginalize();
            ROS_DEBUG("end marginalization, %f ms", t_margin.toc());

            std::unordered_map<long, double *> addr_shift;
            for (int i = 0; i <= WINDOW_SIZE; i++)
            {
                if (i == WINDOW_SIZE - 1)
                    continue;
                else if (i == WINDOW_SIZE)
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
                    for (uint32_t k = 0; k < 4; ++k)
                        addr_shift[reinterpret_cast<long>(para_rcv_dt+i*4+k)] = para_rcv_dt+(i-1)*4+k;
                    addr_shift[reinterpret_cast<long>(para_rcv_ddt+i)] = para_rcv_ddt+i-1;
                }
                else
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i];
                    for (uint32_t k = 0; k < 4; ++k)
                        addr_shift[reinterpret_cast<long>(para_rcv_dt+i*4+k)] = para_rcv_dt+i*4+k;
                    addr_shift[reinterpret_cast<long>(para_rcv_ddt+i)] = para_rcv_ddt+i;
                }
            }
            for (int i = 0; i < NUM_OF_CAM; i++)
                addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
            if (ESTIMATE_TD)
            {
                addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
            }
            addr_shift[reinterpret_cast<long>(para_yaw_enu_local)] = para_yaw_enu_local;
            addr_shift[reinterpret_cast<long>(para_anc_ecef)] = para_anc_ecef;
            vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);
            if (last_marginalization_info)
                delete last_marginalization_info;
            last_marginalization_info = marginalization_info;
            last_marginalization_parameter_blocks = parameter_blocks;

        }
    }
    ROS_DEBUG("whole marginalization costs: %f", t_whole_marginalization.toc());

    ROS_DEBUG("whole time for ceres: %f", t_whole.toc());
}

void Estimator::slideWindow()
{
    TicToc t_margin;
    if (marginalization_flag == MARGIN_OLD)
    {
        double t_0 = Headers[0].stamp.toSec();
        back_R0 = Rs[0];
        back_P0 = Ps[0];
        //滑窗向前推
        if (frame_count == WINDOW_SIZE)
        {
            //一帧帧交换
            for (int i = 0; i < WINDOW_SIZE; i++)
            {
                Rs[i].swap(Rs[i + 1]);

                std::swap(pre_integrations[i], pre_integrations[i + 1]);

                dt_buf[i].swap(dt_buf[i + 1]);
                linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
                angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);

                Headers[i] = Headers[i + 1];
                Ps[i].swap(Ps[i + 1]);
                Vs[i].swap(Vs[i + 1]);
                Bas[i].swap(Bas[i + 1]);
                Bgs[i].swap(Bgs[i + 1]);

                // GNSS related
                gnss_meas_buf[i].swap(gnss_meas_buf[i+1]);
                gnss_ephem_buf[i].swap(gnss_ephem_buf[i+1]);
                for (uint32_t k = 0; k < 4; ++k)
                    para_rcv_dt[i*4+k] = para_rcv_dt[(i+1)*4+k];
                para_rcv_ddt[i] = para_rcv_ddt[i+1];
            }
            Headers[WINDOW_SIZE] = Headers[WINDOW_SIZE - 1];
            Ps[WINDOW_SIZE] = Ps[WINDOW_SIZE - 1];
            Vs[WINDOW_SIZE] = Vs[WINDOW_SIZE - 1];
            Rs[WINDOW_SIZE] = Rs[WINDOW_SIZE - 1];
            Bas[WINDOW_SIZE] = Bas[WINDOW_SIZE - 1];
            Bgs[WINDOW_SIZE] = Bgs[WINDOW_SIZE - 1];

            // GNSS related
            gnss_meas_buf[WINDOW_SIZE].clear();
            gnss_ephem_buf[WINDOW_SIZE].clear();

            delete pre_integrations[WINDOW_SIZE];
            //等IMU值进入重新预积分
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();

            //初始化相关操作
            if (true || solver_flag == INITIAL)
            {
                map<double, ImageFrame>::iterator it_0;
                //all_image_frame 和初始化相关
                it_0 = all_image_frame.find(t_0);
                delete it_0->second.pre_integration;
                it_0->second.pre_integration = nullptr;
                for (map<double, ImageFrame>::iterator it = all_image_frame.begin(); it != it_0; ++it)
                {
                    if (it->second.pre_integration)
                        delete it->second.pre_integration;
                    it->second.pre_integration = NULL;
                }

                all_image_frame.erase(all_image_frame.begin(), it_0);
                all_image_frame.erase(t_0);

            }
            slideWindowOld();
        }
    }
    else
    {
        if (frame_count == WINDOW_SIZE)
        {
            //合并预计分的值
            for (unsigned int i = 0; i < dt_buf[frame_count].size(); i++)
            {
                double tmp_dt = dt_buf[frame_count][i];
                Vector3d tmp_linear_acceleration = linear_acceleration_buf[frame_count][i];
                Vector3d tmp_angular_velocity = angular_velocity_buf[frame_count][i];

                pre_integrations[frame_count - 1]->push_back(tmp_dt, tmp_linear_acceleration, tmp_angular_velocity);

                dt_buf[frame_count - 1].push_back(tmp_dt);
                linear_acceleration_buf[frame_count - 1].push_back(tmp_linear_acceleration);
                angular_velocity_buf[frame_count - 1].push_back(tmp_angular_velocity);
            }

            Headers[frame_count - 1] = Headers[frame_count];
            Ps[frame_count - 1] = Ps[frame_count];
            Vs[frame_count - 1] = Vs[frame_count];
            Rs[frame_count - 1] = Rs[frame_count];
            Bas[frame_count - 1] = Bas[frame_count];
            Bgs[frame_count - 1] = Bgs[frame_count];

            // GNSS related
            gnss_meas_buf[frame_count-1] = gnss_meas_buf[frame_count];
            gnss_ephem_buf[frame_count-1] = gnss_ephem_buf[frame_count];
            for (uint32_t k = 0; k < 4; ++k)
                para_rcv_dt[(frame_count-1)*4+k] = para_rcv_dt[frame_count*4+k];
            para_rcv_ddt[frame_count-1] = para_rcv_ddt[frame_count];
            gnss_meas_buf[frame_count].clear();
            gnss_ephem_buf[frame_count].clear();

            delete pre_integrations[WINDOW_SIZE];
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();

            slideWindowNew();
        }
    }
}

// real marginalization is removed in solve_ceres()
void Estimator::slideWindowNew()
{
    sum_of_front++;
    f_manager.removeFront(frame_count);
}

// real marginalization is removed in solve_ceres()
//地图点逆深度保存在第一个观察到的帧，如果边缘化掉最后一点，则赋给下一帧
void Estimator::slideWindowOld()
{
    sum_of_back++;
    //更新特征点逆深度的起始帧
    bool shift_depth = solver_flag == NON_LINEAR ? true : false;
    if (shift_depth)
    {
        Matrix3d R0, R1;
        Vector3d P0, P1;
        //back_R0 被移除的帧的位姿
        R0 = back_R0 * ric[0];   //被移除相机姿态
        R1 = Rs[0] * ric[0];     //当前最老的相机姿态
        P0 = back_P0 + back_R0 * tic[0];        //被移除的相机的位置
        P1 = Ps[0] + Rs[0] * tic[0];        //当前最老的相机位置
        f_manager.removeBackShiftDepth(R0, P0, R1, P1);
    }
    else
        f_manager.removeBack();
}
