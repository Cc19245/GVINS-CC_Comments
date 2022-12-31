#pragma once

#include <typeinfo>
#include "parameters.h"
#include "feature_manager.h"
#include "utility/utility.h"
#include "utility/tic_toc.h"
#include "initial/solve_5pts.h"
#include "initial/initial_sfm.h"
#include "initial/initial_alignment.h"
#include "initial/initial_ex_rotation.h"
#include "initial/gnss_vi_initializer.h"
#include <std_msgs/Header.h>

#include <ceres/ceres.h>
#include "factor/imu_factor.h"
#include "factor/pose_local_parameterization.h"
#include "factor/projection_factor.h"
#include "factor/projection_td_factor.h"
#include "factor/marginalization_factor.h"
#include "factor/gnss_psr_dopp_factor.hpp"
#include "factor/gnss_dt_ddt_factor.hpp"
#include "factor/gnss_dt_anchor_factor.hpp"
#include "factor/gnss_ddt_smooth_factor.hpp"
#include "factor/pos_vel_factor.hpp"
#include "factor/pose_anchor_factor.h"

#include <opencv2/core/eigen.hpp>

#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include <gnss_comm/gnss_spp.hpp>

using namespace gnss_comm;

class Estimator
{
  public:
    Estimator();

    void setParameter();

    // interface
    void processIMU(double t, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
    void processGNSS(const std::vector<ObsPtr> &gnss_mea);
    void inputEphem(EphemBasePtr ephem_ptr);
    void inputIonoParams(double ts, const std::vector<double> &iono_params);
    void inputGNSSTimeDiff(const double t_diff);

    void processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, const std_msgs::Header &header);

    // internal
    void clearState();
    bool initialStructure();
    bool visualInitialAlign();
    // GNSS related
    bool GNSSVIAlign();

    void updateGNSSStatistics();

    bool relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l);
    void slideWindow();
    void solveOdometry();
    void slideWindowNew();
    void slideWindowOld();
    void optimization();
    void vector2double();
    void double2vector();
    bool failureDetection();

    enum SolverFlag
    {
        INITIAL,
        NON_LINEAR
    };

    enum MarginalizationFlag
    {
        MARGIN_OLD = 0,
        MARGIN_SECOND_NEW = 1
    };

    SolverFlag solver_flag;
    MarginalizationFlag  marginalization_flag;
    Vector3d g;
    // MatrixXd Ap[2];
    // VectorXd bp[2];

    Matrix3d ric[NUM_OF_CAM];
    Vector3d tic[NUM_OF_CAM];

    Vector3d Ps[(WINDOW_SIZE + 1)];
    Vector3d Vs[(WINDOW_SIZE + 1)];
    Matrix3d Rs[(WINDOW_SIZE + 1)];
    Vector3d Bas[(WINDOW_SIZE + 1)];
    Vector3d Bgs[(WINDOW_SIZE + 1)];
    double td;

    Matrix3d back_R0, last_R, last_R0;
    Vector3d back_P0, last_P, last_P0;
    std_msgs::Header Headers[(WINDOW_SIZE + 1)];

    IntegrationBase *pre_integrations[(WINDOW_SIZE + 1)];
    Vector3d acc_0, gyr_0;

    vector<double> dt_buf[(WINDOW_SIZE + 1)];
    vector<Vector3d> linear_acceleration_buf[(WINDOW_SIZE + 1)];
    vector<Vector3d> angular_velocity_buf[(WINDOW_SIZE + 1)];

    // GNSS related
    bool gnss_ready;   // gnss是否初始化
    Eigen::Vector3d anc_ecef;
    Eigen::Matrix3d R_ecef_enu;
    double yaw_enu_local;
    std::vector<ObsPtr> gnss_meas_buf[(WINDOW_SIZE+1)];   // 过滤之后的观测信息
    std::vector<EphemBasePtr> gnss_ephem_buf[(WINDOW_SIZE+1)];  // 过滤之后的星历信息
    std::vector<double> latest_gnss_iono_params;   // 电离层参数，必须是8个数

    //; 注意下面存储的键都是以卫星编号作为索引的，也就是同一个卫星的所有的观测信息，都是存在一个数组中的
    std::map<uint32_t, std::vector<EphemBasePtr>> sat2ephem;  // 《卫星编号，<星历信息>》
    std::map<uint32_t, std::map<double, size_t>> sat2time_index;  // 《卫星编号，<卫星时间，存储的这个卫星的所有观测的个数》
    std::map<uint32_t, uint32_t> sat_track_status;

    double para_anc_ecef[3];   // ECEF锚点的平移
    double para_yaw_enu_local[1];  // yaw角度的偏置
    double para_rcv_dt[(WINDOW_SIZE+1)*4]; // 接收机的时间零偏，4大卫星系统有4个不同的时间零偏
    double para_rcv_ddt[WINDOW_SIZE+1]; // 接受机的时间零偏的变化率，只和接收机有关，因此4大卫星系统共用一个时间零偏变化率
    // GNSS statistics
    double diff_t_gnss_local;   // PPS触发时间和实际的VI传感器被触发时间之间的差值
    Eigen::Matrix3d R_enu_local;
    Eigen::Vector3d ecef_pos, enu_pos, enu_vel, enu_ypr;

    int frame_count;
    int sum_of_outlier, sum_of_back, sum_of_front, sum_of_invalid;

    FeatureManager f_manager;
    MotionEstimator m_estimator;
    InitialEXRotation initial_ex_rotation;

    bool first_imu;
    bool is_valid, is_key;
    bool failure_occur;

    vector<Vector3d> point_cloud;
    vector<Vector3d> margin_cloud;
    vector<Vector3d> key_poses;
    double initial_timestamp;

    double para_Pose[WINDOW_SIZE + 1][SIZE_POSE];
    double para_SpeedBias[WINDOW_SIZE + 1][SIZE_SPEEDBIAS];
    double para_Feature[NUM_OF_F][SIZE_FEATURE];
    double para_Ex_Pose[NUM_OF_CAM][SIZE_POSE];
    double para_Td[1][1];

    MarginalizationInfo *last_marginalization_info;
    vector<double *> last_marginalization_parameter_blocks;

    map<double, ImageFrame> all_image_frame;
    IntegrationBase *tmp_pre_integration;

    bool first_optimization;
};
