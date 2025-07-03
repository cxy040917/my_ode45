#include "ode45.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

// 数学常量
const double PI = 3.141592653589793;
const double DEG_TO_RAD = PI / 180.0;
const double RAD_TO_DEG = 180.0 / PI;

// 导弹参数
struct MissileParams {
    double m0 = 320.0;     // 初始质量 (kg)
    double Jz = 350.0;     // 转动惯量 (kg.m²)
    double P = 2000.0;     // 发动机推力 (N)
    double ms = 0.46;      // 质量秒消耗量 (kg/s)
    double Sref = 0.45;    // 参考面积 (m²)
    double Lref = 2.5;     // 参考长度 (m)
    double g0 = 9.8;       // 重力加速度 (m/s²)
    double rou0 = 1.2495;  // 初始大气密度 (kg/m³)
    double T0 = 288.15;    // 初始温度 (K)
};

// 目标参数 - 简化：固定目标位置
struct TargetParams {
    std::vector<double> position0 = { 2000.0, 0.0, 2000.0 }; // 初始位置 [x, y, z] (m)
    std::vector<double> velocity = { 15.0, 10.0, 0.0 };     // 速度为0，固定目标
};

// PID 控制器参数 - 增大增益
struct PIDParams {
    double K_P_theta = 20.0;   // 增大比例增益
    double K_D_theta = 1.0;   // 增大微分增益
    double K_I_theta = 0.1;   // 增大积分增益
    double K_P_varphi = 20.0;
    double K_D_varphi = 1.0;
    double K_I_varphi = 0.1;
};

// 导引头参数
struct SeekerParams {
    double FOV = 60.0;               // 视场角 (度)
    double target_actual_size = 2.0; // 目标实际尺寸 (m)
    double reference_distance = 1000.0; // 基准距离 (m)
    double reference_area = 100.0;   // 基准面积
};

// 仿真状态
struct SimulationState {
    std::vector<double> error_theta;
    std::vector<double> error_varphi;
    std::vector<double> delta_z;
    std::vector<double> delta_y;
    std::vector<double> alpha;
    std::vector<double> bleta;
    std::vector<double> q1;
    std::vector<double> q2;
    std::vector<double> target_x;
    std::vector<double> target_y;
    std::vector<double> target_z;
    double sum1 = 0.0;
    double sum2 = 0.0;
};

// 目标运动函数
std::vector<double> target_move(const std::vector<double>& position0,
    const std::vector<double>& velocity,
    double t) {
    return {
        position0[0] + velocity[0] * t,
        position0[1] + velocity[1] * t,
        position0[2] + velocity[2] * t
    };
}

// 视线角计算
std::vector<double> calculate_light_angle(double x, double y, double z,
    double x_t, double y_t, double z_t) {
    double dx = x_t - x;
    double dy = y_t - y;
    double dz = z_t - z;

    double distance_xy = sqrt(dx * dx + dz * dz);
    double q1 = atan2(dy, distance_xy); // 俯仰角
    double q2 = atan2(dz, dx);          // 偏航角

    return { q1, q2 };
}

// 导引头数据获取
struct SeekerData {
    std::vector<double> target_center;
    double target_area;
    std::vector<double> image_center;
    double image_area;
    bool is_target_detected;
};

// 更新导引头函数签名，添加时间参数用于调试
SeekerData getSeekerData(const std::vector<double>& missile_pos,
    const std::vector<double>& target_pos,
    double seeker_elevation_deg,
    double seeker_azimuth_deg,
    double t) {  // 添加时间参数
    SeekerData data;
    double FOV = 60.0; // 视场角 (度)

    // 计算相对位置
    std::vector<double> target_relative(3);
    for (int i = 0; i < 3; i++) {
        target_relative[i] = target_pos[i] - missile_pos[i];
    }

    double distance = sqrt(target_relative[0] * target_relative[0] +
        target_relative[1] * target_relative[1] +
        target_relative[2] * target_relative[2]);

    // 构建变换矩阵
    double az_rad = seeker_azimuth_deg * DEG_TO_RAD;
    double el_rad = seeker_elevation_deg * DEG_TO_RAD;

    // 方位旋转矩阵
    double R_az[3][3] = {
        {cos(az_rad), 0, sin(az_rad)},
        {0, 1, 0},
        {-sin(az_rad), 0, cos(az_rad)}
    };

    // 俯仰旋转矩阵
    double R_el[3][3] = {
        {cos(el_rad), -sin(el_rad), 0},
        {sin(el_rad), cos(el_rad), 0},
        {0, 0, 1}
    // 修正方位旋转矩阵（可能方向定义错误）
    //double R_az[3][3] = {
    //    {cos(az_rad), 0, sin(az_rad)},  // 修改第三列符号
    //    {0, 1, 0},
    //    {-sin(az_rad), 0, cos(az_rad)} // 修改第一行第三列和第三行第一列
    //};

    //// 修正俯仰旋转矩阵
    //double R_el[3][3] = {
    //    {cos(el_rad), 0, sin(el_rad)},  // 调整元素位置
    //    {0, 1, 0},
    //    {-sin(el_rad), 0, cos(el_rad)}
    };

    // 组合旋转矩阵 R_seeker = R_el * R_az
    double R_seeker[3][3] = { {0} };
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                R_seeker[i][j] += R_el[i][k] * R_az[k][j];
                //R_seeker[i][j] += R_az[i][k] * R_el[k][j];
            }
        }
    }

    // 坐标变换
    std::vector<double> target_seeker(3, 0.0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            target_seeker[i] += R_seeker[i][j] * target_relative[j];
        }
    }

    // 图像中心
    data.image_center = { 0.0, 0.0 };
    data.image_area = FOV * FOV;

    // 初始化
    data.target_center = { 0.0, 0.0 };
    data.target_area = 0.0;
    data.is_target_detected = false;

    // 目标在导引头前方
    if (target_seeker[0] > 0) {
        double horizontal_angle = atan2(target_seeker[2], target_seeker[0]) * RAD_TO_DEG;
        double vertical_angle = atan2(target_seeker[1],
            sqrt(target_seeker[0] * target_seeker[0] +
                target_seeker[2] * target_seeker[2])) * RAD_TO_DEG;
        //double vertical_angle = atan2(target_seeker[1], target_seeker[0]) * RAD_TO_DEG;

        // 检查是否在视场内
        if (fabs(horizontal_angle) <= FOV / 2 && fabs(vertical_angle) <= FOV / 2) {
            data.is_target_detected = true;
            data.target_center = { horizontal_angle, vertical_angle };

            // 计算目标面积 (与距离平方成反比)
            double target_actual_size = 2.0;
            double reference_distance = 1000.0;
            double reference_area = 100.0;

            double distance_ratio = reference_distance / std::max(distance, 1.0);
            data.target_area = reference_area * distance_ratio * distance_ratio;

            // 限制面积范围
            double min_dot_size = 5.0;
            double max_dot_size = 5000.0;
            data.target_area = std::max(min_dot_size, std::min(data.target_area, max_dot_size));

        }
        if (target_seeker[0] <= 0) {
            std::cout << "[GUIDANCE] t=" << t << " 目标不在导引头前方! X_seeker: "
                << target_seeker[0] << "\n";
        }
        else if (!(fabs(horizontal_angle) <= FOV / 2 || fabs(vertical_angle) <= FOV / 2)) {
            std::cout << "[GUIDANCE] t=" << t << " 目标在导引头前方但不在视场内! "
                << "水平角: " << horizontal_angle << "° (FOV: ±" << FOV / 2 << "°) "
                << "垂直角: " << vertical_angle << "° (FOV: ±" << FOV / 2 << "°)\n";
        }
    }

    return data;
}
//SeekerData getSeekerData(const std::vector<double>& missile_pos,
//    const std::vector<double>& target_pos,
//    double missile_pitch_deg,
//    double missile_yaw_deg,
//    double t) {
//    SeekerData data;
//    const double FOV = 60.0; // 视场角 (度)
//
//    // 计算相对位置向量 (导弹->目标)
//    const double dx = target_pos[0] - missile_pos[0];
//    const double dy = target_pos[1] - missile_pos[1];
//    const double dz = target_pos[2] - missile_pos[2];
//    const double distance = sqrt(dx * dx + dy * dy + dz * dz);
//
//    // 导弹姿态角度转弧度
//    const double theta = missile_pitch_deg * DEG_TO_RAD; // 俯仰角
//    const double psi = missile_yaw_deg * DEG_TO_RAD;     // 偏航角
//
//    // 第一步：将相对位置向量从地面坐标系转换到导弹坐标系
//    // 使用正确的坐标系转换公式
//    double rel_missile_x = dx * cos(psi) * cos(theta) +
//        dy * sin(theta) +
//        dz * sin(psi) * cos(theta);
//
//    double rel_missile_y = -dx * cos(psi) * sin(theta) +
//        dy * cos(theta) -
//        dz * sin(psi) * sin(theta);
//
//    double rel_missile_z = -dx * sin(psi) +
//        dz * cos(psi);
//
//    // 计算导引头坐标系中的角度
//    double horizontal_angle = atan2(rel_missile_z, rel_missile_x) * RAD_TO_DEG;
//    double vertical_angle = atan2(rel_missile_y, rel_missile_x) * RAD_TO_DEG;
//
//    // 初始化返回数据
//    data.image_center = { 0.0, 0.0 };
//    data.image_area = FOV * FOV;
//    data.target_center = { 0.0, 0.0 };
//    data.target_area = 0.0;
//    data.is_target_detected = false;
//
//    // 动态调整视场角 - 近距离时扩大视场
//    double effective_FOV = FOV;
//    if (distance < 1000.0) effective_FOV = 90.0;
//    if (distance < 500.0) effective_FOV = 120.0;
//    if (distance < 100.0) effective_FOV = 180.0;
//
//    // 检查目标是否在导引头前方（X分量>0）且在视场内
//    if (rel_missile_x > 0.1) {
//        // 检查是否在视场内
//        if (fabs(horizontal_angle) <= effective_FOV / 2 &&
//            fabs(vertical_angle) <= effective_FOV / 2) {
//
//            data.is_target_detected = true;
//            data.target_center = { horizontal_angle, vertical_angle };
//
//            // 计算目标面积（与距离平方成反比）
//            const double reference_distance = 1000.0;
//            const double reference_area = 100.0;
//            const double distance_ratio = reference_distance / std::max(distance, 1.0);
//            data.target_area = reference_area * distance_ratio * distance_ratio;
//
//            // 限制面积范围
//            data.target_area = std::max(5.0, std::min(data.target_area, 5000.0));
//        }
//    }
//
//    // 调试输出 - 当距离小于1000米时打印详细信息
//    if (distance < 1000.0) {
//        std::cout << "[GUIDANCE-DEBUG] t=" << t << " distance=" << distance << "m\n";
//        std::cout << "        missile pos: (" << missile_pos[0] << ", "
//            << missile_pos[1] << ", " << missile_pos[2] << ")\n";
//        std::cout << "        target pos: (" << target_pos[0] << ", "
//            << target_pos[1] << ", " << target_pos[2] << ")\n";
//        std::cout << "        relative vector: (" << dx << ", " << dy << ", " << dz << ")\n";
//        std::cout << "        missile coord: (" << rel_missile_x << ", "
//            << rel_missile_y << ", " << rel_missile_z << ")\n";
//        std::cout << "        angles: horizontal=" << horizontal_angle
//            << "°, vertical=" << vertical_angle << "°\n";
//        std::cout << "        is_target_detected: " << data.is_target_detected << "\n";
//    }
//
//    return data;
//}


// 导弹动力学模型
// 导弹动力学模型
class MissileModel {
public:
    MissileModel(const MissileParams& mp, const TargetParams& tp, const PIDParams& pp, const SeekerParams& sp)
        : missile_params(mp), target_params(tp), pid_params(pp), seeker_params(sp) {
    }

    // ODE 右侧函数
    std::vector<double> operator()(double t, const std::vector<double>& y) {
        // 状态变量
        double v = y[0];        // 速度
        double theta = y[1];    // 弹道倾角 (rad)
        double varphi = y[2];   // 弹道偏角 (rad)
        double x = y[3];        // x 位置
        double y_pos = y[4];    // y 位置 (高度)
        double z = y[5];        // z 位置
        double mass = y[6];     // 质量

        // 目标位置
        auto target_pos = target_move(target_params.position0, target_params.velocity, t);
        target_x.push_back(target_pos[0]);
        target_y.push_back(target_pos[1]);
        target_z.push_back(target_pos[2]);

        // 获取导引头数据
        auto seeker_data = getSeekerData({ x, y_pos, z }, target_pos,
            theta * RAD_TO_DEG, varphi * RAD_TO_DEG, t);  // 添加时间参数用于调试

        double error_theta_i = 0.0;
        double error_varphi_i = 0.0;

        if (seeker_data.is_target_detected) {
            double dis_cal = sqrt(seeker_data.image_area / seeker_data.target_area) * 100;

            // 打印导引头检测信息
            if (debug_counter % 5 == 0) {
                std::cout << "[SUCCESS] t=" << t << " 导引头检测到目标! "
                    << "水平角: " << seeker_data.target_center[0]
                    << "° 垂直角: " << seeker_data.target_center[1] << "°\n";
            }

            // 修正误差计算
            error_theta_i = (seeker_data.target_center[1] - seeker_data.image_center[1]) * DEG_TO_RAD;
            error_varphi_i = (seeker_data.target_center[0] - seeker_data.image_center[0]) * DEG_TO_RAD;

            sum1 += error_theta_i;
            sum2 += error_varphi_i;

            if (dis_cal < 100) {
                error_theta_i = (seeker_data.target_center[1] - seeker_data.image_center[1]) * DEG_TO_RAD * 2;
                error_varphi_i = (seeker_data.target_center[0] - seeker_data.image_center[0]) * DEG_TO_RAD * 2;
            }
        }
        else {
            auto q = calculate_light_angle(x, y_pos, z, target_pos[0], target_pos[1], target_pos[2]);
            q1.push_back(q[0]);
            q2.push_back(q[1]);

            // 打印视线角信息
            if (debug_counter % 5 == 0) {
                std::cout << "[WARNING] t=" << t << " 导引头未检测到目标! "
                    << "视线俯仰角: " << q[0] * RAD_TO_DEG
                    << "° 视线偏航角: " << q[1] * RAD_TO_DEG << "°\n";
            }

            // 修正误差计算
            error_theta_i = q[0] - theta;
            error_varphi_i = q[1] - varphi;
        }

        debug_counter++;

        error_theta.push_back(error_theta_i);
        error_varphi.push_back(error_varphi_i);

        // PID 控制器
        double delta_z_i = 0.0;
        double delta_y_i = 0.0;

        if (error_theta.size() < 2) {
            delta_z_i = pid_params.K_P_theta * error_theta_i;
            delta_y_i = pid_params.K_P_varphi * error_varphi_i;
        }
        else {
            double d_error_theta = error_theta_i - error_theta[error_theta.size() - 2];
            double d_error_varphi = error_varphi_i - error_varphi[error_varphi.size() - 2];

            delta_z_i = pid_params.K_P_theta * error_theta_i +
                pid_params.K_D_theta * d_error_theta +
                pid_params.K_I_theta * sum1;

            delta_y_i = pid_params.K_P_varphi * error_varphi_i +
                pid_params.K_D_varphi * d_error_varphi +
                pid_params.K_I_varphi * sum2;
        }

        // 限制舵偏角
        delta_z_i = std::max(-30.0 * DEG_TO_RAD, std::min(delta_z_i, 30.0 * DEG_TO_RAD));
        delta_y_i = std::max(-30.0 * DEG_TO_RAD, std::min(delta_y_i, 30.0 * DEG_TO_RAD));

        delta_z.push_back(delta_z_i);
        delta_y.push_back(delta_y_i);

        // 瞬时平衡假设
        double alpha_i = 0.24 * delta_z_i;
        double bleta_i = 0.24 * delta_y_i;

        alpha.push_back(alpha_i);
        bleta.push_back(bleta_i);

        double alpha_deg = alpha_i * RAD_TO_DEG;
        double bleta_deg = bleta_i * RAD_TO_DEG;

        // 气动力系数
        double C_X = 0.2 + 0.005 * alpha_deg * alpha_deg;
        double C_Y = 0.25 * alpha_deg + 0.05 * delta_z_i * RAD_TO_DEG;
        double C_Z = -0.25 * bleta_deg - 0.05 * delta_y_i * RAD_TO_DEG;

        // 大气参数
        double T = missile_params.T0 - 0.0065 * y_pos;
        double rou = missile_params.rou0 * pow(T / missile_params.T0, 4.25588);
        double q_dynamic = 0.5 * rou * v * v;

        // 气动力
        double X = C_X * q_dynamic * missile_params.Sref;
        double Y = C_Y * q_dynamic * missile_params.Sref;
        double Z = C_Z * q_dynamic * missile_params.Sref;

        // 推力
        double thrust = (t < 2.0) ? missile_params.P : 0.0;
        double dmass = (t < 2.0) ? -missile_params.ms : 0.0; // 质量变化率 (kg/s)

        // 导弹运动方程组
        double dv = (thrust * cos(alpha_i) * cos(bleta_i) - X - mass * missile_params.g0 * sin(theta)) / mass;
        double dtheta_val = (thrust * sin(alpha_i) + Y - mass * missile_params.g0 * cos(theta)) / (mass * v);
        double dvarphi_val = -(-thrust * cos(alpha_i) * sin(bleta_i) + Z) / (mass * v * cos(theta));
        double dx = v * cos(theta) * cos(varphi);
        double dy = v * sin(theta);
        double dz = v * cos(theta) * sin(varphi); // 已修正z方向

        // 打印关键参数
        if (debug_counter % 5 == 0) {
            auto target_pos = target_move(target_params.position0, target_params.velocity, t);
            double dx_dist = x - target_pos[0];
            double dy_dist = y_pos - target_pos[1];
            double dz_dist = z - target_pos[2];
            double distance = sqrt(dx_dist * dx_dist + dy_dist * dy_dist + dz_dist * dz_dist);

            std::cout << "[DEBUG] t=" << t << " 导弹位置: (" << x << ", " << y_pos << ", " << z << ")"
                << " 目标位置: (" << target_pos[0] << ", " << target_pos[1] << ", " << target_pos[2] << ")"
                << " 距离: " << distance << " m\n";
            std::cout << "       速度: " << v << " m/s"
                << " 俯仰角: " << theta * RAD_TO_DEG << "°"
                << " 偏航角: " << varphi * RAD_TO_DEG << "°"
                << " 质量: " << mass << " kg"
                << " 质量变化率: " << dmass << " kg/s\n";
        }

        return { dv, dtheta_val, dvarphi_val, dx, dy, dz, dmass };
    }

    // 获取仿真状态
    SimulationState getSimulationState() const {
        return {
            error_theta, error_varphi, delta_z, delta_y, alpha, bleta,
            q1, q2, target_x, target_y, target_z, sum1, sum2
        };
    }

private:
    MissileParams missile_params;
    TargetParams target_params;
    PIDParams pid_params;
    SeekerParams seeker_params;

    // 仿真状态变量
    std::vector<double> error_theta;
    std::vector<double> error_varphi;
    std::vector<double> delta_z;
    std::vector<double> delta_y;
    std::vector<double> alpha;
    std::vector<double> bleta;
    std::vector<double> q1;
    std::vector<double> q2;
    std::vector<double> target_x;
    std::vector<double> target_y;
    std::vector<double> target_z;
    double sum1 = 0.0;
    double sum2 = 0.0;
    int debug_counter = 0;
};

// 打印导弹状态数据
void printMissileData(const ode45::ODEResult& result) {
    // 打印表头
    std::cout << "\n导弹状态数据:\n";
    std::cout << std::setw(8) << "时间(s)"
        << std::setw(12) << "速度(m/s)"
        << std::setw(15) << "弹道倾角(°)"
        << std::setw(15) << "弹道偏角(°)"
        << std::setw(12) << "X位置(m)"
        << std::setw(12) << "Y位置(m)"
        << std::setw(12) << "Z位置(m)"
        << std::setw(12) << "质量(kg)"
        << "\n";

    // 打印数据 - 每10步打印一次
    for (size_t i = 0; i < result.t.size(); i += 10) {
        const auto& state = result.y[i];
        std::cout << std::setw(8) << std::setprecision(4) << result.t[i]
            << std::setw(12) << std::setprecision(6) << state[0]
            << std::setw(15) << std::setprecision(6) << state[1] * RAD_TO_DEG
            << std::setw(15) << std::setprecision(6) << state[2] * RAD_TO_DEG
            << std::setw(12) << std::setprecision(2) << state[3]
            << std::setw(12) << std::setprecision(2) << state[4]
            << std::setw(12) << std::setprecision(2) << state[5]
            << std::setw(12) << std::setprecision(2) << state[6]
            << "\n";
    }
}

// 打印附加状态数据
void printAdditionalData(const SimulationState& sim_state) {
    // 打印表头
    std::cout << "\n附加状态数据 (前50步):\n";
    std::cout << std::setw(8) << "Index"
        << std::setw(15) << "俯仰角误差"
        << std::setw(15) << "偏航角误差"
        << std::setw(15) << "俯仰舵偏(°)"
        << std::setw(15) << "偏航舵偏(°)"
        << std::setw(15) << "攻角(°)"
        << std::setw(15) << "侧滑角(°)"
        << "\n";

    // 打印数据 - 仅打印前50行
    size_t print_count = std::min(sim_state.error_theta.size(), static_cast<size_t>(50));
    for (size_t i = 0; i < print_count; i++) {
        std::cout << std::setw(8) << i
            << std::setw(15) << std::setprecision(6) << sim_state.error_theta[i]
            << std::setw(15) << std::setprecision(6) << sim_state.error_varphi[i]
            << std::setw(15) << std::setprecision(6) << sim_state.delta_z[i] * RAD_TO_DEG
            << std::setw(15) << std::setprecision(6) << sim_state.delta_y[i] * RAD_TO_DEG
            << std::setw(15) << std::setprecision(6) << sim_state.alpha[i] * RAD_TO_DEG
            << std::setw(15) << std::setprecision(6) << sim_state.bleta[i] * RAD_TO_DEG
            << "\n";
    }

    if (sim_state.error_theta.size() > 50) {
        std::cout << "..." << (sim_state.error_theta.size() - 50) << " 更多行未显示...\n";
    }
}

int main() {
    // 设置参数
    MissileParams missile_params;
    TargetParams target_params;
    PIDParams pid_params;
    SeekerParams seeker_params;

    // 初始状态
    std::vector<double> initial_conditions = {
        250.0,   // 初始速度 (m/s)
        0.0,     // 初始弹道倾角 (rad)
        0.0,     // 初始弹道偏角 (rad)
        0.0,     // 初始x位置 (m)
        7000.0,  // 初始y位置 (高度) (m)
        0.0,     // 初始z位置 (m)
        missile_params.m0 // 初始质量 (kg)
    };

    // 创建导弹模型
    MissileModel missile_model(missile_params, target_params, pid_params, seeker_params);

    // 配置ODE45选项
    ode45::Options options;
    options.rtol = 1e-6;
    options.atol = 1e-6;
    options.initial_step = 0.01;
    options.max_step = 0.1;

    // 事件函数 - 计算导弹到目标的实际距离
    options.event_fcn = [target_params](double t, const std::vector<double>& y) {
        // 计算目标位置
        auto target_pos = target_move(target_params.position0, target_params.velocity, t);

        // 计算导弹到目标的距离
        double dx = y[3] - target_pos[0];
        double dy = y[4] - target_pos[1];
        double dz = y[5] - target_pos[2];
        double distance = sqrt(dx * dx + dy * dy + dz * dz);

        return std::vector<double>{distance - 10.0}; // 距离小于10m时触发事件
        };

    options.event_directions = { -1 }; // 仅当距离减小时触发
    options.event_terminal = { true }; // 事件触发时终止

    // 设置时间范围
    std::vector<double> tspan = { 0.0, 50.0 };

    // 运行仿真
    std::cout << "开始导弹仿真...\n";
    std::cout << "目标位置: (" << target_params.position0[0] << ", "
        << target_params.position0[1] << ", " << target_params.position0[2] << ")\n";
    std::cout << "目标速度: (" << target_params.velocity[0] << ", "
        << target_params.velocity[1] << ", " << target_params.velocity[2] << ")\n";

    auto result = ode45::solve(missile_model, tspan, initial_conditions, options);
    std::cout << "仿真完成! 共 " << result.t.size() << " 个时间步\n";

    // 打印导弹状态数据
    printMissileData(result);

    // 获取附加状态
    auto sim_state = missile_model.getSimulationState();

    // 打印附加状态数据
    printAdditionalData(sim_state);

    // 打印摘要信息
    if (!result.events.empty()) {
        std::cout << "\n事件信息:\n";
        for (const auto& event : result.events) {
            std::cout << "事件时间: " << event.event_time << " s\n";
            std::cout << "是否终止: " << (event.is_terminal ? "是" : "否") << "\n";

            // 打印事件发生时导弹和目标的位置
            auto target_pos = target_move(target_params.position0, target_params.velocity, event.event_time);
            std::cout << "导弹位置: ("
                << event.event_y[3] << ", "
                << event.event_y[4] << ", "
                << event.event_y[5] << ") m\n";
            std::cout << "目标位置: ("
                << target_pos[0] << ", "
                << target_pos[1] << ", "
                << target_pos[2] << ") m\n";
        }
    }
    else {
        std::cout << "\n未检测到事件\n";

        // 计算最终距离
        if (!result.y.empty()) {
            const auto& final_state = result.y.back();
            auto final_target_pos = target_move(target_params.position0, target_params.velocity, result.t.back());
            double dx = final_state[3] - final_target_pos[0];
            double dy = final_state[4] - final_target_pos[1];
            double dz = final_state[5] - final_target_pos[2];
            double final_distance = sqrt(dx * dx + dy * dy + dz * dz);

            std::cout << "最终导弹-目标距离: " << final_distance << " m\n";
        }
    }

    // 打印最终状态
    if (!result.y.empty()) {
        const auto& final_state = result.y.back();
        auto final_target_pos = target_move(target_params.position0, target_params.velocity, result.t.back());

        std::cout << "\n最终状态:\n";
        std::cout << "时间: " << result.t.back() << " s\n";
        std::cout << "导弹位置: (" << final_state[3] << ", " << final_state[4]
            << ", " << final_state[5] << ") m\n";
        std::cout << "目标位置: (" << final_target_pos[0] << ", " << final_target_pos[1]
            << ", " << final_target_pos[2] << ") m\n";
        std::cout << "速度: " << final_state[0] << " m/s\n";
        std::cout << "质量: " << final_state[6] << " kg\n";
    }

    return 0;
}