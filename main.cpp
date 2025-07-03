//#include "ode45.hpp"
//#include <iostream>
//#include <vector>
//
//// 定义ODE函数(例如简谐振动系统)
//std::vector<double> harmonic_oscillator(double x, const std::vector<double>& y) {
//    return { -y[0] + x + 1 }; // dy1/dt = y2, dy2/dt = -y1
//}
//
//int main() {
//    // 设置求解区间和初始条件
//    std::vector<double> tspan = { 0.0, 1.0 };
//    std::vector<double> y0 = { 1.0 }; // 初始状态
//
//    // 配置求解器选项
//    ode45::Options options;
//    options.rtol = 1e-6;
//    options.atol = 1e-8;
//	options.fixed_step = true; // 使用定步长
//	options.initial_step = 0.1; // 初始步长
//
//    // 求解ODE
//    auto result = ode45::solve(harmonic_oscillator, tspan, y0, options);
//
//    // 输出结果
//    for (size_t i = 0; i < result.t.size(); ++i) {
//		std::cout << "t: " << result.t[i] << ", y: " << result.y[i][0]<<std::endl;
//    }
//
//    return 0;
//}