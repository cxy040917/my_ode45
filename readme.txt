V1.0

ode45命名空间提供了Runge-Kutta(4,5)方法的实现，用于求解常微分方程(ODE)。这是一个变步长求解器，能够自动调整步长以保持计算精度，同时保证计算效率。

ode45::solve函数
ODEResult solve(
    const ODEFunction& func,
    const std::vector<double>& tspan,
    const std::vector<double>& y0,
    const Options& options)
参数说明:
func: 需要求解的ODE函数，签名为std::vector<double>(double t, const std::vector<double>& y)
tspan: 时间区间[t0, tf]，表示求解的时间范围
y0: 初始状态向量
options: 求解器配置选项(参见下面的Options结构体)


返回值:返回一个ODEResult结构体，包含:
t: 时间点向量
y: 每个时间点对应的状态向量
events: 检测到的事件列表(如果提供了事件函数)

可能抛出的异常:
std::invalid_argument: 输入参数无效时抛出
std::runtime_error: 步长过小时抛出


关键数据结构

ODEResult结构体
struct ODEResult {
    std::vector<double> t;            // 时间点集合
    std::vector<std::vector<double>> y; // 每个时间点的状态
    std::vector<EventInfo> events;     // 检测到的事件
};

Options配置选项
struct Options {
    double rtol = 1e-6;               // 相对容差
    double atol = 1e-8;               // 绝对容差
    double initial_step = -1;         // 初始步长(负值表示自动计算)
    double max_step = std::numeric_limits<double>::infinity(); // 最大步长
    
    // 可选的回调函数，每个成功步长后调用
    std::function<bool(double t, const std::vector<double>& y, 
                      const std::vector<double>& yp)> output_fcn = nullptr;
    
    // 可选的事件函数
    EventFunction event_fcn = nullptr;
    std::vector<int> event_directions; // 事件检测方向(-1,0,1)
    std::vector<bool> event_terminal;  // 是否终止事件
};

EventInfo事件信息
struct EventInfo {
    size_t event_index;               // 事件索引
    double event_time;                // 事件发生时间
    std::vector<double> event_y;      // 事件发生时状态
    bool is_terminal;                 // 是否终止求解
    int direction;                   // 检测到的事件方向
};


使用示例

#include "ode45.hpp"
#include <iostream>
#include <vector>

// 定义ODE函数(例如简谐振动系统)
std::vector<double> harmonic_oscillator(double t, const std::vector<double>& y) {
    return {y[1], -y[0]}; // dy1/dt = y2, dy2/dt = -y1
}

int main() {
    // 设置求解区间和初始条件
    std::vector<double> tspan = {0.0, 10.0};
    std::vector<double> y0 = {1.0, 0.0}; // 初始状态
    
    // 配置求解器选项
    ode45::Options options;
    options.rtol = 1e-6;
    options.atol = 1e-8;
    
    // 求解ODE
    auto result = ode45::solve(harmonic_oscillator, tspan, y0, options);
    
    // 输出结果
    for (size_t i = 0; i < result.t.size(); ++i) {
        std::cout << "t = " << result.t[i] << ", y = [";
        for (auto val : result.y[i]) {
            std::cout << val << " ";
        }
        std::cout << "]\n";
    }
    
    return 0;
}


高级功能

事件检测
可以在求解过程中检测特定事件:
// 定义事件函数(检测y[0]过零点)
auto event_fcn = [](double t, const std::vector<double>& y) {
    return std::vector<double>{y[0]}; // 返回需要检测过零点的量
};

ode45::Options options;
options.event_fcn = event_fcn;
options.event_directions = {0}; // 检测所有方向的过零
options.event_terminal = {false}; // 检测到事件不终止求解

auto result = ode45::solve(harmonic_oscillator, tspan, y0, options);

// 输出检测到的事件
for (const auto& event : result.events) {
    std::cout << "在 t = " << event.event_time << " 检测到事件\n";
}

进度监控
通过输出函数监控求解进度:
auto output_fcn = [](double t, const std::vector<double>& y, 
                    const std::vector<double>& yp) {
    std::cout << "已完成 t = " << t << " 处的计算\n";
    return true; // 返回false可提前终止求解
};

options.output_fcn = output_fcn;


实现细节
求解器采用以下算法流程:
初始化阶段: 检查输入参数并计算初始步长

主循环:
执行Runge-Kutta(4,5)步进
使用嵌入式方法估计误差
根据误差估计调整步长
检查是否发生事件(如果配置了事件函数)
存储结果并调用输出函数

终止条件: 达到结束时间或检测到终止事件
误差控制通过比较4阶和5阶解的差异来实现步长调整。


误差控制机制
求解器通过以下方式保证精度:
使用相对容差(rtol)和绝对容差(atol)
基于误差估计自动调整步长
拒绝误差超过容差的步长

限制说明
当前仅实现Dormand-Prince(4,5)方法对
不支持刚性问题(刚性问题应考虑使用隐式方法)
不支持质量矩阵和雅可比矩阵



V2.0

添加了定步长选项
用法：
	options.fixed_step = true; // 使用定步长
	options.initial_step = 0.1; // 设定步长






