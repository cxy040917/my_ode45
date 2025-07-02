#ifndef ODE45_HPP
#define ODE45_HPP

#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <limits>
#include <utility>

namespace ode45 {

    // 定义微分方程类型: dy/dt = f(t, y)
    using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

    // 定义事件函数类型
    using EventFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

    // 定义输出函数类型
    using OutputFcn = std::function<bool(double, const std::vector<double>&, const std::vector<double>&)>;

    // 事件信息结构体
    struct EventInfo {
        size_t event_index;      // 事件索引
        double event_time;       // 事件发生时间
        std::vector<double> event_y; // 事件发生时的状态
        bool is_terminal;        // 是否为终止事件
        int direction;           // 事件方向 (0=所有, 1=正, -1=负)
    };

    // 定义结果结构体
    struct ODEResult {
        std::vector<double> t;
        std::vector<std::vector<double>> y;
        std::vector<EventInfo> events;
    };

    // 默认容差
    constexpr double DEFAULT_RTOL = 1e-6;
    constexpr double DEFAULT_ATOL = 1e-8;

    // 求解选项结构体
    struct Options {
        double rtol = DEFAULT_RTOL;
        double atol = DEFAULT_ATOL;
        double max_step = std::numeric_limits<double>::infinity();
        double initial_step = 0.0;
        EventFunction event_fcn = nullptr;
        OutputFcn output_fcn = nullptr;
        std::vector<int> event_directions; // 每个事件的方向 (0=所有, 1=正, -1=负)
        std::vector<bool> event_terminal; // 每个事件是否为终止事件
    };

    // 主求解函数
    ODEResult solve(
        const ODEFunction& func,
        const std::vector<double>& tspan,
        const std::vector<double>& y0,
        const Options& options = Options());

    // 辅助函数
    namespace internal {
        std::vector<double> rk45_step(
            const ODEFunction& func,
            double t,
            const std::vector<double>& y,
            double h);

        double estimate_error(
            const std::vector<double>& y,
            const std::vector<double>& y_high,
            const std::vector<double>& y_low,
            double rtol,
            double atol);

        double compute_initial_step(
            const ODEFunction& func,
            double t0,
            const std::vector<double>& y0,
            double rtol,
            double atol);

        // @brief 检测事件的条件：
        // 1. 初始状态不触发（|val_old| > tolerance）
        // 2. 值函数符号变化（val_old * val_new <= 0）
        // 3. 方向匹配（direction == 0 或与斜率一致）
        std::vector<EventInfo> check_events(
            const EventFunction& event_fcn,
            const std::vector<int>& directions,
            const std::vector<bool>& terminal,
            double t_old, double t_new,
            const std::vector<double>& y_old, const std::vector<double>& y_new);

        std::pair<double, std::vector<double>> interpolate(
            double t_old, double t_new,
            const std::vector<double>& y_old, const std::vector<double>& y_new,
            double t_interp);
    }

} // namespace ode45

#endif // ODE45_HPP