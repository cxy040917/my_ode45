#include "ode45.hpp"
#include <iostream>
#include <numeric>
#include <cassert>

namespace ode45 {

    ODEResult solve(
        const ODEFunction& func,
        const std::vector<double>& tspan,
        const std::vector<double>& y0,
        const Options& options) {

        // 检查输入参数
        if (tspan.size() != 2) {
            throw std::invalid_argument("tspan must have exactly 2 elements");
        }
        if (y0.empty()) {
            throw std::invalid_argument("Initial state y0 cannot be empty");
        }
        if (tspan[0] >= tspan[1]) {
            throw std::invalid_argument("tspan[0] must be less than tspan[1]");
        }

        const double t0 = tspan[0];
        const double tf = tspan[1];
        const size_t n = y0.size();

        ODEResult result;
        result.t.push_back(t0);
        result.y.push_back(y0);

        // 如果有输出函数，调用它
        if (options.output_fcn) {
            std::vector<double> yp = func(t0, y0);
            if (!options.output_fcn(t0, y0, yp)) {
                return result; // 用户请求终止
            }
        }

        // 计算初始步长
        double h = options.initial_step > 0 ? options.initial_step :
            internal::compute_initial_step(func, t0, y0, options.rtol, options.atol);
        h = std::min(h, options.max_step);

        double t = t0;
        std::vector<double> y = y0;

        // 如果有事件函数，检查初始状态
        if (options.event_fcn) {
            auto initial_events = internal::check_events(
                options.event_fcn, options.event_directions, options.event_terminal,
                t, t, y, y);
            for (const auto& event : initial_events) {
                result.events.push_back(event);
                if (event.is_terminal) {
                    return result; // 初始状态就满足终止条件
                }
            }
        }

        // 主循环
        while (t < tf) {
            // 确保最后一步不会超过tf
            if (t + h > tf) {
                h = tf - t;
            }

            // 尝试RK45步
            auto y_new = internal::rk45_step(func, t, y, h);

            // 估计误差
            auto y_low = internal::rk45_step(func, t, y, h / 2);
            y_low = internal::rk45_step(func, t + h / 2, y_low, h / 2);

            double error = internal::estimate_error(y, y_new, y_low, options.rtol, options.atol);

            // 接受或拒绝步长
            if (error <= 1.0) { // 步长可接受
                double t_old = t;
                auto y_old = y;

                t += h;
                y = y_new;

                // 检查事件
// 在ode45.cpp的solve函数中
                if (options.event_fcn) {
                    auto events = internal::check_events(
                        options.event_fcn, options.event_directions, options.event_terminal,
                        t_old, t, y_old, y);

                    for (const auto& event : events) {
                        // 修改后的插值调用
                        auto interpolated = internal::interpolate(
                            t_old, t, y_old, y, event.event_time);
                        double t_event = interpolated.first;
                        auto y_event = interpolated.second;

                        EventInfo event_info = event;
                        event_info.event_time = t_event;
                        event_info.event_y = y_event;
                        result.events.push_back(event_info);

                        if (event.is_terminal) {
                            result.t.push_back(t_event);
                            result.y.push_back(y_event);

                            if (options.output_fcn) {
                                std::vector<double> yp = func(t_event, y_event);
                                options.output_fcn(t_event, y_event, yp);
                            }
                            return result;
                        }
                    }
                }

                // 存储结果
                result.t.push_back(t);
                result.y.push_back(y);

                // 调用输出函数
                if (options.output_fcn) {
                    std::vector<double> yp = func(t, y);
                    if (!options.output_fcn(t, y, yp)) {
                        return result; // 用户请求终止
                    }
                }

                // 调整下一步的步长
                if (error > 0) {
                    h *= 0.9 * std::pow(1.0 / error, 1.0 / 5.0);
                    h = std::min(h, options.max_step);
                }
            }
            else { // 步长不可接受，减小步长重试
                h *= 0.9 * std::pow(1.0 / error, 1.0 / 5.0);
            }

            // 防止步长过小
            if (h < 16 * std::numeric_limits<double>::epsilon() * std::abs(t)) {
                throw std::runtime_error("Step size became too small");
            }
        }

        return result;
    }

    namespace internal {

        std::vector<double> rk45_step(
            const ODEFunction& func,
            double t,
            const std::vector<double>& y,
            double h) {

            const size_t n = y.size();

            // Runge-Kutta 4(5) coefficients
            constexpr double a2 = 1.0 / 4.0;
            constexpr double a3 = 3.0 / 8.0, b31 = 3.0 / 32.0, b32 = 9.0 / 32.0;
            constexpr double a4 = 12.0 / 13.0, b41 = 1932.0 / 2197.0, b42 = -7200.0 / 2197.0, b43 = 7296.0 / 2197.0;
            constexpr double a5 = 1.0, b51 = 439.0 / 216.0, b52 = -8.0, b53 = 3680.0 / 513.0, b54 = -845.0 / 4104.0;
            constexpr double a6 = 1.0 / 2.0, b61 = -8.0 / 27.0, b62 = 2.0, b63 = -3544.0 / 2565.0, b64 = 1859.0 / 4104.0, b65 = -11.0 / 40.0;

            constexpr double c1 = 25.0 / 216.0, c3 = 1408.0 / 2565.0, c4 = 2197.0 / 4104.0, c5 = -1.0 / 5.0;
            constexpr double d1 = 16.0 / 135.0, d3 = 6656.0 / 12825.0, d4 = 28561.0 / 56430.0, d5 = -9.0 / 50.0, d6 = 2.0 / 55.0;

            // 计算各个阶段的k值
            auto k1 = func(t, y);

            std::vector<double> y_temp(n);
            for (size_t i = 0; i < n; ++i) {
                y_temp[i] = y[i] + h * a2 * k1[i];
            }
            auto k2 = func(t + h * a2, y_temp);

            for (size_t i = 0; i < n; ++i) {
                y_temp[i] = y[i] + h * (b31 * k1[i] + b32 * k2[i]);
            }
            auto k3 = func(t + h * a3, y_temp);

            for (size_t i = 0; i < n; ++i) {
                y_temp[i] = y[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
            }
            auto k4 = func(t + h * a4, y_temp);

            for (size_t i = 0; i < n; ++i) {
                y_temp[i] = y[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
            }
            auto k5 = func(t + h * a5, y_temp);

            for (size_t i = 0; i < n; ++i) {
                y_temp[i] = y[i] + h * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
            }
            auto k6 = func(t + h * a6, y_temp);

            // 计算高阶解(4阶)
            std::vector<double> y_next(n);
            for (size_t i = 0; i < n; ++i) {
                y_next[i] = y[i] + h * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i]);
            }

            return y_next;
        }

        double estimate_error(
            const std::vector<double>& y,
            const std::vector<double>& y_high,
            const std::vector<double>& y_low,
            double rtol,
            double atol) {

            double error = 0.0;
            const size_t n = y.size();

            for (size_t i = 0; i < n; ++i) {
                double scale = atol + rtol * std::max(std::abs(y[i]), std::abs(y_high[i]));
                double e = std::abs(y_high[i] - y_low[i]) / scale;
                error += e * e;
            }

            return std::sqrt(error / n);
        }

        double compute_initial_step(
            const ODEFunction& func,
            double t0,
            const std::vector<double>& y0,
            double rtol,
            double atol) {

            const size_t n = y0.size();

            // 计算初始导数
            auto f0 = func(t0, y0);

            // 计算y的尺度
            double scale = 0.0;
            for (size_t i = 0; i < n; ++i) {
                double s = atol + rtol * std::abs(y0[i]);
                scale += (y0[i] / s) * (y0[i] / s);
            }
            scale = std::sqrt(scale / n);

            // 计算导数的尺度
            double dscale = 0.0;
            for (size_t i = 0; i < n; ++i) {
                double s = atol + rtol * std::abs(y0[i]);
                dscale += (f0[i] / s) * (f0[i] / s);
            }
            dscale = std::sqrt(dscale / n);

            if (scale < 1e-10 || dscale < 1e-10) {
                return 1e-6;
            }

            double h0 = 0.01 * scale / dscale;

            // 使用h0计算y1和f1
            std::vector<double> y1(n);
            for (size_t i = 0; i < n; ++i) {
                y1[i] = y0[i] + h0 * f0[i];
            }
            auto f1 = func(t0 + h0, y1);

            // 计算二阶导数的尺度
            double ddscale = 0.0;
            for (size_t i = 0; i < n; ++i) {
                double s = atol + rtol * std::abs(y0[i]);
                ddscale += ((f1[i] - f0[i]) / s) * ((f1[i] - f0[i]) / s);
            }
            ddscale = std::sqrt(ddscale / n) / h0;

            // 计算初始步长
            double h1;
            if (std::max(dscale, ddscale) <= 1e-15) {
                h1 = std::max(1e-6, h0 * 1e-3);
            }
            else {
                h1 = std::pow(0.01 / std::max(dscale, ddscale), 1.0 / 3.0);
            }

            return std::min(100.0 * h0, h1);
        }

        std::vector<EventInfo> check_events(
            const EventFunction& event_fcn,
            const std::vector<int>& directions,
            const std::vector<bool>& terminal,
            double t_old, double t_new,
            const std::vector<double>& y_old, const std::vector<double>& y_new) {

            std::vector<EventInfo> events;

            // 计算旧状态和新状态的事件值
            auto val_old = event_fcn(t_old, y_old);
            auto val_new = event_fcn(t_new, y_new);

            const size_t num_events = val_old.size();

            // 检查每个事件
            for (size_t i = 0; i < num_events; ++i) {
                int direction = (i < directions.size()) ? directions[i] : 0;
                bool is_terminal = (i < terminal.size()) ? terminal[i] : false;

                // 检查是否发生事件
                if (std::abs(val_old[i]) > 1e-12 && val_old[i] * val_new[i] <= 0.0) { // 符号变化
                    // 检查方向是否符合要求
                    double slope = (val_new[i] - val_old[i]) / (t_new - t_old);
                    if ((direction == 0) ||
                        (direction == 1 && slope > 0) ||
                        (direction == -1 && slope < 0)) {

                        // 使用线性插值估计事件时间
                        double t_event = t_old - val_old[i] * (t_new - t_old) / (val_new[i] - val_old[i]);

                        EventInfo event;
                        event.event_index = i;
                        event.event_time = t_event;
                        event.is_terminal = is_terminal;
                        event.direction = direction;

                        events.push_back(event);
                    }
                }
            }

            return events;
        }

        std::pair<double, std::vector<double>> internal::interpolate(
            double t_old, double t_new,
            const std::vector<double>& y_old, const std::vector<double>& y_new,
            double t_interp) {

            if (t_interp <= t_old) return { t_old, y_old };
            if (t_interp >= t_new) return { t_new, y_new };

            // 线性插值
            double theta = (t_interp - t_old) / (t_new - t_old);
            std::vector<double> y_interp(y_old.size());
            for (size_t i = 0; i < y_old.size(); ++i) {
                y_interp[i] = y_old[i] + theta * (y_new[i] - y_old[i]);
            }

            return { t_interp, y_interp };
        }

    } // namespace internal
} // namespace ode45