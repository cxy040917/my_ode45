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

    // ����΢�ַ�������: dy/dt = f(t, y)
    using ODEFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

    // �����¼���������
    using EventFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

    // ���������������
    using OutputFcn = std::function<bool(double, const std::vector<double>&, const std::vector<double>&)>;

    // �¼���Ϣ�ṹ��
    struct EventInfo {
        size_t event_index;      // �¼�����
        double event_time;       // �¼�����ʱ��
        std::vector<double> event_y; // �¼�����ʱ��״̬
        bool is_terminal;        // �Ƿ�Ϊ��ֹ�¼�
        int direction;           // �¼����� (0=����, 1=��, -1=��)
    };

    // �������ṹ��
    struct ODEResult {
        std::vector<double> t;
        std::vector<std::vector<double>> y;
        std::vector<EventInfo> events;
    };

    // Ĭ���ݲ�
    constexpr double DEFAULT_RTOL = 1e-6;
    constexpr double DEFAULT_ATOL = 1e-8;

    // ���ѡ��ṹ��
    struct Options {
        double rtol = DEFAULT_RTOL;
        double atol = DEFAULT_ATOL;
        double max_step = std::numeric_limits<double>::infinity();
        double initial_step = 0.0;
        EventFunction event_fcn = nullptr;
        OutputFcn output_fcn = nullptr;
        std::vector<int> event_directions; // ÿ���¼��ķ��� (0=����, 1=��, -1=��)
        std::vector<bool> event_terminal; // ÿ���¼��Ƿ�Ϊ��ֹ�¼�
    };

    // ����⺯��
    ODEResult solve(
        const ODEFunction& func,
        const std::vector<double>& tspan,
        const std::vector<double>& y0,
        const Options& options = Options());

    // ��������
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

        // @brief ����¼���������
        // 1. ��ʼ״̬��������|val_old| > tolerance��
        // 2. ֵ�������ű仯��val_old * val_new <= 0��
        // 3. ����ƥ�䣨direction == 0 ����б��һ�£�
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