from lab1 import *
import numpy as np
from matplotlib import pyplot as plt


def main():
    voltage_values = np.arange(-80.0e-3, 80.0e-3, 10.0e-3)
    m_values = [gating_steady_state(
        alpha_function=alpha_m,
        beta_function=beta_m,
        voltage=voltage
    ) for voltage in voltage_values]

    h_values = [gating_steady_state(
        alpha_function=alpha_h,
        beta_function=beta_h,
        voltage=voltage
    ) for voltage in voltage_values]

    tau_values = [tau(
        alpha_function=alpha_m,
        beta_function=beta_m,
        voltage=voltage
    ) for voltage in voltage_values]

    ax1 = plt.subplot(311)
    ax1.plot(voltage_values, m_values, 'b-', label='m')
    ax1.grid()
    ax1.legend(loc='upper right')

    ax2 = plt.subplot(312)
    ax2.plot(voltage_values, tau_values, 'g-', label='tau')
    ax2.grid()
    ax2.legend(loc='upper right')

    ax3 = plt.subplot(313)
    ax3.plot(voltage_values, h_values, 'r-', label='h')
    ax3.grid()
    ax3.legend(loc='upper right')

    plt.xlabel('Voltage V')
    plt.ylabel('m / h / tau')
    plt.show()

    dt = 0.1e-3
    time_frame = 0.05
    time_values = np.arange(0.0, time_frame, dt)
    num_of_states = 10000
    num_of_plots = 16
    fig, subplots = plt.subplots(num_of_plots, 1, sharex=True, sharey=True, figsize=(15, 10))

    for num, (voltage, ax) in enumerate(zip(voltage_values, subplots)):
        m_values = []
        recent_state = np.zeros(num_of_states)
        for _ in enumerate(time_values):
            recent_state = gate_state_change(
                recent_state=recent_state,
                alpha_function=alpha_m,
                beta_function=beta_m,
                voltage=voltage,
                dt=dt
            )

            m_values.append(active_ratio(recent_state))

        ax.plot(time_values, m_values, 'b-', label='Evolution of Sodium Channel at %f' % voltage)
        ax.plot(time_values, np.repeat(np.mean(np.asarray(m_values)), time_values.shape[0]),
                'r-',
                label='Mean of Sodium Channel at %f' % voltage)
        ax.legend(loc='upper right')

    plt.show()

    fig, subplots = plt.subplots(num_of_plots, 1, sharex=True, sharey=True, figsize=(15, 10))
    for num, (voltage, ax) in enumerate(zip(voltage_values, subplots)):
        recent_state = np.zeros((4, num_of_states))
        open_channels = []
        for _ in time_values:
            recent_state = sodium_channel(
                recent_state=recent_state,
                voltage=voltage,
                dt=dt
            )
            open_channels.append(active_ratio(evaluate_open_channel(recent_state)))

        print('Mean', np.mean(np.asarray(open_channels)), 'at', voltage)
        ax.plot(time_values, open_channels, 'b-', label='Evolution of Sodium Channel at %f' % voltage)
        ax.plot(time_values, np.repeat(np.mean(np.asarray(open_channels)), time_values.shape[0]),
                'r-',
                label='Mean of Sodium Channel at %f' % voltage)
        ax.legend(loc='upper right')
    plt.show()


if __name__ == '__main__':
    main()


