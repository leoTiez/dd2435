from scipy.integrate import ode
from lab2 import *
import matplotlib.pyplot as plt
import numpy as np


VERBOSITY = 0


def main():
    flux_param_dict = {
        'pfk': {
            'limiting_rate': 100 / float(180),
            'michaelis_const': 8.,
            'cooperativity': 2.5,
            'effect': 1 / (0.1 ** 2.5),
            'dissociation_const': 3e-3
        },
        'aldolase': {
            'limiting_rate': 60e-3,
            'michaelis_const': 10e-3,
            'cooperativity': 1,
            'effect': None,
            'dissociation_const': None
        }
    }

    def concentrations(time, conc):
        return flux_production(conc[0], conc[1], flux_param_dict=flux_param_dict)

    def concentrations_f6p_influx(time, conc):
        return flux_production(conc[0], conc[1], flux_param_dict=flux_param_dict, f6b_influx=6e-3)

    initial_values = [1., 1e-3]
    ode_simulation = ode(concentrations).set_integrator('vode', method='bdf', order=15)\
        .set_initial_value(initial_values)

    ode_simulation_influx = ode(concentrations_f6p_influx).set_integrator('vode', method='bdf', order=15) \
        .set_initial_value(initial_values)

    time = 2000
    f6p_conc = []
    fbp_conc = []

    f6p_conc_influx = []
    fbp_conc_influx = []
    for _ in range(time):
        ode_simulation.integrate(ode_simulation.t + 1)
        f6p_conc.append(ode_simulation.y[0])
        fbp_conc.append(ode_simulation.y[1])

        if VERBOSITY > 0:
            print('Result for time', ode_simulation.t, ':', ode_simulation.y)

    for _ in range(time):
        ode_simulation_influx.integrate(ode_simulation_influx.t + 1)
        f6p_conc_influx.append(ode_simulation_influx.y[0])
        fbp_conc_influx.append(ode_simulation_influx.y[1])

        if VERBOSITY > 0:
            print('Result for time', ode_simulation_influx.t, ':', ode_simulation_influx.y)

    time_range = np.arange(0, time, 1)
    plt.plot(time_range, f6p_conc, 'b-', label='F6P concentration low influx')
    plt.plot(time_range, fbp_conc, 'g-', label='FBP concentration low influx')
    plt.plot(time_range, f6p_conc_influx, 'm--', label='F6P concentration high influx')
    plt.plot(time_range, fbp_conc_influx, 'c--', label='FBP concentration high influx')
    plt.legend(loc='upper right')
    plt.title('Concentration of F6P and FBP over time')
    plt.xlabel('Time')
    plt.ylabel('Concentration')
    plt.show()


if __name__ == '__main__':
    main()

