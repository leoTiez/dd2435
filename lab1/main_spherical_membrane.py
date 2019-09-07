from lab1 import *
import numpy as np
from matplotlib import pyplot as plt


def main():
    # Multiply by surface
    capacitance = 0.01
    diameter_soma = 100e-6
    diameter_dendrite = 1e-6
    soma_area = calc_surface(diameter_soma / 2.)
    dendrite_area = calc_surface(diameter_dendrite / 2.)
    capacitance_soma = capacitance * soma_area
    capacitance_dendrite = capacitance * dendrite_area
    voltage = -50e-3

    # Given parameters
    temperature = 293.
    k_conduct = 4.00e-9
    na_conduct = 0.12e-9
    cl_conduct = 0.40e-9

    k_concen_in = 400.
    na_concen_in = 50.
    cl_concen_in = 40.

    k_concen_out = 10.
    na_concen_out = 460.
    cl_concen_out = 5.

    k_current = ghk_current(
        voltage=voltage,
        permeability=k_conduct,
        concen_in=k_concen_in,
        concen_out=k_concen_out,
        valence=1.
    )

    na_current = ghk_current(
        voltage=voltage,
        permeability=na_conduct,
        concen_in=na_concen_in,
        concen_out=na_concen_out,
        valence=1.
    )

    cl_current = ghk_current(
        voltage=voltage,
        permeability=cl_conduct,
        concen_in=cl_concen_in,
        concen_out=cl_concen_out,
        valence=-1.
    )

    # Calculate current with constant voltage of -50mV
    total_current = calc_membrane_current(k_current, na_current, cl_current, 0.0, soma_area)

    print('Membrane current for a soma\n', total_current)

    rest_potential = ghk_voltage(
        k_conduct=k_conduct,
        na_conduct=na_conduct,
        cl_conduct=cl_conduct,
        k_concen_in=k_concen_in,
        na_concen_in=na_concen_in,
        cl_concen_in=cl_concen_in,
        k_concen_out=k_concen_out,
        na_concen_out=na_concen_out,
        cl_concen_out=cl_concen_out,
        temperature=temperature
    )

    print('Resting potential', rest_potential)

    # Calculate total chord conductance
    total_conductance = chord_conductance(
        current=total_current,
        voltage=voltage,
        equilibrium=rest_potential
    )
    print('Total membran conductance = I / (Vm - Er)\n', total_conductance)

    # No additional inject current or set it 0.015
    injected_current = 0.0
    membrane_current = calc_membrane_current(k_current, na_current, cl_current, injected_current, soma_area)
    time_frame = 0.050
    time_step = 0.1e-3
    time_values = np.arange(0.0, time_frame, time_step)
    voltage_values = []
    # Calculate evolution of voltage over time
    for _ in range(int(time_frame / time_step)):
        voltage = euler_integration_voltage(
            previous_voltage=voltage,
            capacitance=capacitance_soma,
            current=membrane_current,
            dt=time_step
        )
        voltage_values.append(voltage)

        k_current = ghk_current(
            voltage=voltage,
            permeability=k_conduct,
            concen_in=k_concen_in,
            concen_out=k_concen_out,
            valence=1.
        )

        na_current = ghk_current(
            voltage=voltage,
            permeability=na_conduct,
            concen_in=na_concen_in,
            concen_out=na_concen_out,
            valence=1.
        )

        cl_current = ghk_current(
            voltage=voltage,
            permeability=cl_conduct,
            concen_in=cl_concen_in,
            concen_out=cl_concen_out,
            valence=-1.
        )

        membrane_current = calc_membrane_current(k_current, na_current, cl_current, injected_current, soma_area)

    # Create plots
    plt.plot(time_values, voltage_values, 'b-', label='voltage')
    plt.plot(
        time_values,
        np.repeat(
            (voltage_values[0] - voltage_values[-1]) / np.e + voltage_values[-1],
            len(voltage_values)
        ),
        'k--', label='37% of the original value = tau')
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel('Time t')
    plt.ylabel('Voltage V')
    plt.show()

    # Changing permeability conditions within specific time frames
    # to model spike like behavior
    voltage_values_changed_cond = []
    voltage = -50e-3
    k_current = ghk_current(
        voltage=voltage,
        permeability=k_conduct,
        concen_in=k_concen_in,
        concen_out=k_concen_out,
        valence=1.
    )

    na_current = ghk_current(
        voltage=voltage,
        permeability=na_conduct,
        concen_in=na_concen_in,
        concen_out=na_concen_out,
        valence=1.
    )

    cl_current = ghk_current(
        voltage=voltage,
        permeability=cl_conduct,
        concen_in=cl_concen_in,
        concen_out=cl_concen_out,
        valence=-1.
    )

    membrane_current = calc_membrane_current(k_current, na_current, cl_current, injected_current, soma_area)

    # loop over time and set permeabilities depending ont the time step
    k_conduct = 4.e-9
    na_conduct = 0.12e-9
    for t in time_values:
        if 0.010 <= t < 0.015:
            k_conduct = 4.0e-9
            na_conduct = 6.0e-9
        elif 0.015 <= t < 0.025:
            k_conduct = 4.0e-9
            na_conduct = 0.12e-9
        elif 0.025 <= t < 0.030:
            k_conduct = 40.0e-9
            na_conduct = 0.12e-9
        elif 0.030 <= t:
            k_conduct = 4.0e-9
            na_conduct = 0.12e-9

        voltage = euler_integration_voltage(
            previous_voltage=voltage,
            capacitance=capacitance_soma,
            current=membrane_current,
            dt=time_step
        )
        voltage_values_changed_cond.append(voltage)

        k_current = ghk_current(
            voltage=voltage,
            permeability=k_conduct,
            concen_in=k_concen_in,
            concen_out=k_concen_out,
            valence=1.
        )

        na_current = ghk_current(
            voltage=voltage,
            permeability=na_conduct,
            concen_in=na_concen_in,
            concen_out=na_concen_out,
            valence=1.
        )

        cl_current = ghk_current(
            voltage=voltage,
            permeability=cl_conduct,
            concen_in=cl_concen_in,
            concen_out=cl_concen_out,
            valence=-1.
        )

        membrane_current = calc_membrane_current(k_current, na_current, cl_current, injected_current, soma_area)

    # Plot results
    plt.plot(time_values, voltage_values_changed_cond, 'b-', label='voltage')
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel('Time t')
    plt.ylabel('Voltage V')
    plt.show()

    injected_current = 0.0
    voltage = -50e-3
    k_conduct = 4.0e-9
    na_conduct = 0.12e-9
    k_current = ghk_current(
        voltage=voltage,
        permeability=k_conduct,
        concen_in=k_concen_in,
        concen_out=k_concen_out,
        valence=1.
    )

    na_current = ghk_current(
        voltage=voltage,
        permeability=na_conduct,
        concen_in=na_concen_in,
        concen_out=na_concen_out,
        valence=1.
    )

    cl_current = ghk_current(
        voltage=voltage,
        permeability=cl_conduct,
        concen_in=cl_concen_in,
        concen_out=cl_concen_out,
        valence=-1.
    )

    # Instead of assuming a soma with a diameter of 100microm
    # we assume a dendrite with a diameter of 1mircom
    membrane_current = calc_membrane_current(k_current, na_current, cl_current, injected_current, dendrite_area)
    print('Membrane current for a dendrite\n', membrane_current)
    time_frame = 0.050
    time_step = 0.1e-3
    time_values = np.arange(0.0, time_frame, time_step)
    voltage_values_dendrite = []
    for _ in range(int(time_frame / time_step)):
        voltage = euler_integration_voltage(
            previous_voltage=voltage,
            capacitance=capacitance_dendrite,
            current=membrane_current,
            dt=time_step
        )
        voltage_values_dendrite.append(voltage)

        k_current = ghk_current(
            voltage=voltage,
            permeability=k_conduct,
            concen_in=k_concen_in,
            concen_out=k_concen_out,
            valence=1.
        )

        na_current = ghk_current(
            voltage=voltage,
            permeability=na_conduct,
            concen_in=na_concen_in,
            concen_out=na_concen_out,
            valence=1.
        )

        cl_current = ghk_current(
            voltage=voltage,
            permeability=cl_conduct,
            concen_in=cl_concen_in,
            concen_out=cl_concen_out,
            valence=-1.
        )

        membrane_current = calc_membrane_current(k_current, na_current, cl_current, injected_current, dendrite_area)

    # Plot results
    plt.plot(time_values, voltage_values_dendrite, 'b-', label='voltage')
    plt.plot(
        time_values,
        np.repeat(
            (voltage_values_dendrite[0] - voltage_values_dendrite[-1]) / np.e + voltage_values_dendrite[-1],
            len(voltage_values_dendrite)
        ),
        'k--', label='37% of the original value = tau')
    plt.legend(loc='upper right')
    plt.grid()
    plt.xlabel('Time t')
    plt.ylabel('Voltage V')
    plt.show()


if __name__ == '__main__':
    main()

