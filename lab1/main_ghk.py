from lab1 import *
import numpy as np
import matplotlib.pyplot as plt


def main():
    print('\nSolutions to Goldman-Hodgkin-Katz equation tasks')
    temperature = 293.
    k_conduct = 4.00e-9
    na_conduct = 0.12e-9
    cl_conduct = 0.40e-9

    k_concen_in = 400
    na_concen_in = 50
    cl_concen_in = 40

    k_concen_out = 10
    na_concen_out = 460
    cl_concen_out = 5

    # Use given values
    print('Voltage for given Values')
    print(ghk_voltage(
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
    ))

    # Interchange inside and outside concentrations
    print('\nVoltage for inverted concentrations')
    print(ghk_voltage(
        k_conduct=k_conduct,
        na_conduct=na_conduct,
        cl_conduct=cl_conduct,
        k_concen_in=k_concen_out,
        na_concen_in=na_concen_out,
        cl_concen_in=cl_concen_out,
        k_concen_out=k_concen_in,
        na_concen_out=na_concen_in,
        cl_concen_out=cl_concen_in,
        temperature=temperature
    ))

    # only potassium permeability non-zero
    print('\nVoltage if only permeability of potassium non-zero')
    print(ghk_voltage(
        k_conduct=k_conduct,
        na_conduct=0.,
        cl_conduct=0.,
        k_concen_in=k_concen_in,
        na_concen_in=na_concen_in,
        cl_concen_in=cl_concen_in,
        k_concen_out=k_concen_out,
        na_concen_out=na_concen_out,
        cl_concen_out=cl_concen_out,
        temperature=temperature
    ))

    # only sodium permeability non-zero
    print('\nVoltage if only permeability of sodium non-zero')
    print(ghk_voltage(
        k_conduct=0.,
        na_conduct=na_conduct,
        cl_conduct=0.,
        k_concen_in=k_concen_in,
        na_concen_in=na_concen_in,
        cl_concen_in=cl_concen_in,
        k_concen_out=k_concen_out,
        na_concen_out=na_concen_out,
        cl_concen_out=cl_concen_out,
        temperature=temperature
    ))

    # Calculate the I-V dependence
    # dimension of I is A / m2 and is the  current density
    voltage_values = np.arange(-80.0e-3, 80.0e-3, 5.0e-3)
    k_current_values = []
    na_current_values = []
    cl_current_values = []
    current_values = []
    error_margin = 1e-8
    for voltage in voltage_values:
        print_flag = -70e-3 + error_margin > voltage > -70e-3 - error_margin
        print_flag = print_flag or 0. + error_margin > voltage > 0. - error_margin
        if print_flag:
            print('\nVoltage', voltage)
        # Calculate current of potassium
        k_current = ghk_current(
            voltage=voltage,
            permeability=k_conduct,
            concen_in=k_concen_in,
            concen_out=k_concen_out,
            valence=1.
        )
        k_current_values.append(k_current)
        if print_flag:
            print('\nPotassium current')
            print(k_current)

        # Calculate current of sodium
        na_current = ghk_current(
            voltage=voltage,
            permeability=na_conduct,
            concen_in=na_concen_in,
            concen_out=na_concen_out,
            valence=1.
        )
        na_current_values.append(na_current)

        if print_flag:
            print('\nSodium current')
            print(na_current)

        # Calculate current of chloride
        cl_current = ghk_current(
            voltage=voltage,
            permeability=cl_conduct,
            concen_in=cl_concen_in,
            concen_out=cl_concen_out,
            valence=-1.
        )
        cl_current_values.append(cl_current)

        if print_flag:
            print('\nChloride current')
            print(cl_current)

        # Total current
        current_values.append(k_current + na_current + cl_current)
        if print_flag:
            print('\nTotal current', current_values[-1])

    plt.plot(voltage_values, current_values)
    plt.ylabel('Current I')
    plt.xlabel('Voltage V')
    plt.title('I-V relationship of the membrane')
    plt.grid()
    plt.show()


if __name__ == '__main__':
    main()

