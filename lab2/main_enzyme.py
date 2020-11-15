from lab2 import *
import numpy as np
import matplotlib.pyplot as plt


def main():
    eq_const = calc_equilibrium_constant(cal_to_joule(-3.5e3))
    print('Equilibrium const', eq_const)

    concentrations = np.arange(0, 20, 0.1)

    fluxes = []
    for conc in concentrations:
        fluxes.append(flux(conc, 0))

    plt.title('Concentration-activity profile')
    plt.xlabel('Concentration')
    plt.ylabel('Reaction rate j')
    plt.plot(concentrations, fluxes)
    plt.show()

    mod_concentrations = np.arange(0, 1e-3, 2e-4)
    for mod_conc in mod_concentrations:
        fluxes = []
        for conc in concentrations:
            fluxes.append(flux(conc, mod_conc))
        plt.plot(concentrations, fluxes, label='%f' % mod_conc)
    plt.legend(loc='upper right')
    plt.xlabel('Concentration')
    plt.ylabel('Reaction rate j')
    plt.title('Concentration-activity profile over different concentrations of FBP')
    plt.show()


if __name__ == '__main__':
    main()