import numpy as np
from scipy import constants as cs


GAS_CONSTANT = cs.value(u'molar gas constant')
STOICHIOMETRIC_MATRIX = np.asarray([[1, -1, 0], [0, 1, -1]])


def cal_to_joule(cal):
    """
    Conversion function from Calories to Joule (or kCal to kJ)
    :param cal: Energy in Calories
    :return: Energy in Joules
    """
    return cal * 4.1868


def calc_equilibrium_constant(thermo_pot, temperature=293.):
    """
    Caclulate the equilibrium constant
    :param thermo_pot: Thermodynamical potential or the difference of the Gibbs free energy
    :param temperature: Temperature
    :return: Equilibrium constant
    """
    return np.exp(-thermo_pot/(GAS_CONSTANT * temperature))


def sigma(subs_conc, michaelis_const=8):
    """
    Normalized substrate concentration
    :param subs_conc: Substrate concentration
    :param michaelis_const: Michaelis constant
    :return: Normalized substrate concentration
    """
    return subs_conc / float(michaelis_const)


def xi(modul_conc, dissociation_const=3e-3):
    """
    Normalized modulator concentration
    :param modul_conc: Modulator concentration
    :param dissociation_const: Dissociation constant
    :return: Normalized modulator concentration
    """

    return modul_conc / float(dissociation_const)


def modulator(modul_conc, cooperativity=2.5, effect=1/(0.1**2.5), dissociation_const=3e-3):
    """
    Modulator function
    :param modul_conc: Modulator concentration
    :param cooperativity: Cooperativity which is described by the Hill coefficient
    :param effect: Effect the modulator has on the enzyme
    :param dissociation_const: Dissociation constant
    :return: Modulator term
    """
    xi_value = xi(modul_conc, dissociation_const)
    return (1 + xi_value**cooperativity)/(1 + effect * xi_value**cooperativity)


def flux(
        subs_conc,
        modul_conc,
        limiting_rate=100/float(180),
        michaelis_const=8.,
        cooperativity=2.5,
        effect=1/(0.1**2.5),
        dissociation_const=3e-3,
        use_modulator=True
):
    """
    Calculate the flux j
    :param subs_conc: Substrate concentration
    :param modul_conc: Modulator concentration
    :param limiting_rate: Limiting rate
    :param michaelis_const: Michaelis constant
    :param cooperativity: Cooperativity which is described by the Hill coefficient
    :param effect: Effect the modulator has on the enzyme
    :param dissociation_const: Dissociation constant
    :param use_modulator Flag to determine whether or not to use the modulator
    :return: The flux j
    """
    sigma_value = sigma(subs_conc, michaelis_const=michaelis_const)
    if use_modulator:
        normalized_flux = (sigma_value**cooperativity) / (
                sigma_value**cooperativity
                + modulator(modul_conc,
                            cooperativity=cooperativity,
                            effect=effect,
                            dissociation_const=dissociation_const)
        )
    else:
        normalized_flux = sigma_value / (1 + sigma_value)

    return normalized_flux * limiting_rate


def flux_production(f6p_conc, fbp_conc, flux_param_dict, f6b_influx=0.6e-3, stoichiometric_matrix=None):
    """
    Computes the change in concentration over time
    :param fluxes: Fluxes for F6P and FBP
    :return:
    """
    pfk_param = flux_param_dict['pfk']
    pfk_flux = flux(
        f6p_conc,
        fbp_conc,
        limiting_rate=pfk_param['limiting_rate'],
        michaelis_const=pfk_param['michaelis_const'],
        cooperativity=pfk_param['cooperativity'],
        effect=pfk_param['effect'],
        dissociation_const=pfk_param['dissociation_const'],
        use_modulator=True
    )

    aldolase_param = flux_param_dict['aldolase']
    aldolase_flux = flux(
        fbp_conc,
        None,
        limiting_rate=aldolase_param['limiting_rate'],
        michaelis_const=aldolase_param['michaelis_const'],
        cooperativity=aldolase_param['cooperativity'],
        effect=aldolase_param['effect'],
        dissociation_const=aldolase_param['dissociation_const'],
        use_modulator=False
    )

    fluxes = np.asarray([f6b_influx, pfk_flux, aldolase_flux])
    if stoichiometric_matrix is None:
        stoichiometric_matrix = STOICHIOMETRIC_MATRIX

    return stoichiometric_matrix.dot(fluxes)


