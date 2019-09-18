import numpy as np
from scipy import constants as cs

FARADY_CONSTANT = cs.value(u'Faraday constant')
GAS_CONSTANT = cs.value(u'molar gas constant')
ERROR_MARGIN = 1e-8


def calc_surface(radius=1):
    """
    Calculate the sperical surface area
    :param radius: radius of the spere
    :return: surface
    """
    return 4 * np.pi * radius**2


def calc_membrane_current(i_k, i_na, i_cl, i_inject, surface):
    """
    Calculate the membrane current based on the ionic currents and the surface area
    :param i_k: Potassium current
    :param i_na: Sodium current
    :param i_cl: Chloride current
    :param i_inject: Injected current
    :param surface: Surface area
    :return: The denormalised membrane current in Ampere
    """
    return (i_inject - i_k - i_na- i_cl) * surface


def alpha_m_test(voltage):
    """
    Test function for the opening variable alpha of the m particle
    :param voltage: Potential
    :return: opening probability
    """
    return -0.1 * (voltage + 25.) / (np.exp(0.1 * voltage + 2.5) - 1.)


def beta_m_test(voltage):
    """
    Test function for the closing variable beta of the m particle
    :param voltage: Potential
    :return: closing probability
    """
    return 4.0 * np.exp(voltage / 18.0)


def alpha_m(voltage):
    """
    Opening variable alpha of the m particle as given in the assignment
    :param voltage: Potential
    :return: Opening probability
    """
    if -0.035 - ERROR_MARGIN < voltage < -0.035 + ERROR_MARGIN:
        return 10e3
    else:
        return -10e5 * (voltage + 0.035) / (np.exp(-(voltage + 0.035) / 0.010) - 1.)


def beta_m(voltage):
    """
    Closing behavior of the m particle as given in the assignment
    :param voltage: Potential
    :return: Closing probability
    """
    return 4000. * np.exp(-(voltage + 0.060) / 0.018)


def alpha_h(voltage):
    """
    Opening behavior of the h particle as given in the assignment
    :param voltage: Potential
    :return: Opening probability
    """
    return 12. * np.exp(-voltage * .020)


def beta_h(voltage):
    """
    Closing behavior of the h particle as given in the assingment
    :param voltage: Potential
    :return: Closing probability
    """
    return 180. / (np.exp(-(voltage + .030) / .010) + 1.)


def gating_steady_state(alpha_function, beta_function, voltage):
    """
    Steady state value of the gating variable
    :param alpha_function: Function modelling the opening variable alpha for the particle
    :param beta_function: Function modelling the closing variable alpha for the particle
    :param voltage: Potential
    :return: Steady state value of the particle
    """
    alpha = alpha_function(voltage)
    beta = beta_function(voltage)
    return alpha / (alpha + beta)


def gating_dynamics(recent_state, alpha_function, beta_function, voltage, dt):
    """
    Function desribing the dynamics using the Euler approximation
    :param recent_state: Recent state of the particle
    :param alpha_function: Function modelling the opening variable alpha for the particle
    :param beta_function:  Function modelling the closing variable alpha for the particle
    :param voltage: Potential
    :param dt: Time differential
    :return: New gating state
    """
    alpha = alpha_function(voltage)
    beta = beta_function(voltage)
    return recent_state + (alpha * (1 - recent_state) - beta * recent_state) * dt


def gate_state_change(recent_state, alpha_function, beta_function, voltage, dt):
    """
    Describing the state change of the partile (being either closed or open
    :param recent_state: State of the particle (binary, either closed or open
    :param alpha_function: Function modelling the opening variable alpha for the particle
    :param beta_function: Function modelling the closing variable alpha for the particle
    :param voltage: Potential
    :param dt: Time differential
    :return: New state of the gate (either closed or open)
    """
    alpha_dt = alpha_function(voltage) * dt
    beta_dt = beta_function(voltage) * dt
    probabilities = np.random.rand(recent_state.shape[0])
    next_state_1 = np.asarray(probabilities < alpha_dt, dtype=np.int) * np.asarray(recent_state == 0, dtype=np.int)
    next_state_0 = np.asarray(probabilities < beta_dt, dtype=np.int) * np.asarray(recent_state == 1, dtype=np.int)
    return recent_state + next_state_1 - next_state_0


def sodium_channel(recent_state, voltage, dt):
    """
    Modelling of the dynamics of a whole sodium channel
    :param recent_state: Recent state of the sodium channel (described by for binary values
                    3 m particles and 1 h particle)
    :param voltage: Potential
    :param dt: Time differential
    :return: New state of the sodium channel
    """
    new_state = np.zeros(recent_state.shape)
    for num, m_particle in enumerate(recent_state[:3, :]):
        new_state[num, :] = gate_state_change(
            recent_state=m_particle,
            alpha_function=alpha_m,
            beta_function=beta_m,
            voltage=voltage,
            dt=dt
        )

    new_state[3, :] = gate_state_change(
        recent_state=recent_state[3, :],
        alpha_function=alpha_h,
        beta_function=beta_h,
        voltage=voltage,
        dt=dt
    )

    return new_state


def evaluate_open_channel(recent_state):
    """
    Check whether channel is open or not
    :param recent_state: State of the channel
    :return: Open or closed channel
    """
    return np.all(recent_state == 1, axis=0).astype(np.int)


def active_ratio(recent_state):
    """
    Calculates ratio of active channels in comparison to all channels
    :param recent_state: channel states
    :return: ratio of active channels
    """
    return recent_state[recent_state == 1].shape[0] / float(recent_state.shape[0])


def xi(valence, voltage, temperature):
    """
    Calculates Xi as defined in GHK-I
    :param valence: Valence of the ion
    :param voltage: Potential
    :param temperature: Temperature
    :return: Xi value
    """
    return valence * voltage * FARADY_CONSTANT / (GAS_CONSTANT * temperature)


def equilibrium_potential(concen_in, concen_out, valence, temperature=293.):
    """
    Calculates equilibrium constant
    :param concen_in: Concentration inside the membrane
    :param concen_out: Concentration outside the membrane
    :param valence: Valence of the ion
    :param temperature: temperature
    :return: Equilibrium potential
    """
    first_fact = GAS_CONSTANT * temperature / (valence * FARADY_CONSTANT)
    return first_fact * np.log(concen_out / concen_in)


def chord_conductance(current, voltage, equilibrium):
    """
    Calculates chord conductance
    :param current: Recent current
    :param voltage: Recent potential
    :param equilibrium: Equilibrium potential
    :return: Chord conductance
    """
    return current / (voltage - equilibrium)


def tau(alpha_function, beta_function, voltage):
    """
    Calculates time constant tau of an specific particle (m, n, h)
    :param alpha_function: Function modelling the opening variable alpha for the particle
    :param beta_function: Function modelling the closing variable alpha for the particle
    :param voltage: Potential
    :return: Time constant tau
    """
    return 1 / (alpha_function(voltage) + beta_function(voltage))


def euler_integration_voltage(previous_voltage, capacitance, current, dt):
    """
    Euler integration of the voltage
    :param previous_voltage: Potential at the previous time step
    :param capacitance: Capacitance
    :param current: Current
    :param dt: Time differential
    :return: New voltage
    """
    return previous_voltage + (current / capacitance) * dt


def ghk_voltage(
        k_conduct,
        k_concen_in,
        k_concen_out,
        na_conduct,
        na_concen_in,
        na_concen_out,
        cl_conduct,
        cl_concen_in,
        cl_concen_out,
        temperature=293.
):
    """
    Voltage after the GHK model
    :param k_conduct: Potassium permeability
    :param k_concen_in: Potassium concentration inside the membrane
    :param k_concen_out: Potassium concentration outside the membrane
    :param na_conduct: Sodium permeability
    :param na_concen_in: Sodium concentration inside the membrane
    :param na_concen_out: Sodium concentration outside the membrane
    :param cl_conduct: Chloride permeability
    :param cl_concen_in: Chloride concentration inside the membrane
    :param cl_concen_out: Chloride concentration outside the membrane
    :param temperature: Temperature
    :return: Potential
    """
    pre_fact = (GAS_CONSTANT * temperature) / FARADY_CONSTANT
    inner_log = (k_conduct * k_concen_out + na_conduct * na_concen_out + cl_conduct * cl_concen_in) \
                / (k_conduct * k_concen_in + na_conduct * na_concen_in + cl_conduct * cl_concen_out)
    return pre_fact * np.log(inner_log)


def ghk_current(
        voltage,
        permeability,
        concen_in,
        concen_out,
        valence,
        temperature=293.,
):
    """
    Ion current after the GHK model
    :param voltage: Potential
    :param permeability: Permeability for the ion
    :param concen_in: Concentration of the ion inside the membrane
    :param concen_out: Concentration of the ion outside the membrane
    :param valence: Valence of the ion
    :param temperature: Temperature
    :return: Ion current
    """
    if -1 * ERROR_MARGIN < voltage < ERROR_MARGIN:
        return permeability * valence * FARADY_CONSTANT * (concen_in - concen_out)
    else:
        xi_value = xi(valence, voltage, temperature)
        inner_parenth = (concen_in - concen_out * np.exp(-xi_value)) / (1 - np.exp(-xi_value))
        return permeability * valence * FARADY_CONSTANT * xi_value * inner_parenth

