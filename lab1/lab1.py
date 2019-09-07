import numpy as np
from scipy import constants as cs

FARADY_CONSTANT = cs.value(u'Faraday constant')
GAS_CONSTANT = cs.value(u'molar gas constant')


def calc_surface(radius=1):
    return 4 * np.pi * radius**2


def calc_membrane_current(i_k, i_na, i_cl, i_inject, surface):
    return (i_inject - i_k - i_na- i_cl) * surface


def alpha_m_test(voltage):
    return -0.1 * (voltage + 25.) / (np.exp(0.1 * voltage + 2.5) - 1.)


def beta_m_test(voltage):
    return 4.0 * np.exp(voltage / 18.0)


def alpha_m(voltage, error_margin=1e-8):
    if -0.035 - error_margin < voltage < -0.035 + error_margin:
        return 10e3
    else:
        return -10e5 * (voltage + 0.035) / (np.exp(-(voltage + 0.035) / 0.010) - 1.)


def beta_m(voltage):
    return 4000. * np.exp(-(voltage + 0.060) / 0.018)


def alpha_h(voltage):
    return 12. * np.exp(-voltage * .020)


def beta_h(voltage):
    return 180. / (np.exp(-(voltage + .030) / .010) + 1.)


def gating_steady_state(alpha_function, beta_function, voltage):
    alpha = alpha_function(voltage)
    beta = beta_function(voltage)
    return alpha / (alpha + beta)


def gating_dynamics(recent_state, alpha_function, beta_function, voltage, dt):
    alpha = alpha_function(voltage)
    beta = beta_function(voltage)
    return recent_state + (alpha * (1 - recent_state) - beta * recent_state) * dt


def gate_state_change(recent_state, alpha_function, beta_function, voltage, dt):
    alpha_dt = alpha_function(voltage) * dt
    beta_dt = beta_function(voltage) * dt
    probabilities = np.random.rand(recent_state.shape[0])
    next_state_1 = np.asarray(probabilities < alpha_dt, dtype=np.int) * np.asarray(recent_state == 0, dtype=np.int)
    next_state_0 = np.asarray(probabilities < beta_dt, dtype=np.int) * np.asarray(recent_state == 1, dtype=np.int)
    return recent_state + next_state_1 - next_state_0


def sodium_channel(recent_state, voltage, dt):
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
    return np.all(recent_state == 1, axis=0).astype(np.int)


def active_ratio(recent_state):
    return recent_state[recent_state == 1].shape[0] / float(recent_state.shape[0])


def xi(valence, voltage, temperature):
    return valence * voltage * FARADY_CONSTANT / (GAS_CONSTANT * temperature)


def equilibrium_potential(concen_in, concen_out, valence, temperature=293.):
    first_fact = GAS_CONSTANT * temperature / (valence * FARADY_CONSTANT)
    return first_fact * np.log(concen_out / concen_in)


def chord_conductance(current, voltage, equilibrium):
    return current / (voltage - equilibrium)


def tau(alpha_function, beta_function, voltage):
    return 1 / (alpha_function(voltage) + beta_function(voltage))


def euler_integration_voltage(previous_voltage, capacitance, current, dt):
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
        error_marging=1e-8
):
    if -1 * error_marging < voltage < error_marging:
        return permeability * valence * FARADY_CONSTANT * (concen_in - concen_out)
    else:
        xi_value = xi(valence, voltage, temperature)
        inner_parenth = (concen_in - concen_out * np.exp(-xi_value)) / (1 - np.exp(-xi_value))
        return permeability * valence * FARADY_CONSTANT * xi_value * inner_parenth

