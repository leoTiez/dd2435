#!/usr/bin/python3
import numpy as np


def nabla_sq_1d(substance):
    """
    Second derivative of the vector function of the species for one dimension
    :param substance: species
    :return: second derivative of the species vector function for one dimension
    """
    substance_xn1 = np.roll(substance, -1)
    substance_xp1 = np.roll(substance, 1)
    return substance_xp1 + substance_xn1 - 2 * substance


def nabla_sq_2d(substance):
    """
    Second derivative of the vector function of the species for two dimensions
    :param substance: species
    :return: second derivative of the species vector function for two dimensions
    """
    substance_xn1 = np.roll(substance, -1, axis=0)
    substance_xp1 = np.roll(substance, 1, axis=0)
    substance_yn1 = np.roll(substance, -1, axis=1)
    substance_yp1 = np.roll(substance, 1, axis=1)
    return substance_xp1 + substance_xn1 + substance_yp1 + substance_yn1 - 4 * substance


def react_diff(
        deviation_a,
        deviation_b,
        interact_a=1.,
        interact_b=-1.,
        nonlin_break=0.1,
        diffusion_coef=1.,
        dt=0.1,
        is_a_substance=True,
        is_1d=True
):
    """
    Equation describing the dynamics of a reaction diffusion system with two
    species given by the assignment
    :param deviation_a: Deviation of the equilibrium of species a
    :param deviation_b: Deviation of the equilibrium of species b
    :param interact_a: Interaction coefficient of species a
    :param interact_b: Interaction coefficient of species b
    :param nonlin_break: Non-linear breakdown of the species
    :param diffusion_coef: Diffusion coefficient
    :param dt: Change in time
    :param is_a_substance: Flag to determine whether it is species a or b
    :param is_1d: Flag to determine whether reaction-diffusion takes place in one dimension
    :return: Updated substance value, expressed as a deviation from the equilibrium
    """
    if is_a_substance:
        substance = deviation_a
    else:
        substance = deviation_b
    if is_1d:
        nabla_sq = nabla_sq_1d(substance)
    else:
        nabla_sq = nabla_sq_2d(substance)
    delta = interact_a * deviation_a + interact_b * deviation_b - nonlin_break * substance**3 \
            + diffusion_coef * nabla_sq
    return substance + delta * dt, delta * dt

