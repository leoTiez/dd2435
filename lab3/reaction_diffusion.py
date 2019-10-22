#!/usr/bin/python3
import numpy as np


def nabla_sq_1d(substance):
    substance_xn1 = np.roll(substance, -1)
    substance_xp1 = np.roll(substance, 1)
    return substance_xp1 + substance_xn1 - 2 * substance


def nabla_sq_2d(substance):
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

