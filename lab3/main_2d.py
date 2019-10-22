#!/usr/bin/python3
from reaction_diffusion import *
import numpy as np
import matplotlib.pyplot as plt


def random_state(
        shape=(50, 50),
        time_step=0.01,
        number_timesteps=3000,
        rand_upper_bound=0.1,
        snap_shot_rate=100
):
    time_array = np.arange(0, number_timesteps * time_step, time_step)

    deviation_a = np.random.rand(shape[0], shape[1]) * rand_upper_bound
    deviation_b = np.random.rand(shape[0], shape[1]) * rand_upper_bound

    for num, _ in enumerate(time_array):
        deviation_a_update, _ = react_diff(
            deviation_a,
            deviation_b,
            diffusion_coef=1.,
            dt=time_step,
            is_1d=False
        )
        deviation_b_update, _ = react_diff(
            deviation_a,
            deviation_b,
            diffusion_coef=3.,
            is_a_substance=False,
            dt=time_step,
            is_1d=False
        )
        deviation_a = deviation_a_update
        deviation_b = deviation_b_update

        if num % snap_shot_rate == 0:
            plt.imshow(deviation_a)
            plt.pause(0.1)
            plt.draw()

    plt.show()


def single_high_a_state(
        shape=(50, 20),
        time_step=0.01,
        num_changed_states=1,
        number_timesteps=3000,
        snap_shot_rate=100,
        initial_value=0.14
):
    time_array = np.arange(0, number_timesteps * time_step, time_step)

    deviation_a = np.zeros(shape)
    deviation_b = np.zeros(shape)

    for _ in range(num_changed_states):
        index = np.random.randint(shape[0]), np.random.randint(shape[1])
        while deviation_a[index] > 0:
            index = np.random.randint(shape[0]), np.random.randint(1)
        deviation_a[index] = initial_value

    for num, _ in enumerate(time_array):
        deviation_a_update, _ = react_diff(
            deviation_a,
            deviation_b,
            diffusion_coef=1.,
            dt=time_step,
            is_1d=False
        )
        deviation_b_update, _ = react_diff(
            deviation_a,
            deviation_b,
            diffusion_coef=3.,
            is_a_substance=False,
            dt=time_step,
            is_1d=False
        )
        deviation_a = deviation_a_update
        deviation_b = deviation_b_update

        if num % snap_shot_rate == 0:
            plt.figure(1)
            plt.imshow(deviation_a)
            plt.pause(0.1)
            plt.draw()

            plt.figure(2)
            plt.imshow(deviation_b)
            plt.pause(0.1)
            plt.draw()

    plt.show()


if __name__ == '__main__':
    random_state(shape=(50, 50), number_timesteps=7000, snap_shot_rate=200, time_step=0.005)
    random_state(shape=(50, 10), number_timesteps=5000, snap_shot_rate=1000)
    single_high_a_state(number_timesteps=5000, snap_shot_rate=1000, num_changed_states=1)
    single_high_a_state(
        shape=(50, 50),
        number_timesteps=7000,
        snap_shot_rate=100,
        num_changed_states=1,
        time_step=0.005
    )

