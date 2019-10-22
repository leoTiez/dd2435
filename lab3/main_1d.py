#!/usr/bin/python3
from reaction_diffusion import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors


def single_state_diff(
        number_of_cells=50,
        time_step=0.1,
        number_timesteps=3000,
        diff_equi_a=0.14,
        save_plots=True
):
    """
    Implements the behaviro if a single cell is not in equilibrium
    :param number_of_cells: Number of cells
    :param time_step: Margin of the time step
    :param number_timesteps: Number of time steps
    :param diff_equi_a: The deviation from equilibrium for the single cell for the concentration of species a
    :param save_plots: Flag to determine whether to save or to plot the figures
    :return: None
    """
    time_array = np.arange(0, number_timesteps * time_step, time_step)
    deviation_a = np.zeros(number_of_cells)
    deviation_b = np.zeros(number_of_cells)
    deviation_a[np.random.randint(number_of_cells)] = diff_equi_a

    history_a = [deviation_a]
    history_b = [deviation_b]
    for _ in time_array:
        deviation_a_update, _ = react_diff(
            deviation_a,
            deviation_b,
            diffusion_coef=1.,
            dt=time_step
        )
        deviation_b_update, _ = react_diff(
            deviation_a,
            deviation_b,
            diffusion_coef=3.,
            is_a_substance=False,
            dt=time_step
        )
        deviation_a = deviation_a_update
        deviation_b = deviation_b_update

        history_a.append(deviation_a)
        history_b.append(deviation_b)

    a_min, a_max = np.asarray(history_a).min(), np.asarray(history_a).max()
    b_min, b_max = np.asarray(history_b).min(), np.asarray(history_b).max()
    fig, (ax1, ax2) = plt.subplots(2, 1)
    pc_a = ax1.pcolor(np.asarray(history_a), vmin=a_min, vmax=a_max)
    fig.colorbar(pc_a, ax=ax1)
    pc_b = ax2.pcolor(np.asarray(history_b), vmin=b_min, vmax=b_max)
    fig.colorbar(pc_b, ax=ax2)

    if save_plots:
        plt.savefig('img/1d-single_state_diff.png')
    else:
        plt.show()


def random_state(
        number_of_cells=50,
        time_step=0.1,
        number_timesteps=3000,
        rand_upper_bound=0.1,
        save_plots=True
):
    """
    Sets all cells to a random initial state
    :param number_of_cells: Number of cells
    :param time_step: Margin of the time step
    :param number_timesteps: Number of time steps
    :param rand_upper_bound: Upper bound of the random values for the species
            that are assigned to the cells
    :param save_plots: Flag to determine whether to save or to plot the figures
    :return: None
    """
    time_array = np.arange(0, number_timesteps * time_step, time_step)

    deviation_a = np.random.rand(number_of_cells) * rand_upper_bound
    deviation_b = np.random.rand(number_of_cells) * rand_upper_bound
    history_a = [deviation_a]
    history_b = [deviation_b]
    for _ in time_array:
        deviation_a_update, _ = react_diff(
            deviation_a,
            deviation_b,
            diffusion_coef=1.,
            dt=time_step
        )
        deviation_b_update, _ = react_diff(
            deviation_a,
            deviation_b,
            diffusion_coef=3.,
            is_a_substance=False,
            dt=time_step
        )
        deviation_a = deviation_a_update
        deviation_b = deviation_b_update

        history_a.append(deviation_a)
        history_b.append(deviation_b)

    a_min, a_max = np.asarray(history_a).min(), np.asarray(history_a).max()
    b_min, b_max = np.asarray(history_b).min(), np.asarray(history_b).max()
    fig, (ax1, ax2) = plt.subplots(2, 1)
    pc_a = ax1.pcolor(np.asarray(history_a), vmin=a_min, vmax=a_max)
    fig.colorbar(pc_a, ax=ax1)
    pc_b = ax2.pcolor(np.asarray(history_b), vmin=b_min, vmax=b_max)
    fig.colorbar(pc_b, ax=ax2)

    if save_plots:
        plt.savefig('img/1d-random_state.png')
    else:
        plt.show()


def diffusion_coeff_change(
        number_of_cells=50,
        time_step=0.1,
        number_timesteps=3000,
        rand_upper_bound=0.1,
        diff_a_range_upper_bound=5,
        diff_b_range_upper_bound=4,
        diff_range_step=0.5,
        save_plots=True
):
    """
    Diffusion coefficients are gradually changed to determine impact on the system
    :param number_of_cells: Number of cells
    :param time_step: Margin of the time step
    :param number_timesteps: Number of the time steps
    :param rand_upper_bound: Upper bound of the random initial values
    :param diff_a_range_upper_bound: Upper bound of the diffusion coefficients for species a
    :param diff_b_range_upper_bound: Upper bound of the diffusion coefficients for species b
    :param diff_range_step: Step size between the diffusion coefficients
    :param save_plots: Flag to determine whether to save or to plot figures
    :return: None
    """
    time_array = np.arange(0, number_timesteps * time_step, time_step)

    diff_a_range = np.arange(0, diff_a_range_upper_bound, diff_range_step)
    diff_b_range = np.arange(0, diff_b_range_upper_bound, diff_range_step)
    for diff_a in diff_a_range:
        for diff_b in diff_b_range:
            deviation_a = np.random.rand(number_of_cells) * rand_upper_bound
            deviation_b = np.random.rand(number_of_cells) * rand_upper_bound
            history_a = [deviation_a]
            history_b = [deviation_b]
            for _ in time_array:
                deviation_a_update, _ = react_diff(
                    deviation_a,
                    deviation_b,
                    diffusion_coef=diff_a,
                    dt=time_step
                )
                deviation_b_update, _ = react_diff(
                    deviation_a,
                    deviation_b,
                    diffusion_coef=diff_b,
                    is_a_substance=False,
                    dt=time_step
                )
                deviation_a = deviation_a_update
                deviation_b = deviation_b_update

                history_a.append(deviation_a)
                history_b.append(deviation_b)

            a_min, a_max = np.asarray(history_a).min(), np.asarray(history_a).max()
            b_min, b_max = np.asarray(history_b).min(), np.asarray(history_b).max()
            fig, (ax1, ax2) = plt.subplots(2, 1)
            pc_a = ax1.pcolor(np.asarray(history_a), vmin=a_min, vmax=a_max)
            fig.colorbar(pc_a, ax=ax1)
            pc_b = ax2.pcolor(np.asarray(history_b), vmin=b_min, vmax=b_max)
            fig.colorbar(pc_b, ax=ax2)

            if save_plots:
                plt.savefig('img/1d-diffusion_coeff_change-{0}-{1}.png'.format(diff_a, diff_b))
            else:
                plt.show()


def change_interact_b(
        number_of_cells=50,
        time_step=0.1,
        number_timesteps=3000,
        rand_upper_bound=0.1,
        interact_b_start=-1,
        interact_b_end=-2,
        interact_range_step=-0.1,
        save_plots=True
):
    """
    Change influence of the concentration b
    :param number_of_cells: Number of cells
    :param time_step: Marging of the time step
    :param number_timesteps: Number of time steps
    :param rand_upper_bound: Upper bound for random initial states
    :param interact_b_start: Start value of the interaction coefficients of b
    :param interact_b_end: End value of the interaction coefficients of b
    :param interact_range_step: Step size of the interaction coefficients of b
    :param save_plots: Flag to determine whether to save or to plot figures
    :return: None
    """
    time_array = np.arange(0, number_timesteps * time_step, time_step)
    interact_b_range = np.arange(interact_b_start, interact_b_end, interact_range_step)

    for interact_b in interact_b_range:
        deviation_a = np.random.rand(number_of_cells) * rand_upper_bound
        deviation_b = np.random.rand(number_of_cells) * rand_upper_bound
        history_a = [deviation_a]
        history_b = [deviation_b]

        for _ in time_array:
            deviation_a_update, _ = react_diff(
                deviation_a,
                deviation_b,
                diffusion_coef=1.,
                interact_b=interact_b,
                dt=time_step
            )
            deviation_b_update, _ = react_diff(
                deviation_a,
                deviation_b,
                diffusion_coef=3.,
                is_a_substance=False,
                dt=time_step
            )
            deviation_a = deviation_a_update
            deviation_b = deviation_b_update

            history_a.append(deviation_a)
            history_b.append(deviation_b)

        a_min, a_max = np.asarray(history_a).min(), np.asarray(history_a).max()
        b_min, b_max = np.asarray(history_b).min(), np.asarray(history_b).max()
        fig, (ax1, ax2) = plt.subplots(2, 1)
        pc_a = ax1.pcolor(np.asarray(history_a), vmin=a_min, vmax=a_max)
        fig.colorbar(pc_a, ax=ax1)
        pc_b = ax2.pcolor(np.asarray(history_b), vmin=b_min, vmax=b_max)
        fig.colorbar(pc_b, ax=ax2)

        if save_plots:
            plt.savefig('img/1d-change_interact_b-{0}.png'.format(interact_b))
        else:
            plt.show()


def change_interact_a(
        number_of_cells=50,
        time_step=0.1,
        number_timesteps=3000,
        rand_upper_bound=0.1,
        interact_a_start=1,
        interact_a_end=2,
        interact_range_step=0.1,
        save_plots=True
):
    """
    Change influence of the concentration b
    :param number_of_cells: Number of cells
    :param time_step: Margin of time step
    :param number_timesteps: Number of time steps
    :param rand_upper_bound: Upper bound of random initial values
    :param interact_a_start: Start value of the interaction coefficients of a
    :param interact_a_end: End value of the interaction coefficients of a
    :param interact_range_step: Step size of the interaction coefficients of a
    :param save_plots: Flag to determine whether to save or to plot figures
    :return: None
    """
    time_array = np.arange(0, number_timesteps * time_step, time_step)

    interact_a_range = np.arange(interact_a_start, interact_a_end, interact_range_step)

    for interact_a in interact_a_range:
        deviation_a = np.random.rand(number_of_cells) * rand_upper_bound
        deviation_b = np.random.rand(number_of_cells) * rand_upper_bound
        history_a = [deviation_a]
        history_b = [deviation_b]
        for _ in time_array:
            deviation_a_update, _ = react_diff(
                deviation_a,
                deviation_b,
                interact_a=interact_a,
                diffusion_coef=1.,
                dt=time_step
            )
            deviation_b_update, _ = react_diff(
                deviation_a,
                deviation_b,
                diffusion_coef=3.,
                is_a_substance=False,
                dt=time_step
            )
            deviation_a = deviation_a_update
            deviation_b = deviation_b_update

            history_a.append(deviation_a)
            history_b.append(deviation_b)

        a_min, a_max = np.asarray(history_a).min(), np.asarray(history_a).max()
        b_min, b_max = np.asarray(history_b).min(), np.asarray(history_b).max()
        fig, (ax1, ax2) = plt.subplots(2, 1)
        pc_a = ax1.pcolor(np.asarray(history_a), vmin=a_min, vmax=a_max)
        fig.colorbar(pc_a, ax=ax1)
        pc_b = ax2.pcolor(np.asarray(history_b), vmin=b_min, vmax=b_max)
        fig.colorbar(pc_b, ax=ax2)

        if save_plots:
            plt.savefig('img/1d-change_interact_a-{0}.png'.format(interact_a))
        else:
            plt.show()


def phase_portrait(
        rand_upper_bound=0.1,
        interact_a_start=0,
        interact_a_end=2,
        interact_b_start=0,
        interact_b_end=2,
        interact_range_step=0.1,
        trials=3,
        save_plots=True
):
    """
    Plots the phase portrait of the diffusion reaction system
    :param rand_upper_bound: Upper bound of random initial values
    :param interact_a_start: Start value of the interaction coefficients of a
    :param interact_a_end: End value of the interaction coefficients of a
    :param interact_b_start: Start value of the interaction coefficients of b
    :param interact_b_end: End value of the interaction coefficients of b
    :param interact_range_step: Step size of the interaction coefficients
    :param trials: Number of trials since initial state is random
    :param save_plots: Flag to determine whether to save or to plot figures
    :return: None
    """

    interact_a_range = np.arange(interact_a_start, interact_a_end, interact_range_step)
    interact_b_range = np.arange(interact_b_start, interact_b_end, interact_range_step)

    A, B = np.meshgrid(interact_a_range, interact_b_range)

    phase_plane_vecs_a = []
    phase_plane_vecs_b = []
    plt.clf()
    for color in list(mcolors.BASE_COLORS)[:trials]:
        deviation_a = np.random.rand(1)*rand_upper_bound
        deviation_b = np.random.rand(1)*rand_upper_bound
        for interact_a in interact_a_range:
            for interact_b in interact_b_range:
                _, differential_a = react_diff(
                    deviation_a,
                    deviation_b,
                    interact_a=interact_a,
                    interact_b=interact_b
                )
                _, differential_b = react_diff(
                    deviation_a,
                    deviation_b,
                    diffusion_coef=3.,
                    is_a_substance=False
                )

                phase_plane_vecs_a.append(differential_a)
                phase_plane_vecs_b.append(differential_b)

        plt.quiver(A, B, phase_plane_vecs_a, phase_plane_vecs_b, color=color,
                   label='a: {0}, b: {1}'.format(deviation_a, deviation_b))

    plt.legend(loc='upper right')
    if save_plots:
        plt.savefig('img/1d-phase_portrait.png')
    else:
        plt.show()


if __name__ == '__main__':
    single_state_diff()
    random_state()
    diffusion_coeff_change()
    change_interact_b()
    change_interact_a()
    phase_portrait()

