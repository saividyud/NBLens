if __name__ == '__main__':

    #%% Importing modules
    import sys
    sys.path.append('.')

    import numpy as np
    import pandas as pd
    import numexpr

    import IRSMicroLensing.IRSFunctions as IRSF
    import IRSMicroLensing.IRSCaustics as IRSC

    import csv
    import os

    import argparse

    import matplotlib.pyplot as plt

    #%% Setting plot parameters
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['figure.titlesize'] = 20
    plt.rcParams['figure.titleweight'] = 'bold'
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['figure.labelsize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['mathtext.fontset'] = 'cm'

    #%% Defining some useful functions
    def log_array(min_pow, max_pow):
        powers = np.arange(min_pow, max_pow + 1)
        base_values = np.array([1, 3])
        
        # Use broadcasting to create the sequence
        pos_arr = np.concatenate([base_values * (10.0**p) for p in powers])

        neg_arr = -pos_arr[::-1]

        return np.concatenate([neg_arr, np.array([0]), pos_arr])

    def colors_by_log(arr):
        colors = []
        for val in arr:
            if val < 0:
                colors.append('blue')
            elif val == 0:
                colors.append('green')
            elif val > 0:
                colors.append('red')
        return colors

    #%% Defining arguments
    parser = argparse.ArgumentParser(description='Plot microlensing simulations for triple lens systems.')
    parser.add_argument('--s1', type=str, default='1.0', help='Separation of big planet.')
    parser.add_argument('--s2', type=str, default='1.0', help='Separation of small planet.')
    parser.add_argument('--alpha', type=float, default='0.0', help='Angle of small planet.')
    args = parser.parse_args()
    s1_str = args.s1
    s2_str = args.s2
    alpha = args.alpha

    s1 = numexpr.evaluate(s1_str).item()
    s2 = numexpr.evaluate(s2_str).item()

    q1 = 1e-3
    q2 = 3e-4

    #%% Importing simulation data
    triple_directory = './Unity/Simulations/Part_5_Triple_Collection'
    binary_directory = './Unity/Simulations/Part_4_Binary_Collection'
    triple_file = f'{triple_directory}/triple_{s1:.2e}_{s2:.2e}_{alpha:.0f}.pkl'
    triple_sim = IRSC.caustic_reader(triple_file)
    binary_file = f'{binary_directory}/binary_{q1:.0e}_{s1:.2e}.pkl'
    binary_sim = IRSC.caustic_reader(binary_file)

    #%% Extracting important parameters from simulations
    pixels = triple_sim.pixels
    ang_width = triple_sim.ang_width
    ang_res = triple_sim.ang_res

    triple_magnifications = triple_sim.magnifications
    binary_magnifications = binary_sim.magnifications

    #%% Plotting magnification map
    triple_sim.plot()
    plt.savefig(f'./Unity Analysis/Part 5/Figures/triple_magnification_map_{s1:.2e}_{s2:.2e}_{alpha:.0f}.png', dpi=100)

    binary_sim.plot()
    plt.savefig(f'./Unity Analysis/Part 5/Figures/binary_magnification_map_{q1:.0e}_{s1:.2e}.png', dpi=100)

    #%% Plotting fractional difference map
    frac_diff_map = (triple_magnifications - binary_magnifications) / binary_magnifications

    fig = plt.figure()
    fig.set_constrained_layout(True)
    ax = fig.add_subplot()

    fig.suptitle('Fractional Difference Between Convolved Magnification Maps')

    plot = ax.imshow(frac_diff_map, 
                    cmap='gray',
                    extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2],
    )
    plot = ax.contour(np.flip(frac_diff_map, axis=0),
            levels=[-0.30, -0.10, -0.03, -0.01, -0.003, -0.001, 0, 0.001, 0.003, 0.01, 0.03, 0.10, 0.30],
            colors=['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'green', 'red', 'red', 'red', 'red', 'red', 'red'],
            linewidths=[1.3, 1.1, 0.9, 0.7, 0.5, 0.3, 0.5, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3],
            extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2]
    )

    bar = plt.colorbar(plot)
    bar.set_label('Fractional Difference')

    ax.set_xlabel(r'X [$\theta_E$]')
    ax.set_ylabel(r'Y [$\theta_E$]')

    plt.savefig(f'./Unity Analysis/Part 5/Figures/frac_diff_map_{s1:.2e}_{s2:.2e}_{alpha:.0f}.png', dpi=100)