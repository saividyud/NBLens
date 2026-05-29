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
    # plt.rcParams['font.family'] = 'Times New Roman'
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

    #%% Importing simulation data
    multiplanet_directory = './Unity/Simulations/Part_6_Multiplanet_Collection'
    sun_jupiter_file = f'{multiplanet_directory}/sun_jupiter.pkl'
    sun_jupiter_sim = IRSC.caustic_reader(sun_jupiter_file)
    sun_jupiter_saturn_file = f'{multiplanet_directory}/sun_jupiter_saturn.pkl'
    sun_jupiter_saturn_sim = IRSC.caustic_reader(sun_jupiter_saturn_file)
    sun_jupiter_saturn_uranus_file = f'{multiplanet_directory}/sun_jupiter_saturn_uranus.pkl'
    sun_jupiter_saturn_uranus_sim = IRSC.caustic_reader(sun_jupiter_saturn_uranus_file)
    sun_jupiter_saturn_uranus_earth_file = f'{multiplanet_directory}/sun_jupiter_saturn_uranus_earth.pkl'
    sun_jupiter_saturn_uranus_earth_sim = IRSC.caustic_reader(sun_jupiter_saturn_uranus_earth_file)

    #%% Extracting important parameters from simulations
    pixels = sun_jupiter_sim.pixels
    ang_width = sun_jupiter_sim.ang_width
    ang_res = sun_jupiter_sim.ang_res

    sun_jupiter_magnifications = sun_jupiter_sim.magnifications
    sun_jupiter_saturn_magnifications = sun_jupiter_saturn_sim.magnifications
    sun_jupiter_saturn_uranus_magnifications = sun_jupiter_saturn_uranus_sim.magnifications
    sun_jupiter_saturn_uranus_earth_magnifications = sun_jupiter_saturn_uranus_earth_sim.magnifications

    #%% Plotting magnification map
    sun_jupiter_sim.plot()
    plt.savefig(f'./Unity Analysis/Part 6/Figures/sun_jupiter_magnification_map.png', dpi=100)

    sun_jupiter_saturn_sim.plot()
    plt.savefig(f'./Unity Analysis/Part 6/Figures/sun_jupiter_saturn_magnification_map.png', dpi=100)

    sun_jupiter_saturn_uranus_sim.plot()
    plt.savefig(f'./Unity Analysis/Part 6/Figures/sun_jupiter_saturn_uranus_magnification_map.png', dpi=100)

    sun_jupiter_saturn_uranus_earth_sim.plot()
    plt.savefig(f'./Unity Analysis/Part 6/Figures/sun_jupiter_saturn_uranus_earth_magnification_map.png', dpi=100)

    #%% Plotting fractional difference maps
    frac_diff_jupiter_saturn = (sun_jupiter_saturn_magnifications - sun_jupiter_magnifications) / sun_jupiter_magnifications
    frac_diff_jupiter_saturn_uranus = (sun_jupiter_saturn_uranus_magnifications - sun_jupiter_saturn_magnifications) / sun_jupiter_saturn_magnifications
    frac_diff_jupiter_saturn_uranus_earth = (sun_jupiter_saturn_uranus_earth_magnifications - sun_jupiter_saturn_uranus_magnifications) / sun_jupiter_saturn_uranus_magnifications

    fig = plt.figure()
    fig.set_constrained_layout(True)
    ax = fig.add_subplot()

    fig.suptitle('Fractional Difference Between Sun-Jupiter and Sun-Jupiter-Saturn Magnification Maps')

    plot = ax.imshow(frac_diff_jupiter_saturn, 
                    cmap='gray',
                    extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2],
    )
    plot = ax.contour(np.flip(frac_diff_jupiter_saturn, axis=0),
            levels=[-0.30, -0.10, -0.03, -0.01, -0.003, -0.001, 0, 0.001, 0.003, 0.01, 0.03, 0.10, 0.30],
            colors=['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'green', 'red', 'red', 'red', 'red', 'red', 'red'],
            linewidths=[1.3, 1.1, 0.9, 0.7, 0.5, 0.3, 0.5, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3],
            extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2]
    )

    bar = plt.colorbar(plot)
    bar.set_label('Fractional Difference')

    ax.set_xlabel(r'X [$\theta_E$]')
    ax.set_ylabel(r'Y [$\theta_E$]')

    plt.savefig(f'./Unity Analysis/Part 6/Figures/frac_diff_jupiter_saturn.png', dpi=100)

    fig = plt.figure()
    fig.set_constrained_layout(True)
    ax = fig.add_subplot()

    fig.suptitle('Fractional Difference Between Sun-Jupiter-Saturn and Sun-Jupiter-Saturn-Uranus Magnification Maps')

    plot = ax.imshow(frac_diff_jupiter_saturn_uranus, 
                    cmap='gray',
                    extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2],
    )
    plot = ax.contour(np.flip(frac_diff_jupiter_saturn_uranus, axis=0),
            levels=[-0.30, -0.10, -0.03, -0.01, -0.003, -0.001, 0, 0.001, 0.003, 0.01, 0.03, 0.10, 0.30],
            colors=['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'green', 'red', 'red', 'red', 'red', 'red', 'red'],
            linewidths=[1.3, 1.1, 0.9, 0.7, 0.5, 0.3, 0.5, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3],
            extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2]
    )

    bar = plt.colorbar(plot)
    bar.set_label('Fractional Difference')

    ax.set_xlabel(r'X [$\theta_E$]')
    ax.set_ylabel(r'Y [$\theta_E$]')

    plt.savefig(f'./Unity Analysis/Part 6/Figures/frac_diff_jupiter_saturn_uranus.png', dpi=100)

    fig = plt.figure()
    fig.set_constrained_layout(True)
    ax = fig.add_subplot()

    fig.suptitle('Fractional Difference Between Sun-Jupiter-Saturn-Uranus and Sun-Jupiter-Saturn-Uranus-Earth Magnification Maps')

    plot = ax.imshow(frac_diff_jupiter_saturn_uranus_earth, 
                    cmap='gray',
                    extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2],
    )
    plot = ax.contour(np.flip(frac_diff_jupiter_saturn_uranus_earth, axis=0),
            levels=[-0.30, -0.10, -0.03, -0.01, -0.003, -0.001, 0, 0.001, 0.003, 0.01, 0.03, 0.10, 0.30],
            colors=['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'green', 'red', 'red', 'red', 'red', 'red', 'red'],
            linewidths=[1.3, 1.1, 0.9, 0.7, 0.5, 0.3, 0.5, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3],
            extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2]
    )

    bar = plt.colorbar(plot)
    bar.set_label('Fractional Difference')

    ax.set_xlabel(r'X [$\theta_E$]')
    ax.set_ylabel(r'Y [$\theta_E$]')

    plt.savefig(f'./Unity Analysis/Part 6/Figures/frac_diff_jupiter_saturn_uranus_earth.png', dpi=100)