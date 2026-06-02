if __name__ == '__main__':

    #%% Importing modules
    import sys
    sys.path.append('.')

    import numpy as np
    import pandas as pd
    import numexpr

    import IRSMicroLensing.IRSFunctions as IRSF
    import IRSMicroLensing.IRSCaustics as IRSC

    import astropy.units as u
    import astropy.constants as c

    import csv
    import os

    import argparse

    import matplotlib.pyplot as plt

    #%% Defining arguments
    parser = argparse.ArgumentParser(description='Plot microlensing simulations for triple lens systems.')
    parser.add_argument('--radius', type=float, default=1.0, help='Radius of source profile (in Rsun).')
    args = parser.parse_args()
    radius = args.radius

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

    #%% Defining simulation parameters
    D_s = 8 * u.kpc
    D_l = 4 * u.kpc
    theta_ein, r_ein = IRSF.IRSFunctions.e_ring([D_s, D_l], 1)
    print(f'Einstein ring radius: {r_ein}')
    print(f'Einstein ring angle: {theta_ein}')

    #%% Importing simulation data
    print('Importing Sun-Jupiter simulation data...')
    multiplanet_directory = './Unity/Simulations/Part_6_Multiplanet_Collection'
    sun_jupiter_file = f'{multiplanet_directory}/sun_jupiter.pkl'
    sun_jupiter_sim = IRSC.caustic_reader(sun_jupiter_file)
    print(f'Pixels: {sun_jupiter_sim.pixels}')
    print(f'Angular width: {sun_jupiter_sim.ang_width}')
    print(f'Annulus lower bound: {sun_jupiter_sim.y_minus}')
    print(f'Annulus upper bound: {sun_jupiter_sim.y_plus}')
    print(f'Thickness: {sun_jupiter_sim.y_plus - sun_jupiter_sim.y_minus}')
    print(f'Number of rays in theta: {sun_jupiter_sim.num_theta}')
    print(f'Number of rays in r: {sun_jupiter_sim.num_r}')

    print('Importing Sun-Jupiter-Saturn simulation data...')
    sun_jupiter_saturn_file = f'{multiplanet_directory}/sun_jupiter_saturn.pkl'
    sun_jupiter_saturn_sim = IRSC.caustic_reader(sun_jupiter_saturn_file)
    print(f'Pixels: {sun_jupiter_saturn_sim.pixels}')
    print(f'Angular width: {sun_jupiter_saturn_sim.ang_width}')
    print(f'Annulus lower bound: {sun_jupiter_saturn_sim.y_minus}')
    print(f'Annulus upper bound: {sun_jupiter_saturn_sim.y_plus}')
    print(f'Thickness: {sun_jupiter_saturn_sim.y_plus - sun_jupiter_saturn_sim.y_minus}')
    print(f'Number of rays in theta: {sun_jupiter_saturn_sim.num_theta}')
    print(f'Number of rays in r: {sun_jupiter_saturn_sim.num_r}')

    print('Importing Sun-Jupiter-Saturn-Uranus simulation data...')
    sun_jupiter_saturn_uranus_file = f'{multiplanet_directory}/sun_jupiter_saturn_uranus.pkl'
    sun_jupiter_saturn_uranus_sim = IRSC.caustic_reader(sun_jupiter_saturn_uranus_file)
    print(f'Pixels: {sun_jupiter_saturn_uranus_sim.pixels}')
    print(f'Angular width: {sun_jupiter_saturn_uranus_sim.ang_width}')
    print(f'Annulus lower bound: {sun_jupiter_saturn_uranus_sim.y_minus}')
    print(f'Annulus upper bound: {sun_jupiter_saturn_uranus_sim.y_plus}')
    print(f'Thickness: {sun_jupiter_saturn_uranus_sim.y_plus - sun_jupiter_saturn_uranus_sim.y_minus}')
    print(f'Number of rays in theta: {sun_jupiter_saturn_uranus_sim.num_theta}')
    print(f'Number of rays in r: {sun_jupiter_saturn_uranus_sim.num_r}')

    print('Importing Sun-Jupiter-Saturn-Uranus-Earth simulation data...')
    sun_jupiter_saturn_uranus_earth_file = f'{multiplanet_directory}/sun_jupiter_saturn_uranus_earth.pkl'
    sun_jupiter_saturn_uranus_earth_sim = IRSC.caustic_reader(sun_jupiter_saturn_uranus_earth_file)
    print(f'Pixels: {sun_jupiter_saturn_uranus_earth_sim.pixels}')
    print(f'Angular width: {sun_jupiter_saturn_uranus_earth_sim.ang_width}')
    print(f'Annulus lower bound: {sun_jupiter_saturn_uranus_earth_sim.y_minus}')
    print(f'Annulus upper bound: {sun_jupiter_saturn_uranus_earth_sim.y_plus}')
    print(f'Thickness: {sun_jupiter_saturn_uranus_earth_sim.y_plus - sun_jupiter_saturn_uranus_earth_sim.y_minus}')
    print(f'Number of rays in theta: {sun_jupiter_saturn_uranus_earth_sim.num_theta}')
    print(f'Number of rays in r: {sun_jupiter_saturn_uranus_earth_sim.num_r}')

    #%% Extracting important parameters from simulations
    print('Extracting important parameters from simulations...')
    pixels = sun_jupiter_sim.pixels
    ang_width = sun_jupiter_sim.ang_width
    ang_res = sun_jupiter_sim.ang_res

    print('Extracting magnifications from simulations...')
    sun_jupiter_magnifications = sun_jupiter_sim.magnifications
    sun_jupiter_saturn_magnifications = sun_jupiter_saturn_sim.magnifications
    sun_jupiter_saturn_uranus_magnifications = sun_jupiter_saturn_uranus_sim.magnifications
    sun_jupiter_saturn_uranus_earth_magnifications = sun_jupiter_saturn_uranus_earth_sim.magnifications

    #%% Plotting magnification map
    # print('Plotting Sun-Jupiter Magnification Map...')
    # sun_jupiter_sim.plot()
    # plt.savefig(f'./Unity Analysis/Part 6/Figures/sun_jupiter_magnification_map.png', dpi=300)

    # print('Plotting Sun-Jupiter-Saturn Magnification Map...')
    # sun_jupiter_saturn_sim.plot()
    # plt.savefig(f'./Unity Analysis/Part 6/Figures/sun_jupiter_saturn_magnification_map.png', dpi=300)

    # print('Plotting Sun-Jupiter-Saturn-Uranus Magnification Map...')
    # sun_jupiter_saturn_uranus_sim.plot()
    # plt.savefig(f'./Unity Analysis/Part 6/Figures/sun_jupiter_saturn_uranus_magnification_map.png', dpi=300)

    # print('Plotting Sun-Jupiter-Saturn-Uranus-Earth Magnification Map...')
    # sun_jupiter_saturn_uranus_earth_sim.plot()
    # plt.savefig(f'./Unity Analysis/Part 6/Figures/sun_jupiter_saturn_uranus_earth_magnification_map.png', dpi=300)

    #%% Defining source profile
    print('Defining source profile...')
    source_radius_solar = radius * c.R_sun.to(u.m).value # meters
    source_radius_angular = source_radius_solar / D_l.to(u.m).value * u.rad # radians
    print(f'Source radius: {source_radius_angular.to(u.mas)}')
    source_radius = source_radius_angular.to(u.mas).value / theta_ein.to(u.mas).value # theta_E
    print(f'Source radius: {source_radius}')

    source_profile = IRSF.IRSFunctions.source_profile(ang_res=sun_jupiter_sim.ang_res, profile_type='LD', rad=source_radius, LD=0.5)

    #%% Convolving magnification maps
    print('Convolving Sun-Jupiter Magnification Map...')
    sun_jupiter_convolved_magnifications = sun_jupiter_sim.convolve(source_profile=source_profile)

    print('Convolving Sun-Jupiter-Saturn Magnification Map...')
    sun_jupiter_saturn_convolved_magnifications = sun_jupiter_saturn_sim.convolve(source_profile=source_profile)

    print('Convolving Sun-Jupiter-Saturn-Uranus Magnification Map...')
    sun_jupiter_saturn_uranus_convolved_magnifications = sun_jupiter_saturn_uranus_sim.convolve(source_profile=source_profile)

    print('Convolving Sun-Jupiter-Saturn-Uranus-Earth Magnification Map...')
    sun_jupiter_saturn_uranus_earth_convolved_magnifications = sun_jupiter_saturn_uranus_earth_sim.convolve(source_profile=source_profile)

    #%% Plotting fractional difference maps
    print('Plotting Fractional Difference Between Sun-Jupiter and Sun-Jupiter-Saturn Magnification Maps...')
    frac_diff_jupiter_saturn = (sun_jupiter_saturn_convolved_magnifications - sun_jupiter_convolved_magnifications) / sun_jupiter_convolved_magnifications

    fig = plt.figure()
    fig.set_constrained_layout(True)
    ax = fig.add_subplot()

    fig.suptitle('Fractional Difference Between Sun-Jupiter and Sun-Jupiter-Saturn Magnification Maps')
    ax.set_title(f'Source Radius: {radius:.2e} $R_\odot$')

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

    plt.savefig(f'./Unity Analysis/Part 6/Figures/frac_diff_jupiter_saturn_{radius:.2e}.png', dpi=300)

    print('Plotting Fractional Difference Between Sun-Jupiter-Saturn and Sun-Jupiter-Saturn-Uranus Magnification Maps...')
    frac_diff_jupiter_saturn_uranus = (sun_jupiter_saturn_uranus_convolved_magnifications - sun_jupiter_saturn_convolved_magnifications) / sun_jupiter_saturn_convolved_magnifications
    
    fig = plt.figure()
    fig.set_constrained_layout(True)
    ax = fig.add_subplot()

    fig.suptitle('Fractional Difference Between Sun-Jupiter-Saturn and Sun-Jupiter-Saturn-Uranus Magnification Maps')
    ax.set_title(f'Source Radius: {radius:.2e} $R_\odot$')

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

    plt.savefig(f'./Unity Analysis/Part 6/Figures/frac_diff_jupiter_saturn_{radius:.2e}.png', dpi=300)

    print('Plotting Fractional Difference Between Sun-Jupiter-Saturn-Uranus and Sun-Jupiter-Saturn-Uranus-Earth Magnification Maps...')
    frac_diff_jupiter_saturn_uranus_earth = (sun_jupiter_saturn_uranus_earth_convolved_magnifications - sun_jupiter_saturn_uranus_convolved_magnifications) / sun_jupiter_saturn_uranus_convolved_magnifications

    fig = plt.figure()
    fig.set_constrained_layout(True)
    ax = fig.add_subplot()

    fig.suptitle('Fractional Difference Between Sun-Jupiter-Saturn-Uranus and Sun-Jupiter-Saturn-Uranus-Earth Magnification Maps')
    ax.set_title(f'Source Radius: {radius:.2e} $R_\odot$')

    plot = ax.imshow(frac_diff_jupiter_saturn_uranus_earth, 
                    cmap='gray',
                    extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2],
    )
    plot = ax.contour(np.flip(frac_diff_jupiter_saturn_uranus_earth, axis=0),
            levels=[-0.30, -0.10, -0.03, -0.01, -0.003, -0.001, -0.0003, -0.0001, 0, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.10, 0.30],
            colors=['blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'green', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red'],
            linewidths=[1.3, 1.1, 0.9, 0.7, 0.6, 0.5, 0.4, 0.3, 0.5, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.1, 1.3],
            extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2]
    )

    bar = plt.colorbar(plot)
    bar.set_label('Fractional Difference')

    ax.set_xlabel(r'X [$\theta_E$]')
    ax.set_ylabel(r'Y [$\theta_E$]')

    plt.savefig(f'./Unity Analysis/Part 6/Figures/frac_diff_jupiter_saturn_{radius:.2e}.png', dpi=300)