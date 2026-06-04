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
    # Distance to source plane
    D_s = 8 * u.kpc

    # Distance to lens plane
    D_l = 4 * u.kpc

    theta_ein, r_ein = IRSF.IRSFunctions.e_ring([D_s, D_l], 1)
    print(f'Einstein ring radius: {r_ein}')
    print(f'Einstein ring angle: {theta_ein}')

    planets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']

    # Defining semimajor axes
    r_mercury = 0.387
    r_venus = 0.723
    r_earth = 1
    r_mars = 1.524
    r_jupiter = 5.2
    r_saturn = 9.5
    r_uranus = 19.2
    r_neptune = 30.07

    planet_semimajor_axes = np.array([r_mercury, r_venus, r_earth, r_mars, r_jupiter, r_saturn, r_uranus, r_neptune]) * u.au # au

    # Calculating in terms of theta_E
    planet_ss = (planet_semimajor_axes / r_ein).to(u.dimensionless_unscaled).value # theta_E

    # Defining masses
    m_mercury = 0.330103e24
    m_venus = 4.86731e24
    m_earth = 5.97217e24
    m_mars = 0.641691e24
    m_jupiter = 1898.125e24
    m_saturn = 568.317e24
    m_uranus = 86.8099e24
    m_neptune = 102.4092e24

    planet_masses = np.array([m_mercury, m_venus, m_earth, m_mars, m_jupiter, m_saturn, m_uranus, m_neptune]) * u.kg # kg

    # Calculating in terms of M_sun
    planet_qs = (planet_masses / c.M_sun).to(u.dimensionless_unscaled).value # dimensionless

    # Defining shears
    planet_shears = []
    for i in range(len(planets)):
        q = planet_qs[i]
        s = planet_ss[i]
        if s < 1:
            planet_shears.append(q * s**2)
        else:
            planet_shears.append(q * (1/s)**2)

    planet_shears = np.array(planet_shears)

    # Defining planet angles
    seed = 2
    np.random.seed(seed)
    planet_alphas = np.random.uniform(0, 360, len(planets))

    # Adding all attributes to a Pandas DataFrame
    planet_attributes = pd.DataFrame({
        'Planet': planets,
        'Semimajor Axis [au]': planet_semimajor_axes,
        's': planet_ss,
        'Mass [kg]': planet_masses,
        'q': planet_qs,
        'Shear': planet_shears,
        'Angle [deg]': planet_alphas
    })

    # Sorting planets by shear
    planet_attributes = planet_attributes.sort_values(by='Shear', ascending=False)

    #%% Importing simulation data
    sims = []

    for i in range(len(planets)):
        print(f'Processing {i+1} planets...')
        current_planets = planet_attributes['Planet'].iloc[0:i+1]
        print(current_planets)

        file_name = current_planets.iloc[0]
        for i, planet in enumerate(current_planets.values):
            if i != 0:
                file_name += f'_{planet}'
        file_name += '.pkl'

        file_path = f'./Unity/Simulations/Part_6_Solar_System_Collection/{file_name}'

        print(f'Simulation file path: {file_path}')

        sim_current_lens = IRSC.caustic_reader(file_path)
        # sim_current_lens.plot()
        # plt.savefig(f'./Unity Analysis/Part 6/Figures/solar_system_magnification_map_{file_name}.png', dpi=300)

        sims.append(sim_current_lens)

    #%% Defining source profile
    print('Defining source profile...')
    source_radius_solar = radius * c.R_sun.to(u.m).value # meters
    source_radius_angular = source_radius_solar / D_l.to(u.m).value * u.rad # radians
    print(f'Source radius: {source_radius_angular.to(u.mas)}')
    source_radius = source_radius_angular.to(u.mas).value / theta_ein.to(u.mas).value # theta_E
    print(f'Source radius: {source_radius} theta_E')

    source_profile = IRSF.IRSFunctions.source_profile(ang_res=sims[0].ang_res, profile_type='LD', rad=source_radius, LD=0.5)

    # Convolving magnification maps with source profile
    convolved_brightnesses = []
    for sim in sims:
        convolved_brightnesses.append(sim.convolve(source_profile=source_profile))

    #%% Plotting fractional difference maps
    for i in range(len(sims)-1):
        print(f'Processing {i+2} - {i+1} planets...')

        current_sim = sims[i]
        next_sim = sims[i+1]

        ang_width = current_sim.ang_width
        ang_res = current_sim.ang_res

        frac_diff_map = (convolved_brightnesses[i+1] - convolved_brightnesses[i]) / convolved_brightnesses[i]

        fig = plt.figure()
        fig.set_constrained_layout(True)
        ax = fig.add_subplot()

        fig.suptitle(f'Fractional Difference Map ({i+2} - {i+1} planets)')
        ax.set_title(f'Source Radius: {radius:.2e} $R_\odot$')

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

        plt.savefig(f'./Unity Analysis/Part 6/Figures/frac_diff_solar_sys_{i+2}_{i+1}_{radius:.2e}.png', dpi=300)

        print(f'Fractional difference map saved successfully for {i+2} - {i+1} planets.')