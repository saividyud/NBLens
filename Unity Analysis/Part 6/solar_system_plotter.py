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

    def half_decade_levels(data, num_steps=6, min_steps=3, keep_zero=True):
        '''Build symmetric half-decade contour levels adapted to the data range.

        Levels follow the base ``[1, 3]`` pattern (..., 1e-2, 3e-2, 1e-1, 3e-1,
        1, 3, ...) so each decade contributes two steps. The largest level is
        anchored to the biggest such step at or below the data's maximum absolute
        value, guaranteeing the outer contours are actually reached. ``num_steps``
        levels (>= ``min_steps``) are then placed per sign, with a 0 contour.
        '''
        num_steps = max(num_steps, min_steps)

        base = (1.0, 3.0)

        def value_at(k):
            # Maps integer index k -> base[1, 3] level, e.g. k=0 -> 1, k=1 -> 3,
            # k=2 -> 10, k=-1 -> 0.3, k=-2 -> 0.1 (floor division handles k<0).
            return base[k % 2] * 10.0 ** (k // 2)

        finite = data[np.isfinite(data)]
        max_abs = np.max(np.abs(finite)) if finite.size else 0.0

        if not np.isfinite(max_abs) or max_abs <= 0:
            # Flat / degenerate map: fall back to a sane default range (top -> 0.3).
            top = -1
        else:
            # Largest [1, 3] step at/below the data max (overshoot, then back off).
            top = int(np.ceil(2 * np.log10(max_abs)))
            while value_at(top) > max_abs:
                top -= 1

        # Descending indices top, top-1, ... -> ascending positive levels.
        pos = np.sort([value_at(top - j) for j in range(num_steps)])
        neg = -pos[::-1]

        if keep_zero:
            return np.concatenate([neg, [0.0], pos])
        else:
            return np.concatenate([neg, pos])

    def linewidths_by_level(levels, min_lw=0.3, max_lw=1.3):
        '''Thicker lines for larger-magnitude contours, symmetric about 0.'''
        pos_levels = sorted({abs(v) for v in levels if v != 0})
        n = len(pos_levels)
        if n > 1:
            lw_map = {lv: min_lw + (max_lw - min_lw) * idx / (n - 1)
                      for idx, lv in enumerate(pos_levels)}
        elif n == 1:
            lw_map = {pos_levels[0]: max_lw}
        else:
            lw_map = {}

        lws = []
        for v in levels:
            if v == 0:
                lws.append((min_lw + max_lw) / 2)
            else:
                lws.append(lw_map[abs(v)])
        return lws

    #%% Defining simulation parameters
    # Distance to source plane
    D_s = 8 * u.kpc

    # Distance to lens plane
    D_l = 4 * u.kpc

    # Mass of star
    M_star = 1 * c.M_sun

    theta_ein, r_ein = IRSF.IRSFunctions.e_ring([D_s, D_l], M_star.to(u.M_sun).value)
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
        # plt.savefig(f'./Unity Analysis/Part 6/Figures/Solar System/solar_system_magnification_map_{file_name}.png', dpi=300)

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
        # Adapt contour levels to this map so lines are always visible, even
        # when subtracting small planets shrinks the fractional deviations.
        levels = half_decade_levels(frac_diff_map, num_steps=3, min_steps=3, keep_zero=False)
        print(f'Contour levels: {levels}')

        plot = ax.contour(np.flip(frac_diff_map, axis=0),
                levels=levels,
                colors=colors_by_log(levels),
                linewidths=linewidths_by_level(levels),
                extent=[-ang_width/2 - ang_res/2, ang_width/2 + ang_res/2, -ang_width/2 - ang_res/2, ang_width/2 + ang_res/2]
        )

        bar = plt.colorbar(plot)
        bar.set_label('Fractional Difference')

        ax.set_xlabel(r'X [$\theta_E$]')
        ax.set_ylabel(r'Y [$\theta_E$]')

        view = 0.1

        ax.set_xlim(-view, view)
        ax.set_ylim(-view, view)

        plt.savefig(f'./Unity Analysis/Part 6/Figures/Solar System/frac_diff_solar_sys_{i+2}_{i+1}_{radius:.2e}.png', dpi=300)

        print(f'Fractional difference map saved successfully for {i+2} - {i+1} planets.')