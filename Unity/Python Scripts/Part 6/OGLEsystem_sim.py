if __name__ == '__main__':
    import pickle
    import time as t
    import sys
    import argparse
    import numexpr
    sys.path.append('.')

    import platform
    import psutil
    import os

    print(platform.machine())
    print(platform.version())
    print(platform.system())
    print(platform.processor())
    print(platform.node())

    memory = psutil.virtual_memory()

    print(memory.total / (1024 ** 3))
    print(memory.used / (1024 ** 3))
    print(memory.available / (1024 ** 3))

    import numpy as np
    import pandas as pd

    import matplotlib.pyplot as plt
    from matplotlib import patches

    import astropy.units as u
    import astropy.constants as c

    from IRSMicroLensing import IRSFunctions as IRSF

    init_time = t.time()
    from IRSMicroLensing import IRSCausticsGPU as IRSC
    print(f'Custom library import time (GPU): {(t.time() - init_time):.3} seconds')

    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--job_id', type=int, required=True)
    args = parser.parse_args()
    job_id = args.job_id

    # Ensure folder exists
    file_directory = './Unity/Simulations/Part_6_OGLE_System_Collection/'
    os.makedirs(file_directory, exist_ok=True)

    # ── Run simulations ──
    batch_start_time = t.time()

    # Distance to source plane
    D_s = 8 * u.kpc

    # Distance to lens plane
    D_l = 1.5 * u.kpc

    # Mass of star
    M_star = 0.5 * c.M_sun

    theta_ein, r_ein = IRSF.IRSFunctions.e_ring([D_s, D_l], M_star.to(u.M_sun).value)
    print(f'Einstein ring radius: {r_ein}')
    print(f'Einstein ring angle: {theta_ein}')

    planets = ['False Mercury', 'False Venus', 'False Earth', 'False Mars', 'b', 'c', 'False Uranus', 'False Neptune']

    # Defining our Jupiter for scale
    r_jupiter = 5.2

    # Defining semimajor axes of known planets
    r_b = 2.3
    r_c = 4.6

    # Defining semimajor axes of false Solar System-like planets
    r_false_uranus = 19.2 * (r_b / r_jupiter)
    r_false_neptune = 30.07 * (r_b / r_jupiter)
    r_false_mercury = 0.387 * (r_b / r_jupiter)
    r_false_venus = 0.723 * (r_b / r_jupiter)
    r_false_earth = 1 * (r_b / r_jupiter)
    r_false_mars = 1.524 * (r_b / r_jupiter)

    planet_semimajor_axes = np.array([r_false_mercury, r_false_venus, r_false_earth, r_false_mars, r_b, r_c, r_false_uranus, r_false_neptune]) * u.au # au

    # Calculating in terms of theta_E
    planet_ss = (planet_semimajor_axes / r_ein).to(u.dimensionless_unscaled).value # theta_E

    # Defining our Jupiter mass for scale
    m_jupiter = 1898.125e24

    # Defining masses of known planets
    m_b = 0.71 * c.M_jup.to(u.kg).value
    m_c = 0.27 * c.M_jup.to(u.kg).value

    # Defining masses of false Solar System-like planets
    m_false_uranus = 86.8099e24 * (m_b / m_jupiter)
    m_false_neptune = 102.4092e24 * (m_b / m_jupiter)
    m_false_mercury = 0.330103e24 * (m_b / m_jupiter)
    m_false_venus = 4.86731e24 * (m_b / m_jupiter)
    m_false_earth = 5.97217e24 * (m_b / m_jupiter)
    m_false_mars = 0.641691e24 * (m_b / m_jupiter)

    planet_masses = np.array([m_false_mercury, m_false_venus, m_false_earth, m_false_mars, m_b, m_c, m_false_uranus, m_false_neptune]) * u.kg # kg

    # Calculating in terms of M_sun
    planet_qs = (planet_masses / M_star).to(u.dimensionless_unscaled).value # dimensionless

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
    seed = 3
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
    print(f'Planet attributes: {planet_attributes}')

    # Defining lens attributes
    lens_att = []
    lens_att.append([0, 0, 1])

    for i in range(len(planets)):
        lens_att.append([planet_attributes['s'].iloc[i] * np.cos(np.deg2rad(planet_attributes['Angle [deg]'].iloc[i])), planet_attributes['s'].iloc[i] * np.sin(np.deg2rad(planet_attributes['Angle [deg]'].iloc[i])), planet_attributes['q'].iloc[i]])

    lens_att = np.array(lens_att)
    
    print(f'Lens attributes: {lens_att}')

    if job_id == 0:
        fig = plt.figure()
        ax = fig.add_subplot()

        fig.suptitle('OGLE System Lens Configuration', y=0.95)

        for i in range(lens_att.shape[0]):
            ax.scatter(lens_att[i, 0], lens_att[i, 1], label=planet_attributes['Planet'].iloc[i-1] if i != 0 else 'Sun', marker='*' if i == 0 else 'o')

            ring_patch = patches.Circle((0, 0), radius=1, fill=None, edgecolor='black', linestyle='--', label='Einstein Ring')
            ax.add_patch(ring_patch)

            ax.set_aspect('equal')

            ax.set_xlabel(r'X [$\theta_E$]')
            ax.set_ylabel(r'Y [$\theta_E$]')

            view = 8
            ax.set_xlim(-view, view)
            ax.set_ylim(-view, view)
            ax.hlines(0, -view, view, colors='black', linestyles='dashed', linewidth=0.1)
            ax.vlines(0, -view, view, colors='black', linestyles='dashed', linewidth=0.1)

            ax.legend(ncol=2, loc='best')

            fig.savefig(f'./Unity Analysis/Part 6/Figures/OGLE_system_lens_configuration.png', dpi=300)

    # Creating magnification maps
    center_of_magnification = np.array([planet_attributes['q'].iloc[0] / ((1 + planet_attributes['q'].iloc[0]) * (planet_attributes['s'].iloc[0] + 1/planet_attributes['s'].iloc[0])), 0])
    print(f'Binary offset: {center_of_magnification}')

    lens_att[:, :2] -= center_of_magnification

    # Iterating over all planets and consecutively creating magnification maps
    current_lens_att = lens_att[0:job_id+2, :]
    current_planets = planet_attributes['Planet'].iloc[0:job_id+1]
    print(current_planets)

    pixels, ang_width, (qs, ss), cusp_points, (max_source_radius, min_source_radius) = IRSF.IRSFunctions._ang_width_calculator(current_lens_att, final_multiplier=1, pixels_in_small_source=20, cm_offset='auto')
    print(f'Pixels = {pixels}, Angular Width = {ang_width}, Max Source Radius = {max_source_radius}, Min Source Radius = {min_source_radius}')

    (y_plus, y_minus), min_mag = IRSF.IRSFunctions._annulus_bounds_calculator(ang_width, qs, ss)
    print(f'Annulus bounds: y+ = {y_plus}, y- = {y_minus}, Min Magnification in Annulus = {min_mag}')

    num_r, num_theta = IRSF.IRSFunctions._num_ray_calculator(pixels, ang_width, y_plus, y_minus, min_mag, delta=0.01, r_theta_ratio=4)
    print(f'Number of rays: num_r = {num_r}, num_theta = {num_theta}, Total number of rays: {(num_r * num_theta):.3e}')

    lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': current_lens_att.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    sim_current_lens = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_current_lens = sim_current_lens.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)

    file_name = current_planets.iloc[0]
    for i, planet in enumerate(current_planets.values):
        if i != 0:
            file_name += f'_{planet}'
    file_name += '.pkl'

    file_path = f'{file_directory}/{file_name}'

    print(f'Simulation file path: {file_path}')

    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(sim_current_lens, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

    print(f'Total simulation time: {(t.time() - batch_start_time):.3} seconds')