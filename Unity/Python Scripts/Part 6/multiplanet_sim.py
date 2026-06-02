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

    import astropy.units as u
    import astropy.constants as c

    from IRSMicroLensing import IRSFunctions as IRSF

    # import multiprocessing as mp
    # mp.set_start_method('spawn')

    init_time = t.time()
    from IRSMicroLensing import IRSCausticsGPU as IRSC
    print(f'Custom library import time (GPU): {(t.time() - init_time):.3} seconds')

    # Ensure folder exists
    file_directory = './Unity/Simulations/Part_6_Multiplanet_Collection/'
    os.makedirs(file_directory, exist_ok=True)

    # ── Run simulations ──
    batch_start_time = t.time()

    # Distance to source plane
    D_s = 8 * u.kpc

    # Distance to lens plane
    D_l = 4 * u.kpc

    theta_ein, r_ein = IRSF.IRSFunctions.e_ring([D_s, D_l], 1)

    r_jup = 5.2 * u.au
    r_saturn = 9.5 * u.au
    r_earth = 1 * u.au
    r_uranus = 19.2 * u.au

    s_jup = (r_jup / r_ein).to(u.dimensionless_unscaled).value
    s_saturn = (r_saturn / r_ein).to(u.dimensionless_unscaled).value
    s_earth = (r_earth / r_ein).to(u.dimensionless_unscaled).value
    s_uranus = (r_uranus / r_ein).to(u.dimensionless_unscaled).value

    q_jup = (c.M_jup / c.M_sun).to(u.dimensionless_unscaled).value
    q_saturn = (5.683e26 * u.kg / c.M_sun).to(u.dimensionless_unscaled).value
    q_earth = (c.M_earth / c.M_sun).to(u.dimensionless_unscaled).value
    q_uranus = (8.681e25 * u.kg / c.M_sun).to(u.dimensionless_unscaled).value

    #==========================================================================
    # Sun-Jupiter Lens Configuration
    #==========================================================================

    alpha_jupiter = 0

    sun_jupiter_lens_att = np.array([
        [0, 0, 1],
        [s_jup * np.cos(alpha_jupiter), s_jup * np.sin(alpha_jupiter), q_jup]
    ])

    center_of_magnification = np.array([q_jup / ((1 + q_jup) * (s_jup + 1/s_jup)), 0])
    print(f'Binary offset: {center_of_magnification}')

    sun_jupiter_lens_att[:, :2] -= center_of_magnification

    pixels, ang_width, (qs, ss), cusp_points, (max_source_radius, min_source_radius) = IRSF.IRSFunctions._ang_width_calculator(sun_jupiter_lens_att, final_multiplier=1, pixels_in_small_source=20, cm_offset='auto')
    print(f'Pixels = {pixels}, Angular Width = {ang_width}, Max Source Radius = {max_source_radius}, Min Source Radius = {min_source_radius}')

    y_plus = 1.2790793218503667
    y_minus = 0.7812645092009943
    min_mag = 1.9432064097360755
    # (y_plus, y_minus), min_mag = IRSF.IRSFunctions._annulus_bounds_calculator(ang_width, qs, ss)
    # print(f'Annulus bounds: y+ = {y_plus}, y- = {y_minus}, Min Magnification in Annulus = {min_mag}')

    num_r, num_theta = IRSF.IRSFunctions._num_ray_calculator(pixels, ang_width, y_plus, y_minus, min_mag, delta=0.01, r_theta_ratio=4)
    print(f'Number of rays: num_r = {num_r}, num_theta = {num_theta}')
    print(f'Total number of rays: {(num_r * num_theta):.3e}')

    lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': sun_jupiter_lens_att.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    file_name = f'sun_jupiter_2.pkl'
    file_path = file_directory + file_name

    print(f'Simulation file path: {file_path}')

    calculator_sun_jupiter = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_sun_jupiter = calculator_sun_jupiter.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)
    
    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(calculator_sun_jupiter, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

    #==========================================================================
    # Sun-Jupiter-Saturn Lens Configuration
    #==========================================================================

    alpha_saturn = np.deg2rad(-45)

    sun_jupiter_saturn_lens_att = np.array([
        [0, 0, 1],
        [s_jup * np.cos(alpha_jupiter), s_jup * np.sin(alpha_jupiter), q_jup],
        [s_saturn * np.cos(alpha_saturn), s_saturn * np.sin(alpha_saturn), q_saturn]
    ])

    center_of_magnification = np.array([q_jup / ((1 + q_jup) * (s_jup + 1/s_jup)), 0])
    print(f'Binary offset: {center_of_magnification}')

    sun_jupiter_saturn_lens_att[:, :2] -= center_of_magnification

    pixels, ang_width, (qs, ss), cusp_points, (max_source_radius, min_source_radius) = IRSF.IRSFunctions._ang_width_calculator(sun_jupiter_saturn_lens_att, final_multiplier=1, pixels_in_small_source=20, cm_offset='auto')
    print(f'Pixels = {pixels}, Angular Width = {ang_width}, Max Source Radius = {max_source_radius}, Min Source Radius = {min_source_radius}')

    y_plus = 1.1319354138159257
    y_minus = 0.8828576046682406
    min_mag = 4.056895881009251
    # (y_plus, y_minus), min_mag = IRSF.IRSFunctions._annulus_bounds_calculator(ang_width, qs, ss)
    # print(f'Annulus bounds: y+ = {y_plus}, y- = {y_minus}, Min Magnification in Annulus = {min_mag}')

    num_r, num_theta = IRSF.IRSFunctions._num_ray_calculator(pixels, ang_width, y_plus, y_minus, min_mag, delta=0.01, r_theta_ratio=4)
    print(f'Number of rays: num_r = {num_r}, num_theta = {num_theta}')
    print(f'Total number of rays: {(num_r * num_theta):.3e}')
    
    lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': sun_jupiter_saturn_lens_att.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }
    
    file_name = f'sun_jupiter_saturn_2.pkl'
    file_path = file_directory + file_name

    print(f'Simulation file path: {file_path}')

    calculator_sun_jupiter_saturn = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_sun_jupiter_saturn = calculator_sun_jupiter_saturn.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)
    
    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(calculator_sun_jupiter_saturn, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

    #==========================================================================
    # Sun-Jupiter-Saturn-Uranus Lens Configuration
    #==========================================================================

    alpha_uranus = -120

    sun_jupiter_saturn_uranus_lens_att = np.array([
        [0, 0, 1],
        [s_jup * np.cos(alpha_jupiter), s_jup * np.sin(alpha_jupiter), q_jup],
        [s_saturn * np.cos(alpha_saturn), s_saturn * np.sin(alpha_saturn), q_saturn],
        [s_uranus * np.cos(alpha_uranus), s_uranus * np.sin(alpha_uranus), q_uranus]
    ])

    center_of_magnification = np.array([q_jup / ((1 + q_jup) * (s_jup + 1/s_jup)), 0])
    print(f'Binary offset: {center_of_magnification}')

    sun_jupiter_saturn_uranus_lens_att[:, :2] -= center_of_magnification

    pixels, ang_width, (qs, ss), cusp_points, (max_source_radius, min_source_radius) = IRSF.IRSFunctions._ang_width_calculator(sun_jupiter_saturn_uranus_lens_att, final_multiplier=1, pixels_in_small_source=20, cm_offset='auto')
    print(f'Pixels = {pixels}, Angular Width = {ang_width}, Max Source Radius = {max_source_radius}, Min Source Radius = {min_source_radius}')

    y_plus = 1.1320555131852073
    y_minus = 0.8828576046682406
    min_mag = 4.056895881009251
    # (y_plus, y_minus), min_mag = IRSF.IRSFunctions._annulus_bounds_calculator(ang_width, qs, ss)
    # print(f'Annulus bounds: y+ = {y_plus}, y- = {y_minus}, Min Magnification in Annulus = {min_mag}')

    print(f'Annulus bounds: y+ = {y_plus}, y- = {y_minus}, Min Magnification in Annulus = {min_mag}')
    
    num_r, num_theta = IRSF.IRSFunctions._num_ray_calculator(pixels, ang_width, y_plus, y_minus, min_mag, delta=0.01, r_theta_ratio=4)
    print(f'Number of rays: num_r = {num_r}, num_theta = {num_theta}')
    print(f'Total number of rays: {(num_r * num_theta):.3e}')

    lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': sun_jupiter_saturn_uranus_lens_att.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    file_name = f'sun_jupiter_saturn_uranus_2.pkl'
    file_path = file_directory + file_name

    print(f'Simulation file path: {file_path}')

    calculator_sun_jupiter_saturn_uranus = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_sun_jupiter_saturn_uranus = calculator_sun_jupiter_saturn_uranus.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)

    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(calculator_sun_jupiter_saturn_uranus, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

    #==========================================================================
    # Sun-Jupiter-Saturn-Uranus-Earth Lens Configuration
    #==========================================================================

    alpha_earth = 30

    sun_jupiter_saturn_uranus_earth_lens_att = np.array([
        [0, 0, 1],
        [s_jup * np.cos(alpha_jupiter), s_jup * np.sin(alpha_jupiter), q_jup],
        [s_saturn * np.cos(alpha_saturn), s_saturn * np.sin(alpha_saturn), q_saturn],
        [s_uranus * np.cos(alpha_uranus), s_uranus * np.sin(alpha_uranus), q_uranus],
        [s_earth * np.cos(alpha_earth), s_earth * np.sin(alpha_earth), q_earth]
    ])

    center_of_magnification = np.array([q_jup / ((1 + q_jup) * (s_jup + 1/s_jup)), 0])
    print(f'Binary offset: {center_of_magnification}')

    sun_jupiter_saturn_uranus_earth_lens_att[:, :2] -= center_of_magnification

    pixels, ang_width, (qs, ss), cusp_points, (max_source_radius, min_source_radius) = IRSF.IRSFunctions._ang_width_calculator(sun_jupiter_saturn_uranus_earth_lens_att, final_multiplier=1, pixels_in_small_source=20, cm_offset='auto')
    print(f'Pixels = {pixels}, Angular Width = {ang_width}, Max Source Radius = {max_source_radius}, Min Source Radius = {min_source_radius}')

    y_plus = 1.132077174511849
    y_minus = 0.8828576046682406
    min_mag = 4.056895881009251
    # (y_plus, y_minus), min_mag = IRSF.IRSFunctions._annulus_bounds_calculator(ang_width, qs, ss)
    # print(f'Annulus bounds: y+ = {y_plus}, y- = {y_minus}, Min Magnification in Annulus = {min_mag}')

    print(f'Annulus bounds: y+ = {y_plus}, y- = {y_minus}, Min Magnification in Annulus = {min_mag}')
    
    num_r, num_theta = IRSF.IRSFunctions._num_ray_calculator(pixels, ang_width, y_plus, y_minus, min_mag, delta=0.01, r_theta_ratio=4)
    print(f'Number of rays: num_r = {num_r}, num_theta = {num_theta}')
    print(f'Total number of rays: {(num_r * num_theta):.3e}')
    
    lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': sun_jupiter_saturn_uranus_earth_lens_att.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }
    
    file_name = f'sun_jupiter_saturn_uranus_earth_2.pkl'
    file_path = file_directory + file_name

    print(f'Simulation file path: {file_path}')

    calculator_sun_jupiter_saturn_uranus_earth = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_sun_jupiter_saturn_uranus_earth = calculator_sun_jupiter_saturn_uranus_earth.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)

    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(calculator_sun_jupiter_saturn_uranus_earth, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')