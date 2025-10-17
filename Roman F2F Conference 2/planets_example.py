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

    import multiprocessing as mp
    mp.set_start_method('spawn')

    init_time = t.time()
    from IRSMicroLensing import IRSCaustics as IRSC
    from IRSMicroLensing import IRSFunctions as IRSF
    print(f'Custom library import time: {(t.time() - init_time):.3} seconds')

    def caustic_writer(simulation, file_directory, file_name):
        # Checking if file directory exists
        if not os.path.exists(file_directory):
            raise FileNotFoundError(f'File directory {file_directory} does not exist.')

        ''' Saves class data to file '''
        init_time = t.time()
        file_path = file_directory + file_name

        print(f'Simulation file path: {file_path}')

        with open(file_path, 'wb') as calculator_file:
            pickle.dump(simulation, calculator_file)

        print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

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

    print(f's_jup = {s_jup:.3f}')
    print(f's_saturn = {s_saturn:.3f}')
    print(f's_earth = {s_earth:.3f}')
    print(f's_uranus = {s_uranus:.3f}')

    print(f'q_jup = {q_jup:.3e}')
    print(f'q_saturn = {q_saturn:.3e}')
    print(f'q_earth = {q_earth:.3e}')
    print(f'q_uranus = {q_uranus:.3e}')

    pixels = 2304
    ang_width = 0.33502742168597877

    #--------------------------------------------------------------------------
    sun_jupiter_lens_att = np.array([
        [0, 0, 1],
        [s_jup, 0, q_jup]
    ])

    y_plus = 1.1316152394421348
    y_minus = 0.8828576046682406
    num_r = 856428
    num_theta = 214107

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

    sim_sun_jupiter = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_sun_jupiter = sim_sun_jupiter.parallel_calculate(cm_offset='auto', rows=5)

    caustic_writer(sim_sun_jupiter, './Roman F2F Conference 2/', 'sun_jupiter_example.pkl')
    
    #--------------------------------------------------------------------------
    alpha_saturn = np.deg2rad(-45)
    s_saturn = 0.8

    sun_jupiter_saturn_lens_att = np.array([
        [0, 0, 1],
        [s_jup, 0, q_jup],
        [s_saturn * np.cos(alpha_saturn), s_saturn * np.sin(alpha_saturn), q_saturn]
    ])

    y_plus = 1.132076972940417
    y_minus = 0.8828576046682406
    num_r = 872580
    num_theta = 218145

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

    sim_sun_jupiter_saturn = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_sun_jupiter_saturn = sim_sun_jupiter_saturn.parallel_calculate(cm_offset='auto', rows=5)

    caustic_writer(sim_sun_jupiter_saturn, './Roman F2F Conference 2/', 'sun_jupiter_saturn_example.pkl')

    #--------------------------------------------------------------------------
    alpha_uranus = -120

    four_lens_att = np.array([
        [0, 0, 1],
        [s_jup, 0, q_jup],
        [s_saturn * np.cos(alpha_saturn), s_saturn * np.sin(alpha_saturn), q_saturn],
        [s_uranus * np.cos(np.deg2rad(alpha_uranus)), s_uranus * np.sin(np.deg2rad(alpha_uranus)), q_uranus]
    ])

    y_plus = 1.132076972940417
    y_minus = 0.8828576046682406
    num_r = 872580
    num_theta = 218145

    lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': four_lens_att.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    sim_four = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_four = sim_four.parallel_calculate(cm_offset='auto', rows=5)

    caustic_writer(sim_four, './Roman F2F Conference 2/', 'four_lens_example.pkl')

    #--------------------------------------------------------------------------
    alpha_earth = 30

    five_lens_att = np.array([
        [0, 0, 1],
        [s_jup, 0, q_jup],
        [s_saturn * np.cos(alpha_saturn), s_saturn * np.sin(alpha_saturn), q_saturn],
        [s_earth * np.cos(np.deg2rad(alpha_earth)), s_earth * np.sin(np.deg2rad(alpha_earth)), q_earth],
        [s_uranus * np.cos(np.deg2rad(alpha_uranus)), s_uranus * np.sin(np.deg2rad(alpha_uranus)), q_uranus]
    ])

    y_plus = 1.132077174511849
    y_minus = 0.8828576046682406
    num_r = 872580
    num_theta = 218145

    lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': five_lens_att.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    sim_five = IRSC.IRSCaustics(annulus_param_dict=lens_parameters)
    magnifications_five = sim_five.parallel_calculate(cm_offset='auto', rows=5)

    caustic_writer(sim_five, './Roman F2F Conference 2/', 'five_lens_example.pkl')
