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

    from IRSMicroLensing import IRSFunctions as IRSF

    import multiprocessing as mp
    mp.set_start_method('spawn')

    init_time = t.time()
    from IRSMicroLensing import IRSCausticsGPU as IRSC
    print(f'Custom library import time (GPU): {(t.time() - init_time):.3} seconds')

    # Initialize parser
    parser = argparse.ArgumentParser(description='Compute binary lens magnification map offsetted by binary lens effective center of magnification.')

    # Adding arguments
    parser.add_argument('-s1', '--sep1', help='Separation of big planet')
    parser.add_argument('-s2', '--sep2', help='Separation of small planet')
    parser.add_argument('-q2', '--massratio2', help='Mass ratio of small planet')
    parser.add_argument('-alpha', '--alpha', help='Angle of small planet')

    args = vars(parser.parse_args())
    print(args)

    print()
    print('=========================================================')
    print('Lens parameters:')

    ''' Preparing lens parameters '''
    # Big planet parameters
    s1_str = args['sep1']
    s1 = np.float64(numexpr.evaluate(s1_str).item())
    s2_str = args['sep2']
    s2 = np.float64(numexpr.evaluate(s2_str).item())
    alpha = np.float64(numexpr.evaluate(args['alpha']).item())
    q1 = 1e-3
    q2 = np.float64(numexpr.evaluate(args['massratio2']).item())

    print(f'q1 = {q1}')
    print(f'q2 = {q2}')
    print(f's1 = {s1}')
    print(f's2 = {s2}')
    print(f'alpha = {alpha}')

    # Defining triple lens attributes
    triple_lens_attributes = np.array([
        [0, 0, 1],
        [s1, 0, q1],
        [s2*np.cos(np.deg2rad(alpha)), s2*np.sin(np.deg2rad(alpha)), q2],
    ])

    # Binary offset
    center_of_magnification = np.array([q1 / ((1 + q1) * (s1 + 1/s1)), 0])
    print(f'Binary offset: {center_of_magnification}')

    triple_lens_attributes[:, :2] -= center_of_magnification

    pixels, ang_width, (qs_list_1, ss_list_1), cusp_points, (max_source_radius, min_source_radius) = IRSF.IRSFunctions._ang_width_calculator(triple_lens_attributes, final_multiplier=1, pixels_in_small_source=20, cm_offset='auto')

    # Setting pixels and angular width to 5000 and 0.4 respectively
    pixels = 5000
    ang_width = 0.4

    (y_plus, y_minus), min_mag = IRSF.IRSFunctions._annulus_bounds_calculator(ang_width, qs=qs_list_1, ss=ss_list_1)
    num_r, num_theta = IRSF.IRSFunctions._num_ray_calculator(pixels, ang_width, y_plus, y_minus, min_mag, delta=0.01, r_theta_ratio=4)

    print(f'Pixels: {pixels}')
    print(f'Angular width: {ang_width}')
    print(f'Annulus lower bound: {y_minus}')
    print(f'Annulus upper bound: {y_plus}')
    print(f'Thickness: {y_plus - y_minus}')
    print(f'Number of rays in theta: {num_theta}')
    print(f'Number of rays in r: {num_r}')
    print(f'Total number of rays: {(num_r * num_theta):.3e}')
    print(f'Maximum source radius: {max_source_radius}')
    print(f'Minimum source radius: {min_source_radius}')

    triple_lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': triple_lens_attributes.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    print(f'Number of rays: {(num_r * num_theta):.4e}')
    print('=========================================================')

    ''' Simulating L lens magnification map '''
    file_directory = f'./Unity/Simulations/Part_7_Triple_Collection/'

    param_dict = triple_lens_parameters
    file_name = f'triple_{s2:.2e}_{q2:.2e}_{alpha:.0f}.pkl'

    file_path = file_directory + file_name

    print(f'Simulation file path: {file_path}')

    # Checking if file directory exists
    if not os.path.exists(file_directory):
        os.mkdir(file_directory)

    # If simulation already exists, skipping
    if os.path.exists(file_path):
        print('Simulation already exists. Skipping...')
        exit()

    print(f'Shooting triple lens:\n')

    calculator = IRSC.IRSCaustics(annulus_param_dict=param_dict)
    magnifications = calculator.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)
    # magnifications = calculator.parallel_calculate(cm_offset='auto', annulus_offset=center_of_magnification, rows=1)

    print('=========================================================')

    ''' Saving class data to file '''
    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(calculator, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

    print('Done')
