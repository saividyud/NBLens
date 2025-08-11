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

    import multiprocessing as mp
    mp.set_start_method('spawn')

    init_time = t.time()
    from IRSMicroLensing import IRSCaustics as IRSC
    from IRSMicroLensing import IRSFunctions as IRSF
    print(f'Custom library import time: {(t.time() - init_time):.3} seconds')

    # Initialize parser
    parser = argparse.ArgumentParser(description='Compute binary lens magnification map offsetted by binary lens effective center of magnification.')

    # Adding arguments
    parser.add_argument('-s1', '--sep1', help='Seperation of big planet')
    parser.add_argument('-q1', '--star_mass_ratio', help='Big planet / star mass ratio')

    args = vars(parser.parse_args())
    print(args)

    print()
    print('=========================================================')
    print('Lens parameters:')

    ''' Preparing lens parameters '''
    # Big planet parameters
    s1 = np.float64(numexpr.evaluate(args['sep1']).item())
    alpha1 = 0
    q1 = np.float64(numexpr.evaluate(args['star_mass_ratio']).item())

    print(f'q1 = {q1}')
    print(f's1 = {s1}')
    print(f'alpha1 = {alpha1}')

    # Defining binary lens attributes
    binary_lens_attributes = np.array([
        [0, 0, 1],
        [s1*np.cos(np.deg2rad(alpha1)), s1*np.sin(np.deg2rad(alpha1)), q1],
    ])

    # Rotation matrix for first planet
    first_planet_DCM = np.array([
        [np.cos(np.deg2rad(alpha1)), -np.sin(np.deg2rad(alpha1))],
        [np.sin(np.deg2rad(alpha1)), np.cos(np.deg2rad(alpha1))]
    ])

    # Binary offset
    delta1_rot = np.array([q1 / ((1 + q1) * (s1 + 1/s1)), 0]).reshape(-1, 1)
    delta1 = np.dot(first_planet_DCM, delta1_rot).reshape(2)
    print(f'Binary offset: {delta1}')

    # Correcting binary lens attributes
    binary_lens_attributes[:, :2] -= delta1

    # Map parameters
    pixels = 2000
    delta = 0.01

    # Reading in pre-calculated values
    file = pd.read_csv('./Unity Analysis/binary_attributes.csv')

    row = file[(file['s1'] == s1) & (file['q1'] == q1)].iloc[0]
    row = np.array(row)

    s1 = row[0]
    q1 = row[1]
    ang_width = row[2]
    thickness = row[3]
    y_plus = row[4]
    y_minus = row[5]
    num_r = int(row[6])
    num_theta = int(row[7])

    print(f'Angular width: {ang_width}')
    print(f'Thickness: {thickness}')
    print(f'Annulus lower bound: {y_minus}')
    print(f'Annulus upper bound: {y_plus}')
    print(f'Number of rays in theta: {num_theta}')
    print(f'Number of rays in r: {num_r}')

    binary_lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': thickness,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': binary_lens_attributes.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    print(f'Number of rays: {(num_r * num_theta):.4e}')
    print('=========================================================')

    ''' Simulating L lens magnification map '''
    file_directory = f'./Unity/Simulations/Binary_Collection/'

    param_dict = binary_lens_parameters
    file_name = f'binary_{q1:.0e}_{s1:.2e}.pkl'

    file_path = file_directory + file_name

    print(f'Simulation file path: {file_path}')

    # Checking if file directory exists
    if not os.path.exists(file_directory):
        raise FileNotFoundError(f'File directory {file_directory} does not exist.')

    print(f'Shooting binary lens:\n')

    calculator = IRSC.IRSCaustics(annulus_param_dict=param_dict)
    magnifications = calculator.parallel_calculate(cm_offset='auto', cpus=10)
    # magnifications = calculator.series_calculate(cm_offset='auto')

    print('=========================================================')

    ''' Saving class data to file '''
    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(calculator, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

    print('Done')
