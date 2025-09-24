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

    # import multiprocessing as mp
    # mp.set_start_method('spawn')

    init_time = t.time()
    from IRSMicroLensing import IRSCaustics as IRSC
    from IRSMicroLensing import IRSFunctions as IRSF
    print(f'Custom library import time: {(t.time() - init_time):.3} seconds')

    print()
    print('=========================================================')
    print('Lens parameters:')

    ''' Preparing lens parameters '''
    # Big planet parameters
    s = 1/0.7
    q = 1e-3

    print(f'q1 = {q}')
    print(f's1 = {s}')

    # Defining binary lens attributes
    lens_attributes = np.array([
        [0, 0, 1]
    ])

    # Binary offset
    binary_offset = np.array([q / ((1 + q) * (s + 1/s)), 0])
    print(f'Binary offset: {binary_offset}')

    # Correcting binary lens attributes
    # lens_attributes[:, :2] -= binary_offset

    # Map parameters
    pixels = 4588
    delta = 0.01

    ang_width = 0.345727450083596
    num_r = 1700804
    num_theta = 425201

    y_plus = 1.136547316616518
    y_minus = 0.8798577810002515
    thickness = y_plus - y_minus

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
        'lens_att': lens_attributes.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    print(f'Number of rays: {(num_r * num_theta):.4e}')
    print('=========================================================')

    ''' Simulating L lens magnification map '''
    file_directory = f'./Unity/Simulations/Tests/'

    param_dict = binary_lens_parameters
    file_name = f'source_size_deviations_single_3.pkl'

    file_path = file_directory + file_name

    print(f'Simulation file path: {file_path}')

    # Checking if file directory exists
    if not os.path.exists(file_directory):
        raise FileNotFoundError(f'File directory {file_directory} does not exist.')

    print(f'Shooting single lens:\n')

    calculator = IRSC.IRSCaustics(annulus_param_dict=param_dict)
    magnifications = calculator.series_calculate(cm_offset='auto', rows=5)
    # magnifications = calculator.parallel_calculate(cm_offset='auto', cpus=6, rows=10)

    print('=========================================================')

    ''' Saving class data to file '''
    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(calculator, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

    print('Done')
