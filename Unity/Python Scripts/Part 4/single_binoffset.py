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

    # Initialize parser
    parser = argparse.ArgumentParser(description='Compute single lens magnification map.')

    # Adding arguments
    parser.add_argument('-s', '--sep', help='Seperation of big planet')
    parser.add_argument('-q', '--star_mass_ratio', help='Big planet / star mass ratio')

    args = vars(parser.parse_args())
    print(args)

    print()
    print('=========================================================')
    print('Lens parameters:')

    ''' Preparing lens parameters '''
    # Big planet parameters
    s_str = args['sep']
    s = np.float64(numexpr.evaluate(s_str).item())
    q = np.float64(numexpr.evaluate(args['star_mass_ratio']).item())

    print(f'q1 = {q}')
    print(f's1 = {s}')    

    # Defining single lens attributes
    single_lens_attributes = np.array([
        [0, 0, 1]
    ])

    # Reading in pre-calculated values
    file = pd.read_csv('./Unity Analysis/Part 4/Data Files/part_4_single_attributes.csv')

    row = file[(file['s'] == s) & (file['q'] == q)].iloc[0]
    row = np.array(row)

    print(row)

    pixels = int(row[2])
    ang_width = row[3]
    num_r = int(row[4])
    num_theta = int(row[5])
    y_plus = row[6]
    y_minus = row[7]

    print(f'Pixels: {pixels}')
    print(f'Angular width: {ang_width}')
    print(f'Annulus lower bound: {y_minus}')
    print(f'Annulus upper bound: {y_plus}')
    print(f'Thickness: {y_plus - y_minus}')
    print(f'Number of rays in theta: {num_theta}')
    print(f'Number of rays in r: {num_r}')

    single_lens_parameters = {
        'pixels': pixels,
        'ang_width': ang_width,
        'thickness': y_plus - y_minus,
        'y_plus': y_plus,
        'y_minus': y_minus,
        'lens_att': single_lens_attributes.tolist(),
        'num_theta': num_theta,
        'num_r': num_r
    }

    print(f'Number of rays: {(num_r * num_theta):.4e}')
    print('=========================================================')

    ''' Simulating single lens magnification map '''
    file_directory = f'./Unity/Simulations/Part_4_Single_Collection/'

    param_dict = single_lens_parameters
    file_name = f'single_{q:.0e}_{s:.2e}.pkl'

    file_path = file_directory + file_name

    print(f'Simulation file path: {file_path}')

    # Checking if file directory exists
    if not os.path.exists(file_directory):
        os.mkdir(file_directory)

    print(f'Shooting single lens:\n')

    calculator = IRSC.IRSCaustics(annulus_param_dict=param_dict)
    magnifications = calculator.series_calculate(cm_offset='auto', rows=1)
    # magnifications = calculator.parallel_calculate(cm_offset='auto', rows=1)

    print('=========================================================')

    ''' Saving class data to file '''
    init_time = t.time()
    with open(file_path, 'wb') as calculator_file:
        pickle.dump(calculator, calculator_file)

    print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')

    print('Done')
