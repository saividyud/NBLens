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
    from IRSMicroLensing import IRSCausticsGPU as IRSC
    print(f'Custom library import time (GPU): {(t.time() - init_time):.3} seconds')

    # ── Parameter grid (must match the bash script) ──
    SS1 = ['0.2', '0.6', '1.0', '1/0.6', '1/0.2']
    SS2 = ['0.2', '0.6', '1.0', '1/0.6', '1/0.2']
    ALPHAS = [0.0, 45.0, 90.0, 135.0, 180.0]
    TOTAL_SIMS = len(SS1) * len(SS2) * len(ALPHAS)  # 125

    def task_id_to_params(task_id):
        """Convert a linear task ID (0-124) to (s1_str, s2_str, alpha)."""
        s1_idx = task_id // 25
        s2_idx = (task_id % 25) // 5
        alpha_idx = task_id % 5
        return SS1[s1_idx], SS2[s2_idx], ALPHAS[alpha_idx]

    # ── Argument parsing ──
    parser = argparse.ArgumentParser(
        description='Compute triple lens magnification map offsetted by binary lens effective center of magnification.'
    )

    # Single-simulation mode (original interface)
    parser.add_argument('-s1', '--sep1', help='Separation of big planet')
    parser.add_argument('-s2', '--sep2', help='Separation of small planet')
    parser.add_argument('-alpha', '--alpha', help='Angle of small planet')

    # Batch mode: run multiple simulations from the 5x5x5 grid
    parser.add_argument('--batch-start', type=int, default=None,
                        help='First task ID in batch (inclusive, 0-124)')
    parser.add_argument('--batch-end', type=int, default=None,
                        help='Last task ID in batch (exclusive, 1-125)')

    args = parser.parse_args()

    # ── Build list of simulations to run ──
    if args.batch_start is not None:
        start = args.batch_start
        end = args.batch_end if args.batch_end is not None else start + 1
        end = min(end, TOTAL_SIMS)
        sim_params = [task_id_to_params(tid) for tid in range(start, end)]
        print(f'\nBatch mode: task IDs {start} to {end - 1} ({len(sim_params)} simulations)')
    elif args.sep1 is not None:
        sim_params = [(args.sep1, args.sep2, np.float64(numexpr.evaluate(args.alpha).item()))]
        print(f'\nSingle mode: s1={args.sep1}, s2={args.sep2}, alpha={sim_params[0][2]}')
    else:
        parser.error('Provide either -s1/-s2/-alpha or --batch-start/--batch-end')

    # ── Pre-load the attribute CSV once ──
    file = pd.read_csv('./Unity Analysis/Part 5/Data Files/part_5_triple_attributes.csv')
    file_directory = './Unity/Simulations/Part_5_Triple_Collection/'
    if not os.path.exists(file_directory):
        os.mkdir(file_directory)

    # ── Run simulations ──
    batch_start_time = t.time()

    for sim_idx, (s1_str, s2_str, alpha) in enumerate(sim_params):
        print()
        print('=========================================================')
        print(f'Simulation {sim_idx + 1}/{len(sim_params)}')
        print('=========================================================')

        s1 = np.float64(numexpr.evaluate(s1_str).item())
        s2 = np.float64(numexpr.evaluate(s2_str).item())
        if not isinstance(alpha, float):
            alpha = np.float64(numexpr.evaluate(str(alpha)).item())
        q1 = 1e-3
        q2 = 3e-4

        print(f'q1 = {q1}')
        print(f'q2 = {q2}')
        print(f's1 = {s1} ({s1_str})')
        print(f's2 = {s2} ({s2_str})')
        print(f'alpha = {alpha}')

        triple_lens_attributes = np.array([
            [0, 0, 1],
            [s1, 0, q1],
            [s2*np.cos(np.deg2rad(alpha)), s2*np.sin(np.deg2rad(alpha)), q2],
        ])

        center_of_magnification = np.array([q1 / ((1 + q1) * (s1 + 1/s1)), 0])
        print(f'Binary offset: {center_of_magnification}')

        triple_lens_attributes[:, :2] -= center_of_magnification

        row = file[(file['s1'] == s1_str) & (file['s2'] == s2_str) & (file['alpha'] == alpha)].iloc[0]
        row = np.array(row)

        print(row)

        pixels = int(row[3])
        ang_width = row[4]
        num_r = int(row[5])
        num_theta = int(row[6])
        y_plus = row[7]
        y_minus = row[8]

        print(f'Pixels: {pixels}')
        print(f'Angular width: {ang_width}')
        print(f'Annulus lower bound: {y_minus}')
        print(f'Annulus upper bound: {y_plus}')
        print(f'Thickness: {y_plus - y_minus}')
        print(f'Number of rays in theta: {num_theta}')
        print(f'Number of rays in r: {num_r}')

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

        file_name = f'triple_{s1:.2e}_{s2:.2e}_{alpha:.0f}.pkl'
        file_path = file_directory + file_name
        print(f'Simulation file path: {file_path}')

        if os.path.exists(file_path):
            print('Simulation already exists. Skipping...')
            continue

        print(f'Shooting triple lens:\n')

        calculator = IRSC.IRSCaustics(annulus_param_dict=triple_lens_parameters)
        magnifications = calculator.series_calculate(cm_offset='auto', annulus_offset=center_of_magnification)

        print('=========================================================')

        init_time = t.time()
        with open(file_path, 'wb') as calculator_file:
            pickle.dump(calculator, calculator_file)

        print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')
        print(f'Simulation {sim_idx + 1}/{len(sim_params)} done')

    total_time = t.time() - batch_start_time
    print('\n=========================================================')
    print(f'All {len(sim_params)} simulations completed in {total_time:.1f}s')
    print(f'Average per simulation: {total_time / len(sim_params):.1f}s')
    print('=========================================================')
