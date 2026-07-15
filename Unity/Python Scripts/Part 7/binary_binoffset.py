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

    # ── Parameter grid (1×1 = 1 simulation) ──
    # Part 7 binary reference: s1 and q1 are both fixed.
    SS1 = ['1.2']
    Q1S = [1e-3]
    TOTAL_SIMS = len(SS1) * len(Q1S)  # 1

    def task_id_to_params(task_id):
        """Convert a linear task ID (0-0) to (s1_str, q1).

        Ordering: s1 varies slowest, then q1 (fastest).
        """
        s1_idx = task_id // len(Q1S)
        q1_idx = task_id % len(Q1S)
        return SS1[s1_idx], Q1S[q1_idx]

    # ── Argument parsing ──
    parser = argparse.ArgumentParser(
        description='Adaptive binary lens magnification map calculator (GPU/CPU).'
    )

    parser.add_argument('-s1', '--sep1', help='Separation of primary lens')
    parser.add_argument('-q1', '--massratio1', help='Mass ratio of primary lens')

    parser.add_argument('--batch-start', type=int, default=None,
                        help='First task ID in batch (inclusive, 0-0)')
    parser.add_argument('--batch-end', type=int, default=None,
                        help='Last task ID in batch (exclusive, 1-1)')

    parser.add_argument('--task-id', type=int, default=None,
                        help='Run a single simulation by task ID (0-0)')

    parser.add_argument('--mode', choices=['gpu', 'cpu', 'auto'], default='auto',
                        help='Compute backend (default: auto-detect)')

    parser.add_argument('--reverse', action='store_true',
                        help='Process batch in reverse order (GPU runs 0->0 to avoid overlap with CPU)')

    args = parser.parse_args()

    # ── Determine compute mode ──
    mode = args.mode
    if mode == 'auto':
        try:
            import cupy as cp
            cp.cuda.Device(0).compute_capability
            mode = 'gpu'
            print('Auto-detected: GPU mode')
        except Exception:
            mode = 'cpu'
            print('Auto-detected: CPU mode')

    print(f'Compute mode: {mode.upper()}')

    # ── Import appropriate library ──
    if mode == 'gpu':
        init_time = t.time()
        from IRSMicroLensing import IRSCausticsGPU as IRSC
        print(f'Custom library import time (GPU): {(t.time() - init_time):.3} seconds')
    else:
        init_time = t.time()
        from IRSMicroLensing import IRSCaustics as IRSC
        print(f'Custom library import time (CPU): {(t.time() - init_time):.3} seconds')

    # ── Build list of simulations to run ──
    if args.task_id is not None:
        sim_params = [task_id_to_params(args.task_id)]
        print(f'\nSingle task mode: task ID {args.task_id}')
    elif args.batch_start is not None:
        start = args.batch_start
        end = args.batch_end if args.batch_end is not None else start + 1
        end = min(end, TOTAL_SIMS)
        sim_params = [task_id_to_params(tid) for tid in range(start, end)]
        if args.reverse:
            sim_params = list(reversed(sim_params))
        print(f'\nBatch mode: task IDs {start} to {end - 1} ({len(sim_params)} simulations)')
        if args.reverse:
            print('Processing in REVERSE order')
    elif args.sep1 is not None:
        sim_params = [(
            args.sep1,
            np.float64(numexpr.evaluate(args.massratio1).item()),
        )]
        print(f'\nSingle mode: s1={args.sep1}, q1={sim_params[0][1]}')
    else:
        parser.error('Provide either --task-id, -s1/-q1, or --batch-start/--batch-end')

    # ── Initializing file directory for saving simulations ──
    file_directory = './Unity/Simulations/Part_7_Binary_Collection/'
    os.makedirs(file_directory, exist_ok=True)

    STALE_LOCK_HOURS = 27

    # ── Run simulations ──
    batch_start_time = t.time()
    completed = 0
    skipped = 0

    for sim_idx, (s1_str, q1) in enumerate(sim_params):
        print()
        print('=========================================================')
        print(f'Simulation {sim_idx + 1}/{len(sim_params)} [{mode.upper()}]')
        print('=========================================================')

        s1 = np.float64(numexpr.evaluate(s1_str).item())
        q1 = np.float64(q1)

        print(f'q1 = {q1}')
        print(f's1 = {s1} ({s1_str})')

        binary_lens_attributes = np.array([
            [0, 0, 1],
            [s1, 0, q1]
        ])

        center_of_magnification = np.array([q1 / ((1 + q1) * (s1 + 1/s1)), 0])
        print(f'Binary offset: {center_of_magnification}')

        binary_lens_attributes[:, :2] -= center_of_magnification

        pixels, ang_width, (qs_list_1, ss_list_1), cusp_points, (max_source_radius, min_source_radius) = IRSF.IRSFunctions._ang_width_calculator(binary_lens_attributes, final_multiplier=1, pixels_in_small_source=20, cm_offset='auto')

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

        binary_lens_parameters = {
            'pixels': pixels,
            'ang_width': ang_width,
            'thickness': y_plus - y_minus,
            'y_plus': y_plus,
            'y_minus': y_minus,
            'lens_att': binary_lens_attributes.tolist(),
            'num_theta': num_theta,
            'num_r': num_r
        }

        print(f'Number of rays: {(num_r * num_theta):.4e}')
        print('=========================================================')

        file_name = f'binary_{q1:.0e}_{s1:.2e}.pkl'
        file_path = file_directory + file_name
        lock_path = file_path + '.running'

        print(f'Simulation file path: {file_path}')

        # Skip if already completed by another job
        if os.path.exists(file_path):
            print('Simulation already exists. Skipping...')
            skipped += 1
            continue

        # Atomically claim this simulation via lock file to prevent
        # GPU and CPU jobs from duplicating the same work.
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, f'{mode} pid={os.getpid()}\n'.encode())
            os.close(fd)
        except FileExistsError:
            try:
                lock_age_hours = (t.time() - os.path.getmtime(lock_path)) / 3600
                if lock_age_hours > STALE_LOCK_HOURS:
                    print(f'Stale lock file ({lock_age_hours:.1f}h old). Reclaiming...')
                    os.remove(lock_path)
                    fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                    os.write(fd, f'{mode} pid={os.getpid()} (reclaimed)\n'.encode())
                    os.close(fd)
                else:
                    print(f'Simulation claimed by another job ({lock_age_hours:.1f}h ago). Skipping...')
                    skipped += 1
                    continue
            except (OSError, FileExistsError):
                print('Simulation claimed by another job. Skipping...')
                skipped += 1
                continue

        try:
            print(f'Shooting binary lens [{mode.upper()}]:\n')

            calculator = IRSC.IRSCaustics(annulus_param_dict=binary_lens_parameters)
            magnifications = calculator.series_calculate(
                cm_offset='auto', annulus_offset=center_of_magnification
            )

            print('=========================================================')

            init_time = t.time()
            with open(file_path, 'wb') as calculator_file:
                pickle.dump(calculator, calculator_file)

            print(f'Saving class data to file: {(t.time() - init_time):.3} seconds')
            print(f'Simulation {sim_idx + 1}/{len(sim_params)} done')
            completed += 1

        finally:
            try:
                os.remove(lock_path)
            except OSError:
                pass

    total_time = t.time() - batch_start_time
    print('\n=========================================================')
    print(f'Finished: {completed} computed, {skipped} skipped, {len(sim_params)} total')
    print(f'Total wall time: {total_time:.1f}s ({total_time/3600:.2f}h)')
    if completed > 0:
        print(f'Average per computed simulation: {total_time / completed:.1f}s')
    print('=========================================================')
