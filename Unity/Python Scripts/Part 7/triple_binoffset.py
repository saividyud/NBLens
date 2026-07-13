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

    # ── Parameter grid (1×6×2×5 = 60 simulations) ──
    # s1 is fixed at 1.2; s2, q2, and alpha vary.
    S1 = '1.2'
    SS2 = ['0.2', '0.5', '0.8', '1/0.8', '1/0.5', '1/0.2']
    Q2S = [1e-4, 3e-4]
    ALPHAS = [0.0, 45.0, 90.0, 135.0, 180.0]
    TOTAL_SIMS = len(SS2) * len(Q2S) * len(ALPHAS)  # 60

    def task_id_to_params(task_id):
        """Convert a linear task ID (0-59) to (s1_str, s2_str, q2, alpha).

        Ordering: s2 varies slowest, then q2, then alpha (fastest).
        """
        s2_idx = task_id // (len(Q2S) * len(ALPHAS))          # // 10
        q2_idx = (task_id % (len(Q2S) * len(ALPHAS))) // len(ALPHAS)  # (% 10) // 5
        alpha_idx = task_id % len(ALPHAS)                      # % 5
        return S1, SS2[s2_idx], Q2S[q2_idx], ALPHAS[alpha_idx]

    # ── Argument parsing ──
    parser = argparse.ArgumentParser(
        description='Adaptive triple lens magnification map calculator (GPU/CPU).'
    )

    parser.add_argument('-s1', '--sep1', help='Separation of big planet')
    parser.add_argument('-s2', '--sep2', help='Separation of small planet')
    parser.add_argument('-q2', '--massratio2', help='Mass ratio of small planet')
    parser.add_argument('-alpha', '--alpha', help='Angle of small planet')

    parser.add_argument('--batch-start', type=int, default=None,
                        help='First task ID in batch (inclusive, 0-59)')
    parser.add_argument('--batch-end', type=int, default=None,
                        help='Last task ID in batch (exclusive, 1-60)')

    parser.add_argument('--task-id', type=int, default=None,
                        help='Run a single simulation by task ID (0-59)')

    parser.add_argument('--mode', choices=['gpu', 'cpu', 'auto'], default='auto',
                        help='Compute backend (default: auto-detect)')

    parser.add_argument('--reverse', action='store_true',
                        help='Process batch in reverse order (GPU runs 59->0 to avoid overlap with CPU)')

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
            args.sep2,
            np.float64(numexpr.evaluate(args.massratio2).item()),
            np.float64(numexpr.evaluate(args.alpha).item()),
        )]
        print(f'\nSingle mode: s1={args.sep1}, s2={args.sep2}, q2={sim_params[0][2]}, alpha={sim_params[0][3]}')
    else:
        parser.error('Provide either --task-id, -s1/-s2/-q2/-alpha, or --batch-start/--batch-end')

    # ── Initializing file directory for saving simulations ──
    file_directory = './Unity/Simulations/Part_7_Triple_Collection/'
    os.makedirs(file_directory, exist_ok=True)

    STALE_LOCK_HOURS = 27

    # ── Run simulations ──
    batch_start_time = t.time()
    completed = 0
    skipped = 0

    for sim_idx, (s1_str, s2_str, q2, alpha) in enumerate(sim_params):
        print()
        print('=========================================================')
        print(f'Simulation {sim_idx + 1}/{len(sim_params)} [{mode.upper()}]')
        print('=========================================================')

        s1 = np.float64(numexpr.evaluate(s1_str).item())
        s2 = np.float64(numexpr.evaluate(s2_str).item())
        if not isinstance(alpha, float):
            alpha = np.float64(numexpr.evaluate(str(alpha)).item())
        q1 = 1e-3
        q2 = np.float64(q2)

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

        file_name = f'triple_{s2:.2e}_{q2:.2e}_{alpha:.0f}.pkl'
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
            print(f'Shooting triple lens [{mode.upper()}]:\n')

            calculator = IRSC.IRSCaustics(annulus_param_dict=triple_lens_parameters)
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
