"""
Benchmark script for GPU magnification map computation.

Runs the same simulation as CPU_GPU_checker.py but focuses on
timing the series_calculate method over multiple trials.

Usage (from project root):
    python "Unity Analysis/Part 5/gpu_benchmark.py"
    python "Unity Analysis/Part 5/gpu_benchmark.py" --trials 5
"""

import sys
sys.path.append('.')

import argparse
import time
import numpy as np

import IRSMicroLensing.IRSCausticsGPU as IRSCGPU

try:
    import cupy as cp
except ImportError:
    cp = None


def build_simulation():
    """Set up the same lens configuration as CPU_GPU_checker.py."""
    q1 = 1e-3
    s1 = 0.2
    q2 = 3e-4
    s2 = 0.6
    alpha = 45

    triple_lens_att = np.array([
        [0, 0, 1],
        [s1, 0, q1],
        [s2 * np.cos(np.radians(alpha)), s2 * np.sin(np.radians(alpha)), q2]
    ])

    center_of_magnification = np.array([
        q2 / ((1 + q2) * (s1 + 1 / s1)), 0
    ])
    triple_lens_att[:, :2] -= center_of_magnification

    params = {
        'lens_att': triple_lens_att.tolist(),
        'pixels': 1000,
        'ang_width': 0.2,
        'num_r': 185900,
        'num_theta': 46475,
        'y_plus': 1.0769985822006392,
        'y_minus': 0.9285063760849294,
        'thickness': 1.0769985822006392 - 0.9285063760849294
    }

    return params, center_of_magnification


def run_trial(params, center_of_magnification, trial_num, total_trials):
    """Run one simulation and return wall-clock time in seconds."""
    if cp is not None:
        cp.get_default_memory_pool().free_all_blocks()
        cp.cuda.Device().synchronize()

    sim = IRSCGPU.IRSCaustics(annulus_param_dict=params)

    start = time.perf_counter()
    sim.series_calculate(
        cm_offset='auto',
        annulus_offset=center_of_magnification,
        print_stats=(trial_num == 1)
    )
    if cp is not None:
        cp.cuda.Device().synchronize()
    elapsed = time.perf_counter() - start

    print(f"  Trial {trial_num}/{total_trials}: {elapsed:.3f}s")
    return elapsed


def main():
    parser = argparse.ArgumentParser(description='GPU magnification map benchmark')
    parser.add_argument('--trials', type=int, default=3,
                        help='Number of timed trials (default: 3)')
    parser.add_argument('--warmup', type=int, default=1,
                        help='Number of warmup runs (default: 1, not timed)')
    args = parser.parse_args()

    params, cm_offset = build_simulation()

    total_rays = params['num_r'] * params['num_theta']
    print("=" * 60)
    print("GPU Magnification Map Benchmark")
    print("=" * 60)
    print(f"  Pixels:      {params['pixels']}x{params['pixels']}")
    print(f"  num_r:       {params['num_r']}")
    print(f"  num_theta:   {params['num_theta']}")
    print(f"  Total rays:  {total_rays:.3e}")
    print(f"  Warmup runs: {args.warmup}")
    print(f"  Timed trials: {args.trials}")
    print("=" * 60)

    # Warmup
    for i in range(args.warmup):
        print(f"\n--- Warmup {i + 1}/{args.warmup} ---")
        run_trial(params, cm_offset, i + 1, args.warmup)

    # Timed trials
    times = []
    for i in range(args.trials):
        print(f"\n--- Trial {i + 1}/{args.trials} ---")
        t = run_trial(params, cm_offset, i + 1, args.trials)
        times.append(t)

    # Results
    times = np.array(times)
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"  Trials:  {args.trials}")
    print(f"  Times:   {', '.join(f'{t:.3f}s' for t in times)}")
    print(f"  Min:     {times.min():.3f}s")
    print(f"  Max:     {times.max():.3f}s")
    print(f"  Mean:    {times.mean():.3f}s")
    print(f"  Std:     {times.std():.3f}s")
    print(f"  Median:  {np.median(times):.3f}s")
    print(f"  Rays/sec (median): {total_rays / np.median(times):.3e}")
    print("=" * 60)


if __name__ == '__main__':
    main()
