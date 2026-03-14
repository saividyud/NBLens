#!/usr/bin/env bash
# Launches both GPU and CPU jobs for the adaptive hybrid approach.
#
# How it works:
#   - GPU job: 1 SLURM job that loops through all 125 sims in REVERSE
#     order (124 -> 0), completing ~1 sim every 1h 45m.
#   - CPU job: 125-task array with up to 40 running concurrently in
#     FORWARD order (0 -> 124), each sim taking ~10h.
#   - Both jobs write to the same output directory and use atomic lock
#     files to avoid duplicating work. Whichever job finishes a sim
#     first wins; the other skips it.
#
# Expected wall time: ~30 hours (vs ~40h CPU-only, ~219h GPU-only).
#
# Run from project root: bash "Unity/Bash Scripts/Part 5/launch_triple.sh"

SCRIPT_DIR="Unity/Bash Scripts/Part 5"

echo "=============================================="
echo "  Adaptive Hybrid Triple Lens Launcher"
echo "=============================================="
echo ""

echo "Submitting GPU series job (1 job, sims 124 -> 0)..."
GPU_JOB=$(sbatch --parsable "${SCRIPT_DIR}/triple_binoffset_gpu.sh")
echo "  GPU job ID: ${GPU_JOB}"

echo "Submitting CPU array job (125 tasks, 40 concurrent, sims 0 -> 124)..."
CPU_JOB=$(sbatch --parsable "${SCRIPT_DIR}/triple_binoffset_cpu.sh")
echo "  CPU array job ID: ${CPU_JOB}"

echo ""
echo "Both jobs submitted. Lock files prevent duplicate work."
echo "Monitor with: squeue -u $USER"
echo "=============================================="
