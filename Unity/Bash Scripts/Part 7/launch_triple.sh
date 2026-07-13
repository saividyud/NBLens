#!/usr/bin/env bash
# Launches both GPU and CPU jobs for the adaptive hybrid approach.
#
# How it works:
#   - GPU job: 1 SLURM job that loops through all 60 sims in REVERSE
#     order (59 -> 0).
#   - CPU job: 60-task array with up to 40 running concurrently in
#     FORWARD order (0 -> 59).
#   - Both jobs write to the same output directory and use atomic lock
#     files to avoid duplicating work. Whichever job finishes a sim
#     first wins; the other skips it.
#
# Run from project root: bash "Unity/Bash Scripts/Part 7/launch_triple.sh"

SCRIPT_DIR="Unity/Bash Scripts/Part 7"

echo "=============================================="
echo "  Adaptive Hybrid Triple Lens Launcher"
echo "=============================================="
echo ""

echo "Submitting GPU series job (1 job, sims 59 -> 0)..."
GPU_JOB=$(sbatch --parsable "${SCRIPT_DIR}/triple_binoffset_gpu.sh")
echo "  GPU job ID: ${GPU_JOB}"

echo "Submitting CPU array job (60 tasks, 40 concurrent, sims 0 -> 59)..."
CPU_JOB=$(sbatch --parsable "${SCRIPT_DIR}/triple_binoffset_cpu.sh")
echo "  CPU array job ID: ${CPU_JOB}"

echo ""
echo "Both jobs submitted. Lock files prevent duplicate work."
echo "Monitor with: squeue -u $USER"
echo "=============================================="
