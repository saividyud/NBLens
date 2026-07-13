#!/usr/bin/env bash
# Local test script for triple_binoffset.sh (no SLURM)
# Run from project root: bash "Unity/Bash Scripts/Part 7/test_triple_binoffset_local.sh"
#
# Usage:
#   ./test_triple_binoffset_local.sh           # Dry run: print params for all 60 tasks
#   ./test_triple_binoffset_local.sh 0         # Run task 0 only
#   ./test_triple_binoffset_local.sh 0 1 2    # Run tasks 0, 1, 2
#   ./test_triple_binoffset_local.sh --dry    # Dry run (same as no args)

# Defining parameters (must match triple_binoffset_parallel.sh; s1 fixed at 1.2)
s1=1.2
ss2=(0.2 0.5 0.8 1/0.8 1/0.5 1/0.2)
q2s=(1e-4 3e-4)
alphas=(0 45 90 135 180)

dry_run=false
task_ids=()

for arg in "$@"; do
    if [[ "$arg" == "--dry" ]]; then
        dry_run=true
    elif [[ "$arg" =~ ^[0-9]+$ ]]; then
        task_ids+=("$arg")
    fi
done

# If no task IDs given, dry run all 60
if [[ ${#task_ids[@]} -eq 0 ]]; then
    dry_run=true
    task_ids=($(seq 0 59))
fi

for SLURM_ARRAY_TASK_ID in "${task_ids[@]}"; do
    # Same index logic as triple_binoffset_parallel.sh
    s2_index=$((SLURM_ARRAY_TASK_ID / 10))
    q2_index=$(((SLURM_ARRAY_TASK_ID % 10) / 5))
    alpha_index=$((SLURM_ARRAY_TASK_ID % 5))

    s2=${ss2[$s2_index]}
    q2=${q2s[$q2_index]}
    alpha=${alphas[$alpha_index]}

    if [[ "$dry_run" == true ]]; then
        echo "Task $SLURM_ARRAY_TASK_ID: s1=$s1 s2=$s2 q2=$q2 alpha=$alpha"
    else
        echo "Running task $SLURM_ARRAY_TASK_ID: s1=$s1 s2=$s2 q2=$q2 alpha=$alpha"
        python "./Unity/Python Scripts/Part 7/triple_binoffset.py" -s1 "$s1" -s2 "$s2" -q2 "$q2" -alpha "$alpha"
    fi
done
