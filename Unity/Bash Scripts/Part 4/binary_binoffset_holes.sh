#!/usr/bin/env bash

#SBATCH --time=60:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=binary_lens_array_4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-13

#SBATCH --output="./Unity/Output Logs/Part_4/binary_%a_2_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters that did not finish running
# File binary_3e-06_3.00e-01.pkl does NOT exist.
# File binary_1e-05_3.00e-01.pkl does NOT exist.
# File binary_3e-05_3.00e-01.pkl does NOT exist.
# File binary_1e-04_3.00e-01.pkl does NOT exist.
# File binary_3e-04_3.00e-01.pkl does NOT exist.
# File binary_1e-03_3.00e-01.pkl does NOT exist.
# File binary_1e-06_4.00e-01.pkl does NOT exist.
# File binary_3e-06_4.00e-01.pkl does NOT exist.
# File binary_1e-05_4.00e-01.pkl does NOT exist.
# File binary_3e-05_4.00e-01.pkl does NOT exist.
# File binary_1e-04_4.00e-01.pkl does NOT exist.
ss=(0.3 0.4)
qs=(1e-6 3e-6 1e-5 3e-5 1e-4 3e-4 1e-3)

# Compute indices
s_index=$((SLURM_ARRAY_TASK_ID / 7))
q_index=$((SLURM_ARRAY_TASK_ID % 7))

# Extract parameters
s=${ss[$s_index]}
q=${qs[$q_index]}

echo "Running $s separation with $q mass ratio with origin of binary offset"

python "./Unity/Python Scripts/Part 4/binary_binoffset.py" -s $s -q $q
