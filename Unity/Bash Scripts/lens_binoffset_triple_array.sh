#!/usr/bin/env bash

#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=triple_lens_binoffset_array
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-4

#SBATCH --output="./Unity/Output Logs/Collection_0.8/triple_1e11_%a_binoffset_output.txt"
#SBATCH --error="./Unity/Error Logs/Collection_0.8/triple_1e11_%a_binoffset_error.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
alphas=(0 45 90 135 180)
pmrs=(0.1 0.01 0.001)

# Compute indices
alpha_index=$((SLURM_ARRAY_TASK_ID))
pmrs_index=2

# Extract parameters
alpha=${alphas[$alpha_index]}
pmr=${pmrs[$pmrs_index]}

echo "Running $alpha degrees with $pmr mass ratio with origin of binary offset"

python "./Unity/Python Scripts/lenses_offsetted.py" -a2 $alpha -pmr $pmr -l triple -o binary_offset,
