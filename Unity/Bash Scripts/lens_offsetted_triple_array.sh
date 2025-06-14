#!/usr/bin/env bash

#SBATCH --time=60:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=triple_offsetted_array_s_vary
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-39

#SBATCH --output="./Unity/Output Logs/Collection_pmr0.001/triple_1e11_tripoffset_%a_output.txt"
#SBATCH --error="./Unity/Error Logs/Collection_pmr0.001/triple_1e11_tripoffset_%a_error.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
seperations=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.9)
alphas=(0 45 90 135 180)
pmrs=(0.1 0.01 0.001)

# Compute indices
# alpha_index=$(( SLURM_ARRAY_TASK_ID / 3 ))
# pmrs_index=$(( SLURM_ARRAY_TASK_ID % 3 ))
seperations_index=$(( SLURM_ARRAY_TASK_ID / 5 ))
alpha_index=$(( SLURM_ARRAY_TASK_ID % 5 ))
pmrs_index=2

# Extract parameters
seperation=${seperations[$seperations_index]}
alpha=${alphas[$alpha_index]}
pmr=${pmrs[$pmrs_index]}

echo "Running $alpha degrees with $pmr mass ratio and $seperation seperation"

python "./Unity/Python Scripts/lenses_offsetted.py" -s2 $seperation -a2 $alpha -pmr $pmr -l triple -o triple_offset
