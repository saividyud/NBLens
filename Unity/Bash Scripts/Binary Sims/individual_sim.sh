#!/usr/bin/env bash

#SBATCH --time=44:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=binary_lens_array
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-2

#SBATCH --output="./Unity/Output Logs/Binary_Collection/binary_%a_output.txt"

# Commands to run
module load mamba
mamba activate .venv

# Defining parameters
if [ $SLURM_ARRAY_TASK_ID -eq 0 ]
then
    s="1/0.4"
    q=1e-03
elif [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
    s="1/0.2"
    q=1e-03
else
    s="1/0.2"
    q=3e-04
fi

echo "Running $s separation with $q mass ratio with origin of binary offset"

python "./Unity/Python Scripts/Binary Sims/binary_binoffset.py" -s $s -q $q
