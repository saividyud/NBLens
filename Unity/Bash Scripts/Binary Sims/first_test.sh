#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=binary_lens
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/Binary_Collection/binary_output.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Python Scripts/Binary Sims/binary_binoffset.py" -s1 0.8 -q1 0.001
