#!/usr/bin/env bash

#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=binary_animator_1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/Analysis/animator_output.txt"

# Commands to run
module load mamba
mamba activate .venv

python "./Unity Analysis/Spring 2026/frac_map_animator.py"