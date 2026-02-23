#!/usr/bin/env bash

#SBATCH --time=44:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=single_lens_array_4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/Part_4/single_output.txt"

# Commands to run
module load mamba
mamba activate .venv

echo "Running 0.9 separation with 1e-3 mass ratio with origin of (0, 0)"

python "./Unity/Python Scripts/Part 4/single_binoffset.py" -s 0.9 -q 1e-3
