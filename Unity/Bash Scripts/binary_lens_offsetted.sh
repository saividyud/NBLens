#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --job-name=binary_lens_1e11_offsetted_45deg
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/binary_lens_1e11_offsetted_45deg_output.txt"
#SBATCH --error="./Unity/Error Logs/binary_lens_1e11_offsetted_45deg_error.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Python Scripts/binary_lens_offsetted.py"