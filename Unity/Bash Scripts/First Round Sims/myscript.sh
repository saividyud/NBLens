#!/usr/bin/env bash

#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --job-name=timing_testing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --output="./Unity/Output Logs/timing_output.txt"
#SBATCH --error="./Unity/Error Logs/timing_error.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Python Scripts/timing_test.py"

git add -A
git commit -m "Timing test job"
git push
