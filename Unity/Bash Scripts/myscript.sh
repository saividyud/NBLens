#!/usr/bin/env bash

#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --job-name=timing_testing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu
#SBATCH --output="./Unity/Analysis 6-4/timing_logging.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Analysis 6-4/timing_test.py"

git add -A
git commit -m "Timing test job"
git push
