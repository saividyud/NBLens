#!/usr/bin/env bash

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --job-name=timing_testing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu
#SBATCH --output="./Unity/Analysis 6-4/bash_logging.txt"

# Commands to run
module load mamba
mamba activate .venv
python "./Unity/Analysis 6-4/timing_test.py"

git add -A
git commit -m "Timing test job"
git push