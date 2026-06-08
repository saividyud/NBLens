#!/usr/bin/env bash

#SBATCH --time=02:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=OGLE_parallel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --array=0-7
#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_6/OGLE_parallel_output_%a.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

# Run simulation
python "./Unity/Python Scripts/Part 6/OGLEsystem_sim.py" --job_id $SLURM_ARRAY_TASK_ID