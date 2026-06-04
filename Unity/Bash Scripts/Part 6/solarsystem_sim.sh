#!/usr/bin/env bash
# Solar system simulation on GPU (no array)
# Submit from project root: sbatch "Unity/Bash Scripts/Part 6/solarsystem_sim.sh"

#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=solarsystem_sim
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_6/solarsystem_sim.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

# Run simulation
python "./Unity/Python Scripts/Part 6/solarsystem_sim.py" --job_id 6