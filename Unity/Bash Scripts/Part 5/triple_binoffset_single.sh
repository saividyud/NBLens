#!/usr/bin/env bash
# Single triple lens simulation on GPU (no array)
# Submit from project root: sbatch "Unity/Bash Scripts/Part 5/triple_binoffset_single.sh"

#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=triple_single_5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_5/triple_single_output.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

# Single simulation parameters (edit these to run different configurations)
s1=0.2
s2=0.2
alpha=45

echo "Running single triple lens: s1=$s1 s2=$s2 alpha=$alpha"

python "./Unity/Python Scripts/Part 5/triple_binoffset.py" -s1 $s1 -s2 $s2 -alpha $alpha
