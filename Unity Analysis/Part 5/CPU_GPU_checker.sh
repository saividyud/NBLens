#!/usr/bin/env bash
# Single triple lens simulation on GPU (no array)
# Submit from project root: sbatch "Unity Analysis/Part 5/CPU_GPU_checker.sh"

#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=triple_cpu_gpu_checker
#SBATCH --mail-type=ALL
#SBATCH --mail-user=senthilnathan.11@osu.edu

#SBATCH --gres=gpu:1

#SBATCH --output="./Unity/Output Logs/Part_5/triple_cpu_gpu_checker_output_%j.txt"

# Commands to run
module load mamba
module load cuda/12.9.0-vvrwkrx
export LD_LIBRARY_PATH=$CUDA_HOME/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
mamba activate .venv

# Clear stale CuPy kernel cache (guards against corrupt .cubin from prior crashes)
rm -rf ~/.cupy/kernel_cache/

echo "Running CPU and GPU checker"

python "./Unity Analysis/Part 5/CPU_GPU_checker.py"
