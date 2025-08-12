#!/bin/bash

#SBATCH --job-name=test_genome_utils
#SBATCH --account=root
#SBATCH --cpus-per-task=64
#SBATCH --output=slurm_out/%x/slurm-%j.out


cd /home/ayat/Repositories/genome_utils
source /home/ayat/.venv/ASODesignPipeline/bin/activate

echo "Starting SLURM job"
date
echo "Running on host: $(hostname)"
echo "Executing command: python test_genome_utils.py"

python test_genome_utils.py

echo "SLURM job finished"
date 