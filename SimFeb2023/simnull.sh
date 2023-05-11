#!/bin/bash
#SBATCH --account=def-jgraham
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=10:00:00
module load r
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
Rscript simnull.R
