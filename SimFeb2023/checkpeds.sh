#!/bin/bash
#SBATCH --account=def-jgraham
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=00:10:00
module load r
Rscript checkpeds.R
