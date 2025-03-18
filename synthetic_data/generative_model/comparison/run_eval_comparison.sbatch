#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=comparison
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=comparison_%j.out

module load R/4.4.1
conda activate bascule-env

Rscript /u/cdslab/ebusca00/scratch/GitHub/bascule_validation/synthetic_data/generative_model/comparison/eval_comparison.R
