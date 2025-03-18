#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=comparison_genmodel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=20G
#SBATCH --time=20:00:00
#SBATCH --output=logs/comparison_genmodel_%j.out

module load R/4.4.1
conda activate bascule-env

Rscript /u/cdslab/ebusca00/scratch/GitHub/bascule_validation/synthetic_data/eval_comparison.R generative_model/all_fits/
