#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=comparison_sigfitest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=20G
#SBATCH --time=20:00:00
#SBATCH --output=logs/comparison_sigfitest_%j.out

module load R/4.4.1
conda activate bascule-env

Rscript /u/cdslab/ebusca00/scratch/GitHub/bascule_validation/synthetic_data/eval_comparison.R SigFitTest
