#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=comparison_sigfitest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=50G
#SBATCH --time=30:00:00
#SBATCH --output=logs/comparison_sigfitest_%j.out

module load R/4.4.1
conda activate basilica-env

Rscript /u/cdslab/ebusca00/scratch/GitHub/bascule_validation/synthetic_data/eval_comparison.R SigFitTest
