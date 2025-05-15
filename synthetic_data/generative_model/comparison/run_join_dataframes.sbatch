#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=join_stats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --output=logs/join_stats_%j.out

module load R/4.4.1
conda activate basilica-env

Rscript /u/cdslab/ebusca00/scratch/GitHub/bascule_validation/synthetic_data/generative_model/comparison/join_stats_dataframes.R 
