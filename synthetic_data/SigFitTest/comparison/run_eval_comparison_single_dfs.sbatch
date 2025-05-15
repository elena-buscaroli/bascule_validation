#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --array=1-1000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=20G
#SBATCH --job-name=single_dfs_sigfitest
#SBATCH --partition=THIN
#SBATCH --output=logs/single_dfs_sigfitest_%A-%a.log

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'

module load R/4.4.1
conda activate basilica-env

Rscript /u/cdslab/ebusca00/scratch/GitHub/bascule_validation/synthetic_data/SigFitTest/comparison/single_comparison_dataframes.R \
  SigFitTest $SLURM_ARRAY_TASK_ID 0
