#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=gen_obj
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --output=run_gen_obj_SigFitTest_%j.out

module load R/4.4.1
conda activate bascule-env

Rscript /u/cdslab/ebusca00/scratch/GitHub/bascule_validation/synthetic_data/gen_bascule_obj_SigFitTest.R
