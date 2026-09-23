#!/bin/bash
#SBATCH --job-name=gene_corr_jk
#SBATCH --mem=20G          # memory per node
#SBATCH --time=10:00:00    # max runtime hh:mm:ss
#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)
#SBATCH --output=gene_corr_jk.%j.out
#SBATCH --error=gene_corr_jk.%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=buona008@umn.edu

module purge
module load R/4.4.2-openblas-rocky8

Rscript gene_zstat_correlation.R file_list.txt corr_ dice_
Rscript gene_corr_no_jk.R file_list.txt corr_njk_ dice_njk_
