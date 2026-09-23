#!/bin/bash

#SBATCH --job-name=gene_corr_jk

#SBATCH --mem=20G          # memory per node

#SBATCH --time=10:00:00    # max runtime hh:mm:ss

#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)

#SBATCH --output=gene_corr_jk.%j.out

#SBATCH --error=gene_corr_jk.%j.err

Rscript gene_corr.R corr_output dice_output