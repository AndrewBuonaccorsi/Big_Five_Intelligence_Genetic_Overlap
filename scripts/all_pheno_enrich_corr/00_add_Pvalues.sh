#!/bin/bash

#SBATCH --job-name=magma_annotate

#SBATCH --mem=20G          # memory per node

#SBATCH --time=10:00:00    # max runtime hh:mm:ss

#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)

#SBATCH --output=magma_anno.%j.out

#SBATCH --error=magma_anno.%j.err

module load R/4.4.2-openblas-rocky8

Rscript add_P.R
