#!/bin/bash
#SBATCH --job-name=personality_corr

#SBATCH --mem=20G          # memory per node

#SBATCH --time=10:00:00    # max runtime hh:mm:ss

#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)

#SBATCH --output=job.%j.out

#SBATCH --error=job.%j.err


module load freetype libpng libtiff libjpeg libwebp pkg-config
module load R/4.4.2-openblas-rocky8

Rscript personality-genetic-corr.R


