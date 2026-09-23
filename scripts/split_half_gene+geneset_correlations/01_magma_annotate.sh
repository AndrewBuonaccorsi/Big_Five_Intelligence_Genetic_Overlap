#!/bin/bash

#SBATCH --job-name=magma_annotate

#SBATCH --mem=20G          # memory per node

#SBATCH --time=10:00:00    # max runtime hh:mm:ss

#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)

#SBATCH --output=magma_anno.%j.out

#SBATCH --error=magma_anno.%j.err

GENELOC='NCBI37.3.gene.loc'

/users/0/buona008/andrew/magma --annotate window=10 --snp-loc /users/0/buona008/andrew/personality_iq_project/split_half/1000gEUR/EUR.bim --gene-loc /users/0/buona008/andrew/personality_iq_project/split_half/${GENELOC} --out snps

