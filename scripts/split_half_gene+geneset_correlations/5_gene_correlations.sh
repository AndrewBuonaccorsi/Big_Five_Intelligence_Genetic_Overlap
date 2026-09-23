#!/bin/bash
#SBATCH --job-name=gene_corr

#SBATCH --mem=10G

#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=21

#SBATCH --output=gene_corr.%j.out

#SBATCH --error=gene_corr.%j.err

module purge
module load R/4.4.2-openblas-rocky8


set -e




FILE_LIST="genes_files.txt"




find merged_output -type f -path "*.genes.out" > $FILE_LIST

Rscript gene_zstat_correlation.R genes_files.txt correlation_output.txt dice_output.txt

echo "Correlation matrix written to zstat_correlation_matrix.txt"


