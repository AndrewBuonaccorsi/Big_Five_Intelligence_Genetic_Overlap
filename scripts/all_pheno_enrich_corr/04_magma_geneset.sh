#!/bin/bash
#SBATCH --job-name=magma_geneset
#SBATCH --mem=20G          # memory per node
#SBATCH --time=20:00:00    # max runtime hh:mm:ss
#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)
#SBATCH --array=0-22              # number of .gmt files in gene_sets minus 1
#SBATCH --output=magma_geneset_%A_%a.out
#SBATCH --error=magma_geneset_%A_%a.err
#module load R/4.4.2-openblas-rocky8

mkdir -p geneset_output

GMT_FILES=(gene_sets/*.gmt)
gmt_file=${GMT_FILES[$SLURM_ARRAY_TASK_ID]}

GMT_BASE="${gmt_file##*/}"
GMT_BASE="${GMT_BASE%.gmt}"
OUT_DIR="geneset_output/${GMT_BASE}"
mkdir -p "$OUT_DIR"
echo "Gene set file: $gmt_file -> output dir: $OUT_DIR"

for file in merged_output/*.genes.raw; do
    BASE="${file##*/}"
    BASE="${BASE%.genes.raw}"
    echo "  BASE = $BASE"

    /users/0/buona008/andrew/magma --gene-results /users/0/buona008/andrew/personality_iq_project/all_phenotypes_enrich_corr/merged_output/"${BASE}".genes.raw --set-annot "$gmt_file" --out "${OUT_DIR}/${BASE}"
done
