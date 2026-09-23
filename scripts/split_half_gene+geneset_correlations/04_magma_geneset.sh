#!/bin/bash
#SBATCH --job-name=magma_geneset
#SBATCH --mem=20G          # memory per node
#SBATCH --time=10:00:00    # max runtime hh:mm:ss
#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)
#SBATCH --output=magma_geneset.%j.out
#SBATCH --error=magma_geneset.%j.err
#module load R/4.4.2-openblas-rocky8

for gmt_file in gene_sets/*.gmt; do
    GMT_BASE="${gmt_file##*/}"
    GMT_BASE="${GMT_BASE%.gmt}"
    OUT_DIR="geneset_output/${GMT_BASE}"
    mkdir -p "$OUT_DIR"
    echo "Gene set file: $gmt_file -> output dir: $OUT_DIR"

    for file in gene_enrichments/*.genes.raw; do
        BASE="${file##*/}"
        BASE="${BASE%.genes.raw}"
        echo "  BASE = $BASE"

        /users/0/buona008/andrew/magma \
            --gene-results /users/0/buona008/andrew/personality_iq_project/split_half/gene_enrichments/"${BASE}".genes.raw \
            --set-annot /users/0/buona008/andrew/personality_iq_project/split_half/gene_sets/"${GMT_BASE}".gmt \
            --out "${OUT_DIR}/${BASE}"
    done
done