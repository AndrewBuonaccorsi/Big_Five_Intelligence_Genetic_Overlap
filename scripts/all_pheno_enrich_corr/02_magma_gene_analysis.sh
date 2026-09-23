#!/bin/bash

#SBATCH --job-name=magma_gene
#SBATCH --mem=100G
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=21
#SBATCH --output=magma_gene.%j.out
#SBATCH --error=magma_gene.%j.err

PHENO_DIR="/users/0/buona008/andrew/all_phenotypes_enrich_corr/phenotypes"

for file in "$PHENO_DIR"/*.sumstats; do

    echo "Processing file: $file"

    BASE="${file##*/}"
    BASE="${BASE%.sumstats}"

    echo "BASE = $BASE"

    mkdir -p gene_output/"${BASE}_dir"

    for CHR in {1..22}; do
    
        echo "  Running MAGMA for chr $CHR..."
        /users/0/buona008/andrew/magma --bfile /users/0/buona008/andrew/enrichment_correlations/g1000_eur/g1000_eur --gene-annot snps.genes.annot --pval "$file" ncol=N --batch "$CHR" chr --out gene_output/"${BASE}_dir/${BASE}.chr${CHR}"
    
    done

    wait

done

wait
