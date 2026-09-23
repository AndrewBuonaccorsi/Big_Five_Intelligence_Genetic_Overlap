#!/bin/bash

#SBATCH --job-name=magma_gene
#SBATCH --mem=100G
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=21
#SBATCH --output=magma_gene.%j.out
#SBATCH --error=magma_gene.%j.err

PHENO_DIR="/projects/standard/leej5/shared/buona008/personality_iq_project/split_half/munged"

ls -1 "$PHENO_DIR"/* || { echo "No .sumstats files found! Exiting."; exit 1; }

for file in "$PHENO_DIR"/*.sumstats; do

    BASE="${file##*/}"          # remove path
    BASE="${BASE%.sumstats}"    # remove suffix

    echo "BASE = $BASE"

    mkdir -p gene_output/"${BASE}_dir"

    for CHR in {1..22}; do
    
        echo "  Running MAGMA for chr $CHR..."
        /projects/standard/leej5/shared/buona008/magma --bfile /projects/standard/leej5/shared/buona008/personality_iq_project/split_half/1000gEUR/g1000_eur synonyms=1000gEUR/g1000_eur.synonyms synonym-dup=drop-dup --gene-annot /projects/standard/leej5/shared/buona008/personality_iq_project/split_half/snps.genes.annot --pval "$file" ncol=N pval=P snp-id=SNP duplicate=first --batch "$CHR" chr --out gene_output/"${BASE}_dir/${BASE}.chr${CHR}"
    
    done


done

wait
