#!/bin/bash

#SBATCH --job-name=sldsc_annotate

#SBATCH --mem=20G          # memory per node

#SBATCH --time=20:00:00    # max runtime hh:mm:ss

#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)

#SBATCH --output=sldsc_annot.%j.out

#SBATCH --error=sldsc_annot.%j.err

GENESET_DIR="gene_sets"
OUTDIR="annotations"

for GENESET in ${GENESET_DIR}/*.txt; do

  BASE=$(basename "${GENESET}" .txt)

  mkdir -p "${OUTDIR}/${BASE}"

  for CHR in {1..22}; do
    python /users/0/buona008/andrew/personality_iq_project/sldsc/ldsc/make_annot.py \
      --gene-set-file "${GENESET}" \
      --gene-coord-file gene_coords.txt \
      --windowsize 10000 \
      --bimfile /users/0/buona008/andrew/personality_iq_project/sldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}.bim \
      --annot-file "${OUTDIR}/${BASE}/${BASE}.${CHR}.annot.gz"
  done

done
