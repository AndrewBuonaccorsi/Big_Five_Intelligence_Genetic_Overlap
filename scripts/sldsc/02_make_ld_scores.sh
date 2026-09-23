#!/bin/bash
#SBATCH --job-name=sldsc_ld
#SBATCH --mem=100G
#SBATCH --time=60:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=sldsc_ld.%A_%a.out
#SBATCH --error=sldsc_ld.%A_%a.err
#SBATCH --array=1-98

GENESET_DIR="gene_sets"
ANNOT_DIR="annotations"
OUTDIR="annotations"

mkdir -p ${OUTDIR}

# pick gene set based on array index
GENESET_PATH=$(ls ${GENESET_DIR}/*.txt | sed -n "${SLURM_ARRAY_TASK_ID}p")
GENESET=$(basename ${GENESET_PATH} .txt)

mkdir -p ${OUTDIR}/${GENESET}

for CHR in {1..22}; do

  python /users/0/buona008/andrew/personality_iq_project/sldsc/ldsc/ldsc.py \
    --l2 \
    --bfile /users/0/buona008/andrew/personality_iq_project/sldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR} \
    --thin-annot \
    --ld-wind-kb 1000 \
    --annot ${ANNOT_DIR}/${GENESET}/${GENESET}.${CHR}.annot.gz \
    --out ${OUTDIR}/${GENESET}/${GENESET}.${CHR} \
    --print-snps w_hm3.clean.snplist

done
