#!/bin/bash
#SBATCH --job-name=sldsc_h2
#SBATCH --mem=100G
#SBATCH --time=60:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=sldsc.%A_%a.out
#SBATCH --error=sldsc.%A_%a.err
#SBATCH --array=1-98

GENESET_DIR="gene_sets"
ANNOT_DIR="annotations"
LDSCORE_DIR="annotations"

BASELINE_DIR="/projects/standard/leej5/edwa0506/religiosity_sldsc/reference/baselineLD_v2.2"
WEIGHTS_DIR="/projects/standard/leej5/edwa0506/religiosity_sldsc/reference/1000G_Phase3_weights_hm3_no_MHC"
FRQ_DIR="/projects/standard/leej5/edwa0506/religiosity_sldsc/reference/1000G_Phase3_frq"

TRAIT_DIR="sumstats"
OUTDIR="sldsc_results"

mkdir -p ${OUTDIR}


GENESET_PATH=$(ls ${GENESET_DIR}/*.txt | sed -n "${SLURM_ARRAY_TASK_ID}p")
GENESET=$(basename ${GENESET_PATH} .txt)

echo "Gene set: ${GENESET}"


for TRAIT_FILE in ${TRAIT_DIR}/*.sumstats; do

    TRAIT=$(basename ${TRAIT_FILE} .sumstats)

    echo "Trait: ${TRAIT}"

    python /users/0/buona008/andrew/personality_iq_project/sldsc/ldsc/ldsc.py --h2 ${TRAIT_FILE} --ref-ld-chr ${BASELINE_DIR}/baselineLD.,${LDSCORE_DIR}/${GENESET}/${GENESET}. --w-ld-chr ${WEIGHTS_DIR}/weights.hm3_noMHC. --overlap-annot --frqfile-chr ${FRQ_DIR}/1000G.EUR.QC. --out ${OUTDIR}/${GENESET}.${TRAIT}

done
