#!/bin/bash
#SBATCH --job-name=enrich_corr_heatmaps
#SBATCH --mem=30G
#SBATCH --time=05:00:00
#SBATCH --output=enrich_corr_heatmaps.out
#SBATCH --error=enrich_corr_heatmaps.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=buona008@umn.edu

module load R/4.4.2-openblas-rocky8

BASE_DIR="/users/0/buona008/andrew/personality_iq_project/all_phenotypes_enrich_corr"
GENESET_DIR="${BASE_DIR}/geneset_output"
OUT_DIR="${BASE_DIR}/heatmaps"

mkdir -p "${OUT_DIR}/unweighted"

cd "$BASE_DIR" || exit 1

for GROUP_PATH in "${GENESET_DIR}"/*/; do
    GROUP=$(basename "$GROUP_PATH")
    echo "==> Gene set group: ${GROUP}"

    export GROUP
    export GROUP_DIR="${GROUP_PATH}"

    export OUT_FILE="${OUT_DIR}/unweighted/${GROUP}_unweighted.pdf"
    Rscript unweighted_heatmap.R

done

echo "Done."
