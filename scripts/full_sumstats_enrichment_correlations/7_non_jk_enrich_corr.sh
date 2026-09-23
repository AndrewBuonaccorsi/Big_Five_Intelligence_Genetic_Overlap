#!/bin/bash
#SBATCH --job-name=non_jk_corr
#SBATCH --mem=20G
#SBATCH --time=05:00:00
#SBATCH --output=non_jk_corr.out
#SBATCH --error=non_jk_corr.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=buona008@umn.edu

module load R/4.4.2-openblas-rocky8

BASE_DIR="/users/0/buona008/andrew/personality_iq_project/jack_enrich_corr"
cd "$BASE_DIR" || exit 1

mkdir -p non_jk_correlations/unweighted

for GROUP_PATH in geneset_output/*/; do
    GROUP=$(basename "$GROUP_PATH")
    echo "==> Gene set group: ${GROUP}"

    export GROUP
    export GROUP_DIR="${BASE_DIR}/${GROUP_PATH}"

    export OUT_FILE="${BASE_DIR}/non_jk_correlations/unweighted/${GROUP}_corr_matrix.csv"
    Rscript unweighted_enrich_corr_njk.R

done

echo "Done."
