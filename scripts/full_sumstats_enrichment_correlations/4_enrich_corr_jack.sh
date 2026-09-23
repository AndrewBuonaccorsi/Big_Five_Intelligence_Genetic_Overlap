#!/bin/bash
#SBATCH --job-name=enrich_corr_jack
#SBATCH --mem=50G
#SBATCH --time=50:00:00
#SBATCH --output=enrich_corr_jack.out
#SBATCH --error=enrich_corr_jack.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=buona008@umn.edu

module load R/4.4.2-openblas-rocky8

BASE_DIR="/users/0/buona008/andrew/personality_iq_project/jack_enrich_corr"
TEMPDIR="${BASE_DIR}/temp_dir"
BLOCKS=$(seq 1 200)

mkdir -p "$TEMPDIR"
cd "$BASE_DIR" || exit 1

for GROUP_DIR in output/*/; do
    GROUP=$(basename "$GROUP_DIR")
    echo "==> Gene set group: ${GROUP}"
    export GROUP

    for BLOCK in $BLOCKS; do
        cp "${GROUP_DIR}"*/*_jk"${BLOCK}".gsa.out   "$TEMPDIR"/ 2>/dev/null
        cp "${GROUP_DIR}"*/*_jk0"${BLOCK}".gsa.out  "$TEMPDIR"/ 2>/dev/null
        cp "${GROUP_DIR}"*/*_jk00"${BLOCK}".gsa.out "$TEMPDIR"/ 2>/dev/null
        
	cd "$TEMPDIR" || exit 1
        export BLOCK
        Rscript ../unweighted_enrich_corr.R
        rm -rf "$TEMPDIR"/*
        cd "$BASE_DIR" || exit 1

    done

done
