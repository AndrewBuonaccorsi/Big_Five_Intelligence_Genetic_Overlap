#!/bin/bash
#SBATCH --job-name=magma_jack_traits
#SBATCH --mem=50G
#SBATCH --time=40:00:00
#SBATCH --array=0-5              # This number of traits minus
#SBATCH --output=magma_jack_%A_%a.out
#SBATCH --error=magma_jack_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=buona008@umn.edu

cd /users/0/buona008/andrew/personality_iq_project/jack_enrich_corr

TRAITS=(agree consc extra neurot open iq)
N_JK=200

trait=${TRAITS[$SLURM_ARRAY_TASK_ID]}
genesraw="gene_enrichments/${trait}.genes.raw"
EXCLUDE_DIR="excluded_genes"

echo "==> Starting trait: ${trait}"
echo "    Using gene results: ${genesraw}"

for gmt_file in gene_sets/*.gmt; do
    GMT_BASE="${gmt_file##*/}"
    GMT_BASE="${GMT_BASE%.gmt}"
    OUTDIR="output/${GMT_BASE}/jack_${trait}"
    mkdir -p "$OUTDIR"
    echo "  ==> Gene set file: ${gmt_file}"

    for k in $(seq -w 1 ${N_JK}); do
        excl="${EXCLUDE_DIR}/exclude_genes_${k}.txt"
        if [ ! -f "$excl" ]; then
            echo "    !! Missing exclusion file: $excl, skipping"
            continue
        fi
        out_prefix="${OUTDIR}/${trait}_jk${k}"
        echo "    - JK ${k}: exclude file ${excl}"
        /users/0/buona008/andrew/magma --gene-results "$genesraw" --set-annot "$gmt_file" --settings gene-exclude="$excl" --out "${out_prefix}"
    done
done

echo "==> Finished trait: ${trait}"
