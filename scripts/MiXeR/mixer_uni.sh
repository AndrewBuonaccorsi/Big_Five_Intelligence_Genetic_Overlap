#!/bin/bash
#SBATCH --job-name=mixer_uni
#SBATCH --output=logs/MIXER_UNI-%A_%a.txt
#SBATCH --error=logs/MIXER_UNI-%A_%a.txt
#SBATCH --open-mode=truncate
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=800M
#SBATCH --array=1-20

set -euo pipefail

# Workdir
cd /projects/standard/leej5/edwa0506/andrew_project/gsa

module load singularity

export MIXER_SIF=/projects/standard/leej5/edwa0506/andrew_project/gsa/mixer.sif
export MIXER_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

REP="rep${SLURM_ARRAY_TASK_ID}"

mkdir -p out logs

LD="reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld"
BIM="reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
EXTRACT="reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.${REP}.snps"

# ------------------------------------------------------------
# Trait list (edit filenames here; paths are under summary_statistics/)
# ------------------------------------------------------------
declare -A SUMSTATS=(
  [OPEN]="summary_statistics/open_schwaba2025.sumstats.gz"
  [CONSC]="summary_statistics/consc_schwaba2025.sumstats.gz"
  [EXTRA]="summary_statistics/extra_schwaba2025.sumstats.gz"
  [AGREE]="summary_statistics/agree_schwaba2025.sumstats.gz"
  [NEURO]="summary_statistics/neurot_schwaba2025.sumstats.gz"
  [IQ]="summary_statistics/iq_savage2018.sumstats.gz"
  [HEIGHT]="summary_statistics/height_yengo2022.sumstats.gz"
)

# ------------------------------------------------------------
# Univariate fits + tests
# ------------------------------------------------------------
for T in "${!SUMSTATS[@]}"; do
  F="${SUMSTATS[$T]}"

  ${MIXER_PY} fit1 \
    --ld-file "${LD}" \
    --bim-file "${BIM}" \
    --threads 16 \
    --extract "${EXTRACT}" \
    --trait1-file "${F}" \
    --out "out/${T}.fit.${REP}"

  ${MIXER_PY} test1 \
    --ld-file "${LD}" \
    --bim-file "${BIM}" \
    --threads 16 \
    --trait1-file "${F}" \
    --load-params "out/${T}.fit.${REP}.json" \
    --out "out/${T}.test.${REP}"
done
