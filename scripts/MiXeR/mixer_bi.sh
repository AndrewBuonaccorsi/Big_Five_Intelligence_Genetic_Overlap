#!/bin/bash
#SBATCH --job-name=mixer_bi
#SBATCH --output=logs/MIXER_BI-%A_%a.txt
#SBATCH --error=logs/MIXER_BI-%A_%a.txt
#SBATCH --open-mode=truncate
#SBATCH -p msismall
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=800M
#SBATCH --array=1-20
#SBATCH --mail-type=BEGIN,END,FAIL  
#SBATCH --mail-user=edwa0506@umn.edu 

set -euo pipefail

cd /projects/standard/leej5/edwa0506/andrew_project/gsa

module load singularity

export MIXER_SIF=/projects/standard/leej5/edwa0506/andrew_project/gsa/mixer.sif
export MIXER_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

REP="rep${SLURM_ARRAY_TASK_ID}"

mkdir -p out logs

LD="reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld"
BIM="reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
EXTRACT="reference/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.${REP}.snps"

declare -A SUMSTATS=(
  [OPEN]="summary_statistics/open_schwaba2025.sumstats.gz"
  [CONSC]="summary_statistics/consc_schwaba2025.sumstats.gz"
  [EXTRA]="summary_statistics/extra_schwaba2025.sumstats.gz"
  [AGREE]="summary_statistics/agree_schwaba2025.sumstats.gz"
  [NEURO]="summary_statistics/neurot_schwaba2025.sumstats.gz"
  [IQ]="summary_statistics/iq_savage2018.sumstats.gz"
  [HEIGHT]="summary_statistics/height_yengo2022.sumstats.gz"
)

# Choose the pairs you want:
PAIRS=(
  "OPEN IQ"
  "CONSC IQ"
  "EXTRA IQ"
  "AGREE IQ"
  "NEURO IQ"
  "OPEN HEIGHT"
  "CONSC HEIGHT"
  "EXTRA HEIGHT"
  "AGREE HEIGHT"
  "NEURO HEIGHT"
  "IQ HEIGHT"
)

for P in "${PAIRS[@]}"; do
  set -- $P
  T1="$1"
  T2="$2"

  F1="${SUMSTATS[$T1]}"
  F2="${SUMSTATS[$T2]}"

  # Requires univariate JSONs to exist:
  J1="out/${T1}.fit.${REP}.json"
  J2="out/${T2}.fit.${REP}.json"

  ${MIXER_PY} fit2 \
    --ld-file "${LD}" \
    --bim-file "${BIM}" \
    --threads 16 \
    --extract "${EXTRACT}" \
    --trait1-file "${F1}" \
    --trait2-file "${F2}" \
    --trait1-params "${J1}" \
    --trait2-params "${J2}" \
    --out "out/${T1}_vs_${T2}.fit.${REP}"

  ${MIXER_PY} test2 \
    --ld-file "${LD}" \
    --bim-file "${BIM}" \
    --threads 16 \
    --trait1-file "${F1}" \
    --trait2-file "${F2}" \
    --load-params "out/${T1}_vs_${T2}.fit.${REP}.json" \
    --out "out/${T1}_vs_${T2}.test.${REP}"
done
