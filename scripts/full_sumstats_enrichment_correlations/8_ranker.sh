#!/usr/bin/env bash
#SBATCH --job-name=geneset_rank_jk
#SBATCH --mem=20G          # memory per node
#SBATCH --time=01:00:00    # max runtime hh:mm:ss
#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)
#SBATCH --output=geneset_rank_jk.%j.out
#SBATCH --error=geneset_rank_jk.%j.err

set -euo pipefail

SCRIPT_DIR=/users/0/buona008/andrew/personality_iq_project/jack_enrich_corr
SELF="${SCRIPT_DIR}/$(basename -- "${BASH_SOURCE[0]}")"
R_SCRIPT="${R_SCRIPT:-${SCRIPT_DIR}/ranker.R}"

BASE="${BASE:-$SCRIPT_DIR}"
NJK="${NJK:-200}"
OUT="${OUT:-results/gene_set_ranking.txt}"
RSCRIPT="${RSCRIPT:-Rscript}"
R_MODULE="${R_MODULE:-}"
TIME="${TIME:-08:00:00}"
ARRAY_TIME="${ARRAY_TIME:-02:00:00}"
COMBINE_TIME="${COMBINE_TIME:-01:00:00}"
MEM="${MEM:-32G}"
SBATCH_EXTRA="${SBATCH_EXTRA:-}"

RESDIR="$BASE/$(dirname "$OUT")"
LISTFILE="$RESDIR/collections.txt"

[[ -f "$R_SCRIPT" ]] || {
  echo "ERROR: cannot find the R script at $R_SCRIPT" >&2
  echo "       Set R_SCRIPT=/path/to/ranker.R if it lives elsewhere." >&2
  exit 1; }

missing=0
for d in "$BASE/output" "$BASE/geneset_output"; do
  [[ -d "$d" ]] || { echo "ERROR: expected directory not found: $d" >&2; missing=1; }
done
if [[ $missing -eq 1 ]]; then
  echo "       Run from jack_enrich_corr, or set BASE=/path/to/jack_enrich_corr" >&2
  exit 1
fi

mkdir -p "$BASE/logs" "$RESDIR"

build_list() {
  comm -12 \
    <(find "$BASE/output"         -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort) \
    <(find "$BASE/geneset_output" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort) \
    > "$LISTFILE"
  wc -l < "$LISTFILE" | tr -d ' '
}

MODE="run"
if [[ "${1:-}" == "--submit" || "${1:-}" == "--submit-array" ]]; then
  MODE="${1#--submit}"; MODE="${MODE:-single}"; shift
fi

if [[ "$MODE" != "run" ]]; then
  command -v sbatch >/dev/null 2>&1 || { echo "ERROR: sbatch not found." >&2; exit 127; }

  if [[ "$MODE" == "single" ]]; then
    # shellcheck disable=SC2086
    jid=$(sbatch --parsable \
            --job-name=ranker --time="$TIME" --mem="$MEM" --cpus-per-task=1 \
            --output="$BASE/logs/ranker_%j.out" \
            $SBATCH_EXTRA "$SELF" "$@")
    echo "submitted job $jid  (parse + combine, --time=$TIME, --mem=$MEM)"
    echo "log: $BASE/logs/ranker_${jid}.out"
  else
    n=$(build_list)
    [[ "$n" -gt 0 ]] || { echo "ERROR: no collections found in both roots." >&2; exit 1; }

    jid=$(sbatch --parsable \
            --job-name=ranker_parse --array="1-${n}" --time="$ARRAY_TIME" \
            --mem="$MEM" --cpus-per-task=1 \
            --output="$BASE/logs/ranker_parse_%A_%a.out" \
            $SBATCH_EXTRA "$SELF" "$@")

    cid=$(sbatch --parsable \
            --job-name=ranker_combine --dependency="afterok:${jid}" \
            --time="$COMBINE_TIME" --mem="$MEM" --cpus-per-task=1 \
            --output="$BASE/logs/ranker_combine_%j.out" \
            $SBATCH_EXTRA "$SELF" --stage combine "$@")
    echo "submitted parse array $jid  ($n collections, --time=$ARRAY_TIME)"
    echo "submitted combine     $cid  (pools all collections after the array succeeds)"
    echo "collection list: $LISTFILE"
    echo "logs: $BASE/logs/ranker_parse_${jid}_*.out"
  fi
  exit 0
fi

EXTRA_ARGS=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  [[ -f "$LISTFILE" ]] || { echo "ERROR: missing $LISTFILE" >&2; exit 1; }
  COLL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LISTFILE")
  [[ -n "$COLL" ]] || { echo "ERROR: no collection at line $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
  EXTRA_ARGS=(--stage parse --collections "$COLL")
fi

if [[ -n "$R_MODULE" ]] && command -v module >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  module load "$R_MODULE"
fi
command -v "$RSCRIPT" >/dev/null 2>&1 || {
  echo "ERROR: '$RSCRIPT' not on PATH. Set RSCRIPT=/path/to/Rscript or R_MODULE=<module>." >&2
  exit 127; }

STAMP="$(date +%Y%m%d_%H%M%S)"
TAG="${SLURM_ARRAY_TASK_ID:-${SLURM_JOB_ID:-local}}"
LOG="$BASE/logs/ranker_${STAMP}_${TAG}.log"

{
  echo "==================================================================="
  echo " gene set ranking: pooled mean Z across IQ and the Big Five"
  echo "-------------------------------------------------------------------"
  echo " started     : $(date)"
  echo " host        : $(hostname)"
  echo " slurm job   : ${SLURM_JOB_ID:-none}${SLURM_ARRAY_TASK_ID:+ (array task $SLURM_ARRAY_TASK_ID)}"
  echo " base        : $BASE"
  echo " R script    : $R_SCRIPT"
  echo " replicates  : $NJK"
  echo " output      : $BASE/$OUT"
  echo " R           : $("$RSCRIPT" --version 2>&1 | head -1)"
  echo " args        : ${EXTRA_ARGS[*]:-} $*"
  echo "==================================================================="
} | tee "$LOG"

set +e
"$RSCRIPT" "$R_SCRIPT" \
    --base "$BASE" --n-jk "$NJK" --out "$OUT" \
    ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} "$@" 2>&1 | tee -a "$LOG"
status=${PIPESTATUS[0]}
set -e

echo "-------------------------------------------------------------------" | tee -a "$LOG"
if [[ $status -eq 0 ]]; then
  echo " finished    : $(date)" | tee -a "$LOG"
  if [[ -f "$BASE/$OUT" ]]; then
    echo " results     : $BASE/$OUT" | tee -a "$LOG"
    echo " top rows:" | tee -a "$LOG"
    head -11 "$BASE/$OUT" | sed 's/^/   /' | tee -a "$LOG"
  fi
else
  echo " FAILED with exit status $status -- see $LOG" | tee -a "$LOG"
fi
echo " log         : $LOG"
exit $status