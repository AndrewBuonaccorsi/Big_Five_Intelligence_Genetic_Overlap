#!/usr/bin/env bash
#SBATCH --job-name=enrich_diff
#SBATCH --mem=20G          # memory per node
#SBATCH --time=01:00:00    # max runtime hh:mm:ss
#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)
#SBATCH --output=enrich_diff.%j.out
#SBATCH --error=enrich_diff.%j.err


set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "$SLURM_SUBMIT_DIR"; fi

SCRIPT_DIR=/users/0/buona008/andrew/personality_iq_project/jack_enrich_corr
R_SCRIPT="${R_SCRIPT:-${SCRIPT_DIR}/enrich_diff.R}"
BASE="${BASE:-$SCRIPT_DIR}"
RSCRIPT="${RSCRIPT:-Rscript}"
R_MODULE="${R_MODULE:-}"

CACHE="$BASE/results/cache"

[[ -f "$R_SCRIPT" ]] || {
  echo "ERROR: cannot find the R script at $R_SCRIPT" >&2
  echo "       Set R_SCRIPT=/path/to/enrich_diff.R if it lives elsewhere." >&2
  exit 1; }

if [[ ! -d "$CACHE" ]] || ! compgen -G "$CACHE/*.beta.rds" >/dev/null; then
  echo "ERROR: no gene set cache in $CACHE" >&2
  echo "       enrich_diff.R reuses the cache ranker.R builds, so that the z" >&2
  echo "       scores match gene_set_ranking.txt exactly. Build it first with:" >&2
  echo "           ./run_ranker.sh --submit-array          # or" >&2
  echo "           ./run_ranker.sh --stage parse --resume" >&2
  exit 1
fi

if ! "${RSCRIPT:-Rscript}" -e "invisible(parse('$R_SCRIPT'))" >/dev/null 2>&1; then
  if command -v "${RSCRIPT:-Rscript}" >/dev/null 2>&1; then
    echo "ERROR: $R_SCRIPT does not parse. The file is probably corrupted." >&2
    echo "       Re-download it and copy it over as a binary file (scp / sftp)," >&2
    echo "       not by pasting into a terminal editor. Check with:" >&2
    echo "           grep -nP '[^\\x09\\x0A\\x20-\\x7E]' $R_SCRIPT   # must print nothing" >&2
    echo "           md5sum $R_SCRIPT" >&2
    "${RSCRIPT:-Rscript}" -e "invisible(parse('$R_SCRIPT'))" 2>&1 | head -5 >&2
    exit 1
  fi
fi

if [[ -n "$R_MODULE" ]] && command -v module >/dev/null 2>&1; then
  module load "$R_MODULE"
fi
command -v "$RSCRIPT" >/dev/null 2>&1 || {
  echo "ERROR: '$RSCRIPT' not on PATH. Set RSCRIPT=/path/to/Rscript or R_MODULE=<module>." >&2
  exit 127; }

mkdir -p "$BASE/logs"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="$BASE/logs/enrich_diff_${STAMP}_${SLURM_JOB_ID:-local}.log"

{
  echo "==================================================================="
  echo " focal trait vs. mean of the other five: jackknife SE identity test"
  echo "-------------------------------------------------------------------"
  echo " started    : $(date)"
  echo " host       : $(hostname)"
  echo " slurm job  : ${SLURM_JOB_ID:-none}"
  echo " base       : $BASE"
  echo " R script   : $R_SCRIPT"
  echo " cache      : $CACHE  ($(ls "$CACHE"/*.beta.rds 2>/dev/null | wc -l | tr -d ' ') collections)"
  echo " R          : $("$RSCRIPT" --version 2>&1 | head -1)"
  echo " args       : $*"
  echo "==================================================================="
} | tee "$LOG"

set +e
"$RSCRIPT" "$R_SCRIPT" --base "$BASE" "$@" 2>&1 | tee -a "$LOG"
status=${PIPESTATUS[0]}
set -e

echo "-------------------------------------------------------------------" | tee -a "$LOG"
if [[ $status -eq 0 ]]; then
  OUT="$BASE/results/enrich_diff"
  echo " finished   : $(date)" | tee -a "$LOG"
  echo " results    : $OUT" | tee -a "$LOG"
  echo "   figures/            $(ls "$OUT/figures" 2>/dev/null | wc -l | tr -d ' ') pdf" | tee -a "$LOG"
  echo "   all_results/        $(ls "$OUT/all_results" 2>/dev/null | wc -l | tr -d ' ') csv" | tee -a "$LOG"
  echo "   significant_results/ $(ls "$OUT/significant_results" 2>/dev/null | wc -l | tr -d ' ') csv" | tee -a "$LOG"
  if [[ -f "$OUT/enrich_diff_diagnostics.txt" ]]; then
    echo " diagnostics:" | tee -a "$LOG"
    sed -n '/SIGNIFICANT DIFFERENCES BY FOCAL TRAIT/,/^$/p' \
        "$OUT/enrich_diff_diagnostics.txt" | sed 's/^/   /' | tee -a "$LOG"
  fi
else
  echo " FAILED with exit status $status -- see $LOG" | tee -a "$LOG"
fi
echo " log        : $LOG"
exit $status