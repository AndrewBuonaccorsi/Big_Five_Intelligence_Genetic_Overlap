# Go to your MiXeR working directory
cd /projects/standard/leej5/edwa0506/andrew_project/gsa

module load singularity

# Path to your local MiXeR Singularity image
export MIXER_SIF=/projects/standard/leej5/edwa0506/andrew_project/gsa/mixer.sif

# Wrapper to call MiXeR inside the container
export MIXER_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

# Wrapper to call MiXeR figures inside the container
export MIXER_FIG="singularity exec --home $PWD:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer_figures.py"

TRAITS=(
  OPEN
  CONSC
  EXTRA
  AGREE
  NEURO
  IQ
  HEIGHT
)


# ----------------------------
# 1) Combine replicates per trait
# ----------------------------

for TRAIT in "${TRAITS[@]}"; do
  $MIXER_FIG combine \
    --json "out/${TRAIT}.fit.rep@.json" \
    --out  "out/${TRAIT}.fit"

  $MIXER_FIG combine \
    --json "out/${TRAIT}.test.rep@.json" \
    --out  "out/${TRAIT}.test"
done


# ----------------------------
# 2) Combine replicates of pairws
# ----------------------------

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
  "OPEN CONSC"
  "OPEN EXTRA"
  "OPEN AGREE"
  "OPEN NEURO"
  "CONSC EXTRA"
  "CONSC AGREE"
  "CONSC NEURO"
  "EXTRA AGREE"
  "EXTRA NEURO"
  "AGREE NEURO"
)


for P in "${PAIRS[@]}"; do
  set -- $P
  T1="$1"
  T2="$2"
  PAIR="${T1}_vs_${T2}"

  $MIXER_FIG combine \
    --json "out/${PAIR}.fit.rep@.json" \
    --out  "out/${PAIR}.fit"

  $MIXER_FIG combine \
    --json "out/${PAIR}.test.rep@.json" \
    --out  "out/${PAIR}.test"
done

# ----------------------------
# 2) One figure that includes ALL traits at once
#    (uses all per-trait combined JSONs)
# ----------------------------
$MIXER_FIG one \
  --json out/OPEN.fit.json out/CONSC.fit.json out/EXTRA.fit.json out/AGREE.fit.json out/NEURO.fit.json out/IQ.fit.json out/HEIGHT.fit.json \
  --out  figures/personality.fit \
  --trait1 OPEN CONSC EXTRA AGREE NEURO IQ HEIGHT \
  --statistic mean std \
  --ext svg

$MIXER_FIG one \
  --json out/OPEN.test.json out/CONSC.test.json out/EXTRA.test.json out/AGREE.test.json out/NEURO.test.json out/IQ.test.json out/HEIGHT.test.json \
  --out  figures/personality.test \
  --trait1 OPEN CONSC EXTRA AGREE NEURO IQ HEIGHT \
  --statistic mean std \
  --ext svg
  
  
  # --------------------------------------------------
# Two-trait (bivariate) plot (figures → figures/)
# --------------------------------------------------

for P in "${PAIRS[@]}"; do
  set -- $P
  T1="$1"
  T2="$2"
  PAIR="${T1}_vs_${T2}"

  $MIXER_FIG two \
    --json-fit  "out/${PAIR}.fit.json" \
    --json-test "out/${PAIR}.test.json" \
    --out       "figures/${PAIR}" \
    --trait1    "${T1}" \
    --trait2    "${T2}" \
    --statistic mean std \
    --ext svg
done