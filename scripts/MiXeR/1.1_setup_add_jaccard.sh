
# I prefer the Jaccard to the pre-existing DICE coefficient
# But I need to alter MIXER to produce the coefficient and its standard error. 
#!/usr/bin/env bash
set -euo pipefail

# Always run from project root (…/gsa), even if invoked elsewhere
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

MIXER_SIF="${MIXER_SIF:-mixer.sif}"
CODE_DIR="mixer_code"

mkdir -p "$CODE_DIR" scripts figures

# ----------------------------
# 1) Copy MiXeR python code out of the container (so we can edit it)
# ----------------------------
singularity exec --home "$ROOT:/home" "$MIXER_SIF" bash -lc "
  rm -rf /home/${CODE_DIR}/bivar_mixer /home/${CODE_DIR}/mixer_figures.py
  cp -r /tools/mixer/precimed/bivar_mixer /home/${CODE_DIR}/
  cp    /tools/mixer/precimed/mixer_figures.py /home/${CODE_DIR}/
"

# ----------------------------
# 2) Patch MiXeR code: add Jaccard as a derived statistic + export it
#    (NO confidence intervals, NO delta-method)
# ----------------------------
cat > scripts/patch_add_jaccard.py <<'PY'
from pathlib import Path
import re

root = Path("mixer_code")
utils_py   = root / "bivar_mixer" / "utils.py"
figures_py = root / "bivar_mixer" / "figures.py"

# ---- utils.py: add derived statistic 'jaccard' next to 'dice'
txt = utils_py.read_text()

dice_pat = re.compile(
    r"\('dice',\s*lambda x:\s*\(2\s*\*\s*x\._pi\[2\]\)\s*/\s*\(x\._pi\[0\]\s*\+\s*x\._pi\[1\]\s*\+\s*2\s*\*\s*x\._pi\[2\]\)\s*\),"
)
m = dice_pat.search(txt)
if not m:
    raise RuntimeError(f"Could not find dice definition in {utils_py}")

if "('jaccard'," not in txt:
    # J = pi12 / (pi1 + pi2 + pi12)
    insert = m.group(0) + "\n             ('jaccard', lambda x: (x._pi[2]) / (x._pi[0] + x._pi[1] + x._pi[2])),"
    txt = txt.replace(m.group(0), insert)

utils_py.write_text(txt)

# ---- figures.py: request jaccard in output (so it appears as "jaccard (mean)" / "jaccard (std)")
txt = figures_py.read_text()

txt = txt.replace(
    "keys = 'dice pi1 pi2 pi12 nc1@p9 nc2@p9 nc12@p9 rho_zero rho_beta rg fraction_concordant_within_shared'.split()",
    "keys = 'dice jaccard pi1 pi2 pi12 nc1@p9 nc2@p9 nc12@p9 rho_zero rho_beta rg fraction_concordant_within_shared'.split()"
)
txt = txt.replace(
    "keys = 'dice pi1 pi2 pi12 nc1@p9 nc2@p9 nc12@p9 rho_zero rho_beta rg_sig2_factor rg fraction_concordant_within_shared'.split()",
    "keys = 'dice jaccard pi1 pi2 pi12 nc1@p9 nc2@p9 nc12@p9 rho_zero rho_beta rg_sig2_factor rg fraction_concordant_within_shared'.split()"
)

# also include it in bivariate_keys lists (export table/CSV)
if "'jaccard'" not in txt:
    txt = txt.replace("'dice',", "'dice', 'jaccard',")

figures_py.write_text(txt)

print("Patched OK:")
print(" -", utils_py)
print(" -", figures_py)
PY

# Apply patch using container python
singularity exec --home "$ROOT:/home" "$MIXER_SIF" bash -lc \
  "python /home/scripts/patch_add_jaccard.py"

# ----------------------------
# 3) Re-run combine for ONE pair so the combined JSONs now include ci['jaccard']
#    (this is what fixes your KeyError)
# ----------------------------
singularity exec --home "$ROOT:/home" "$MIXER_SIF" bash -lc '
PYTHONPATH=/home/mixer_code:/tools/mixer/precimed \
python /home/mixer_code/mixer_figures.py combine \
  --json /home/out/AGREE_vs_IQ.fit.rep@.json \
  --out  /home/out/AGREE_vs_IQ.fit

PYTHONPATH=/home/mixer_code:/tools/mixer/precimed \
python /home/mixer_code/mixer_figures.py combine \
  --json /home/out/AGREE_vs_IQ.test.rep@.json \
  --out  /home/out/AGREE_vs_IQ.test
'

# ----------------------------
# 4) Now run "two" using the NEW combined JSONs
# ----------------------------
singularity exec --home "$ROOT:/home" "$MIXER_SIF" bash -lc '
PYTHONPATH=/home/mixer_code:/tools/mixer/precimed \
python /home/mixer_code/mixer_figures.py two \
  --json-fit  /home/out/AGREE_vs_IQ.fit.json \
  --json-test /home/out/AGREE_vs_IQ.test.json \
  --out       /home/figures/AGREE_vs_IQ \
  --trait1    AGREE \
  --trait2    IQ \
  --statistic mean std \
  --ext svg
'

# ----------------------------
# 5) Verify jaccard columns exist in the CSV header
# ----------------------------
echo "Looking for jaccard columns in CSV header:"
head -n 1 figures/AGREE_vs_IQ.csv | tr '\t' '\n' | grep -E '^jaccard \(mean\)$|^jaccard \(std\)$' || {
  echo "ERROR: jaccard columns not found in header."
  exit 2
}

echo "Success: jaccard is in the CSV."
