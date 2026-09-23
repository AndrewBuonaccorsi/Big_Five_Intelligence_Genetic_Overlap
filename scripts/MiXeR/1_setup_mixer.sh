cd /projects/standard/leej5/edwa0506/andrew_project/gsa

module load singularity
module load git

#-------------------------------------------------------------------------------
# Setup MIXER
#-------------------------------------------------------------------------------

# Get ORAS and put it in the bin
wget https://github.com/oras-project/oras/releases/download/v1.1.0/oras_1.1.0_linux_amd64.tar.gz
tar -xzf oras_1.1.0_linux_amd64.tar.gz
mkdir -p ~/.local/bin
mv oras ~/.local/bin/

# (Temporarily) Add it to my path so I can call oras
export PATH="$HOME/.local/bin:$PATH"
which oras   # sanity check


# Pull the mixer singularity image
oras pull ghcr.io/precimed/gsa-mixer_sif:2.2.1
mv gsa-mixer.sif mixer.sif

# Point MIXER_SIF to the image in this directory
export MIXER_SIF=/projects/standard/leej5/edwa0506/andrew_project/gsa/mixer.sif

# Command wrapper to call MiXeR
# (pwd is expanded *now*, so --home maps this gsa dir into /home in the container)
export MIXER_PY="singularity exec --home $(pwd):/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"


#-------------------------------------------------------------------------------
# Download reference files
#-------------------------------------------------------------------------------

#git clone https://github.com/comorment/mixer.git
#mv mixer/reference .
#rm -rf mixer

#mv reference .
#rm -rf reference

# Github didn't work use dropbox
# https://www.dropbox.com/scl/fo/y5yl2bd5mgplsjwwzsx77/AIFIhSJkzJTFIYhR95TwRVc?rlkey=eydtbzwva5294snzgf6lz0g5f&e=1&st=5g7d4jiv&dl=0
# Download to mac then uploaded with globus

ZIP="reference.zip"
OUT="reference"

mkdir -p "$OUT"

# 1) List all files that SHOULD be in the zip
UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE unzip -Z1 "$ZIP" \
  | grep -v '/$' \
  | sort > zip_files.txt

# 2) Extract everything that can be extracted
UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE unzip -n "$ZIP" -d "$OUT" || true

# 3) List files that were ACTUALLY extracted
find "$OUT" -type f | sed "s|^$OUT/||" | sort > extracted_files.txt

# 4) Files in zip but NOT extracted → need reupload
comm -23 zip_files.txt extracted_files.txt > missing_files.txt

echo "Missing files saved to: missing_files.txt"

# Check it's all ther
ls reference




#-------------------------------------------------------------------------------
# Copy my sumstats (symlink didn't work)
#-------------------------------------------------------------------------------

cd /projects/standard/leej5/edwa0506/andrew_project/gsa

mkdir -p sumstats

cp /projects/standard/leej5/edwa0506/summary_statistics_collection/private_sumstats/cleaned_relig_sumstats_no_attendance.sumstats.gz \
      summary_statistics/cleaned_relig_sumstats_no_attendance.sumstats.gz

cp /projects/standard/leej5/edwa0506/summary_statistics_collection/private_sumstats/ebb_spiritual_person.sumstats.gz \
      summary_statistics/ebb_spiritual_person.sumstats.gz
      




