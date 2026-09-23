#!/bin/bash
#SBATCH --job-name=participation_bias_correction

#SBATCH --mem=20G          # memory per node

#SBATCH --time=10:00:00    # max runtime hh:mm:ss

#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)

#SBATCH --output=pbias.%j.out

#SBATCH --error=pbias.%j.err


# Load R module
module load freetype libpng libtiff libjpeg libwebp pkg-config
module load R/4.4.2-openblas-rocky8
conda activate pbias_env

#Define phenotype set
PHENOLIST='open consc extra agree neurot iq'

#Run LDSC for all phenotypes with participation
for PHENO in $PHENOLIST; do

	cd /users/0/buona008/andrew/personality_iq_project/participation_bias_correction/

	mkdir -p results_${PHENO}

	cd results_${PHENO}

	python2.7 ../ldsc_jackknife/ldsc.py --rg /users/0/buona008/andrew/personality_iq_project/participation_bias_correction/participationprimary_benonisdottir2023.sumstats.gz,/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/${PHENO}_schwaba2025.sumstats.gz --ref-ld ../ldsc_jackknife/UKBB.EUR --w-ld ../ldsc_jackknife/UKBB.EUR --intercept-gencov 0,0 --out res_rg

done

    python2.7 ../ldsc_jackknife/ldsc.py --rg /users/0/buona008/andrew/personality_iq_project/participation_bias_correction/participationprimary_benonisdottir2023.sumstats.gz,/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/iq_savage2018.sumstats.gz --ref-ld ../ldsc_jackknife/UKBB.EUR --w-ld ../ldsc_jackknife/UKBB.EUR --intercept-gencov 0,0 --out res_rg
mkdir -p ../results_ea
cd ../results_ea
   python2.7 ../ldsc_jackknife/ldsc.py --rg /users/0/buona008/andrew/personality_iq_project/participation_bias_correction/participationprimary_benonisdottir2023.sumstats.gz,/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/ea_okbay2022.sumstats.gz --ref-ld ../ldsc_jackknife/UKBB.EUR --w-ld ../ldsc_jackknife/UKBB.EUR --intercept-gencov 0,0 --out res_rg


#Find pairwise LDSC correlations between big five and intelligence
PHENOLIST='open consc extra agree neurot'
for PHENO in $PHENOLIST; do

cd /users/0/buona008/andrew/personality_iq_project/participation_bias_correction

mkdir -p results_${PHENO}_iq

cd results_${PHENO}_iq

python2.7 /users/0/buona008/andrew/personality_iq_project/participation_bias_correction/ldsc_jackknife/ldsc.py --rg /users/0/buona008/andrew/participation_bias_correction/${PHENO}_schwaba2025.sumstats.gz,/users/0/buona008/andrew/participation_bias_correction/iq_savage2018.sumstats.gz --ref-ld ../ldsc_jackknife/UKBB.EUR --w-ld ../ldsc_jackknife/UKBB.EUR --out res_rg

done

for PHENO in $PHENOLIST; do

cd /users/0/buona008/andrew/participation_bias_correction

mkdir -p results_${PHENO}_ea

cd results_${PHENO}_ea

python2.7 /users/0/buona008/andrew/participation_bias_correction/ldsc_jackknife/ldsc.py --rg /users/0/buona008/andrew/participation_bias_correction/${PHENO}_schwaba2025.sumstats.gz,/users/0/buona008/andrew/participation_bias_correction/ea_okbay2022.sumstats.gz --ref-ld ../ldsc_jackknife/UKBB.EUR --w-ld ../ldsc_jackknife/UKBB.EUR --out res_rg

done


# Run the R script
cd /users/0/buona008/andrew/participation_bias_correction/R_scripts
Rscript participation_bias_correction.R


