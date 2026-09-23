#!/bin/bash

#SBATCH --job-name=magma_merge

#SBATCH --mem=20G          # memory per node

#SBATCH --time=10:00:00    # max runtime hh:mm:ss

#SBATCH --cpus-per-task=1  # number of CPUs (adjust as needed)

#SBATCH --output=magma_merge.%j.out

#SBATCH --error=magma_merge.%j.err

#module load R/4.4.2-openblas-rocky8


#find /users/0/buona008/andrew/all_phenotypes_enrich_corr/gene_output -type f -name "*.chr*.batch*_chr.genes.*" \
 #   -exec bash -c '
  #      for f; do
   #         new=$(sed -E "s/\.chr[0-9]+//" <<< "$f")
    #        mv "$f" "$new"
     #   done
    #' bash {} +


mkdir -p merged_output

for folder in gene_output/*_dir; do

    echo "Processing folder: $folder"

    BASE="${folder##*/}"
    BASE="${BASE%_dir}"

    echo "BASE = $BASE"

    mkdir -p merged_output/"${BASE}_dir"

/users/0/buona008/andrew/magma --merge /users/0/buona008/andrew/all_phenotypes_enrich_corr/gene_output/"${BASE}_dir"/"${BASE}" --out merged_output/"${BASE}_dir"

done

wait
