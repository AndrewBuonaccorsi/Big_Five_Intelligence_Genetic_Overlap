# Split-half gene and gene-set correlations

## MAGMA preparation

| Script | Purpose |
| --- | --- |
| [00_add_Pvalues.sh](00_add_Pvalues.sh) | Add a P column to munged sumstats based on their Z-score. |
| [add_P.R](add_P.R) | This script is run by 00_add_Pvalues.sh |
| [01_magma_annotate.sh](01_magma_annotate.sh) | Use NCBI build 37.3 gene locations to map SNPs to genes. |
| [02_magma_gene_analysis.sh](02_magma_gene_analysis.sh) | Create gene scores from summary statistics and snp-gene map. |
| [03_magma_merge.sh](03_magma_merge.sh) | Merge all chromosomes of the above. |
| [04_magma_geneset.sh](04_magma_geneset.sh) | Get gene set enrichments from merged gene scores. |
| [3_magma_jack.sh](3_magma_jack.sh) | Get gene set enrichmetns from merged gene scores for 200 jackknife iterations. Excluded gene blocks can be found in the data directory. |

## Current analyses

| Script | Purpose |
| --- | --- |
| [01_gene_scores.R](01_gene_scores.R) | Gene-score correlations and jackknife uncertainty. |
| [02_gene_sets.R](02_gene_sets.R) | Gene-set correlations, reliability regressions, and the paired figure. |
| [03_analytic_correlations.R](03_analytic_correlations.R) | Analytic latent correlations and comparisons between trait groups. |
| [04_sem.R](04_sem.R) | Correlated-trait and hierarchical SEMs, fit diagnostics, and jackknife uncertainty. |

## Additional scripts

| Script | Purpose |
| --- | --- |
| [5_gene_correlations.sh](5_gene_correlations.sh) | Calculates Pearson gene score correlations for the 200 jackknife iterations. |
| [6_gene_corr_no_jk.sh](6_gene_corr_no_jk.sh) | Calculates Pearson gene score correlations with no excluded genes. |
| [gene_corr_no_jk.R](gene_corr_no_jk.R) | This script is run by 6_gene_corr_no_jk.sh. |
| [gene_zstat_correlation.R](gene_zstat_correlation.R) | This script is run by 5_gene_correlations.sh. |
| [mean_N.R](mean_N.R) | Calculates the mean sample size across all SNPs for a sumstats file. It also returns other quantities like minimum, maximum, and median. |
