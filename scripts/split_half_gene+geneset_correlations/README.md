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

## Correlations and SEM

| Script | Purpose |
| --- | --- |
| [01_gene_scores.R](01_gene_scores.R) | Gene-score correlations and jackknife uncertainty. |
| [02_gene_sets.R](02_gene_sets.R) | Gene-set correlations, reliability regressions, and the paired figure. |
| [03_closed_form_correlations.R](03_closed_form_correlations.R) | Closed-form latent correlations for all 15 trait pairs in gene scores, pooled gene sets, and each gene-set collection; group means and their difference, with direct block-jackknife SEs. |
| [04_sem.R](04_sem.R) | Correlated-trait and hierarchical SEMs, fit diagnostics, and jackknife uncertainty. |

Script 03 writes its tables and `closed_form_correlations.xlsx` workbook to
`output/closed_form_correlations/`.

It reads both pooled `inputs.rds` files and `output/gene_sets/collection_inputs.rds`.
The `collection` column distinguishes `all_genes`, `pooled`, and individual
gene-set collections. Group means require all 10 Big Five pairs or all five
IQ–Big Five pairs; undefined pairs are not dropped. SEs require every one of
the 200 deletions to be defined. Defined-pair and valid-deletion counts are
reported, and out-of-range correlations and means are retained and flagged.
Each group mean and their difference receives a two-sided approximate normal
p-value against zero when its jackknife SE is positive and finite; these
p-values are not adjusted for multiple comparisons. The collection IDs match
the 22 collections in the reliability figure, which labels only eight of them
and displays observed rather than reliability-corrected correlations.

## Additional scripts

| Script | Purpose |
| --- | --- |
| [5_gene_correlations.sh](5_gene_correlations.sh) | Calculates Pearson gene score correlations for the 200 jackknife iterations. |
| [6_gene_corr_no_jk.sh](6_gene_corr_no_jk.sh) | Calculates Pearson gene score correlations with no excluded genes. |
| [gene_corr_no_jk.R](gene_corr_no_jk.R) | This script is run by 6_gene_corr_no_jk.sh. |
| [gene_zstat_correlation.R](gene_zstat_correlation.R) | This script is run by 5_gene_correlations.sh. |
| [mean_N.R](mean_N.R) | Calculates the mean sample size across all SNPs for a sumstats file. It also returns other quantities like minimum, maximum, and median. |
