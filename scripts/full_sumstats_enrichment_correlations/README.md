# Full-sample gene-set enrichment correlations

| Script | Purpose |
| --- | --- |
| [3_magma_jack.sh](3_magma_jack.sh) | Performs gene set enrichment analysis, 200 jackknife iterations. Needs raw gene scores and an excluded genes file, which is in our data directory. |
| [4_enrich_corr_jack.sh](4_enrich_corr_jack.sh) | Calculates enrichment correlations over the 200 jackknife iterations. |
| [04_magma_geneset.sh](04_magma_geneset.sh) | Equivalent to 3_magma_jack.sh but does not exclude any genes. |
| [7_non_jk_enrich_corr.sh](7_non_jk_enrich_corr.sh) | Equivalent to 4_enrich_corr_jack.sh but uses gene set enrichments that did not exclude any genes. |
| [ranker.R](ranker.R) | Takes .gsa.out files for each trait and returns a ranking of their mean z-transformed MAGMA enrichments. |
| [8_ranker.sh](8_ranker.sh) | Runs ranker.R. |
| [unweighted_enrich_corr_njk.R](unweighted_enrich_corr_njk.R) | This is the script run by 7_non_jk_enrich_corr.sh. |
| [unweighted_enrich_corr.R](unweighted_enrich_corr.R) | This is the script run by 4_enrich_corr_jack.sh. |
