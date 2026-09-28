#!/usr/bin/env Rscript

# ==============================================================================
# Setup
# ==============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, data.table, writexl)

dir.create("output/gene_scores", recursive = TRUE, showWarnings = FALSE)
traits <- c("agree", "consc", "extra", "neurot", "open", "iq")
indicators <- c(
  "agr1", "agr2", "con1", "con2", "ext1", "ext2",
  "neu1", "neu2", "ope1", "ope2", "iq_female", "iq_male"
)
files <- file.path("data/gene_enrichments", paste0(c(
  "ReGPC_agr_half_one_no23_dir", "ReGPC_agr_half_two_no23_dir",
  "ReGPC_con_half_one_no23_dir", "ReGPC_con_half_two_no23_dir",
  "ReGPC_ext_half_one_no23_dir", "ReGPC_ext_half_two_no23_dir",
  "ReGPC_neu_half_one_no23_dir", "ReGPC_neu_half_two_no23_dir",
  "ReGPC_ope_half_one_no23_dir", "ReGPC_ope_half_two_no23_dir",
  "iq_female_dir", "iq_male_dir"
), ".genes.raw"))
block_files <- sprintf("data/excluded_genes/exclude_genes_%03d.txt", 1:200)
pairs <- t(combn(1:12, 2))
pair_names <- paste(indicators[pairs[, 1]], indicators[pairs[, 2]], sep = "__")

# ==============================================================================
# Read gene scores and the shared MAGMA deletion blocks
# ==============================================================================

# MAGMA raw files have a variable number of LD columns after column 9 (ZSTAT).
# Extract only the ID, position and score rather than parsing those LD columns.
raw <- map2(files, indicators, \(path, indicator) {
  x <- fread(
    cmd = paste("awk '!/^#/ {print $1, $2, $3, $4, $9}'", shQuote(path)),
    header = FALSE, col.names = c("gene", "chr", "start", "end", "z"),
    showProgress = FALSE
  ) |>
    as_tibble() |>
    mutate(gene = as.character(gene), indicator = indicator)
  stopifnot(!anyDuplicated(x$gene), all(is.finite(x$z)))
  x
}) |>
  list_rbind()

positions <- raw |> distinct(gene, chr, start, end)
stopifnot(!anyDuplicated(positions$gene))
blocks <- map2(block_files, 1:200, \(path, block) {
  tibble(gene = read_lines(path), block = block)
}) |>
  list_rbind()
stopifnot(!anyDuplicated(blocks$gene), all(str_detect(blocks$gene, "^[0-9]+$")))

# Keep the union of genes. Each Pearson correlation uses its own complete pairs.
# Every analyzed gene must occur in exactly one of the shared deletion lists.
data <- raw |>
  select(gene, indicator, z) |>
  pivot_wider(names_from = indicator, values_from = z) |>
  left_join(positions, by = "gene") |>
  left_join(blocks, by = "gene") |>
  arrange(chr, start, end, gene) |>
  select(gene, chr, start, end, block, all_of(indicators))
stopifnot(!anyNA(data$block), setequal(data$block, 1:200))
print(count(data, block))
print(summarise(data, genes = n(), complete_genes = sum(if_all(all_of(indicators), ~ !is.na(.x)))))
scores <- as.matrix(select(data, all_of(indicators)))

# ==============================================================================
# Observed correlations and 200 delete-one-block estimates
# ==============================================================================

correlation <- cor(scores, use = "pairwise.complete.obs", method = "pearson")
jackknife <- do.call(rbind, map(1:200, \(block) {
  cor(scores[data$block != block, ], use = "pairwise.complete.obs")[pairs]
}))
colnames(jackknife) <- pair_names
stopifnot(all(is.finite(correlation)), all(is.finite(jackknife)))
centered <- sweep(jackknife, 2, colMeans(jackknife), "-")
sampling_vcov <- (199 / 200) * crossprod(centered)

correlations <- tibble(
  indicator_1 = indicators[pairs[, 1]],
  indicator_2 = indicators[pairs[, 2]],
  trait_1 = rep(traits, each = 2)[pairs[, 1]],
  trait_2 = rep(traits, each = 2)[pairs[, 2]],
  n_genes = map_int(seq_len(nrow(pairs)), \(i) {
    sum(complete.cases(scores[, pairs[i, ]]))
  }),
  correlation = correlation[pairs],
  standard_error = sqrt(diag(sampling_vcov))
)
reliability <- correlations |> filter(trait_1 == trait_2)

# ==============================================================================
# Save analysis inputs and report tables
# ==============================================================================

sources <- tibble(path = c(files, block_files)) |>
  mutate(md5 = unname(tools::md5sum(path)))
inputs <- list(
  correlation = correlation, jackknife = jackknife, vcov = sampling_vcov,
  traits = traits, indicators = indicators, pairs = pairs,
  block_hashes = unname(tools::md5sum(block_files)), sources = sources
)
saveRDS(inputs, "output/gene_scores/inputs.rds")
write_csv(data, "output/gene_scores/gene_scores.csv")
write_csv(correlations, "output/gene_scores/correlations.csv")
write_csv(reliability, "output/gene_scores/reliability.csv")
write_csv(count(data, block, name = "n_genes"), "output/gene_scores/block_counts.csv")
write_csv(sources, "output/gene_scores/sources.csv")
write_csv(as_tibble(jackknife) |> mutate(block = 1:200, .before = 1),
          "output/gene_scores/jackknife_correlations.csv")
write.csv(correlation, "output/gene_scores/correlation_matrix.csv")
write.csv(sampling_vcov, "output/gene_scores/sampling_covariance.csv")
write_xlsx(list(correlations = correlations, reliability = reliability,
                block_counts = count(data, block, name = "n_genes")),
           "output/gene_scores/gene_scores.xlsx")
writeLines(capture.output(sessionInfo()), "output/gene_scores/session_info.txt")
message("Saved pairwise-complete gene-score correlations and all 200 deletions.")
