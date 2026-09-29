#!/usr/bin/env Rscript

# ==============================================================================
# Setup and collection exclusions
# ==============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, data.table, fixest, ggrepel, ggthemes, writexl)

dir.create("output/gene_sets", recursive = TRUE, showWarnings = FALSE)
traits <- c("agree", "consc", "extra", "neurot", "open", "iq")
indicators <- c(
  "agr1", "agr2", "con1", "con2", "ext1", "ext2",
  "neu1", "neu2", "ope1", "ope2", "iq_female", "iq_male"
)
trait_files <- c(
  "ReGPC_agr_half_one_no23_dir", "ReGPC_agr_half_two_no23_dir",
  "ReGPC_con_half_one_no23_dir", "ReGPC_con_half_two_no23_dir",
  "ReGPC_ext_half_one_no23_dir", "ReGPC_ext_half_two_no23_dir",
  "ReGPC_neu_half_one_no23_dir", "ReGPC_neu_half_two_no23_dir",
  "ReGPC_ope_half_one_no23_dir", "ReGPC_ope_half_two_no23_dir",
  "iq_female_dir", "iq_male_dir"
)
archive <- "data/split_half_geneset_output.tar.gz"
jackknife_dir <- "data/split_half_output/output"
block_files <- sprintf("data/excluded_genes/exclude_genes_%03d.txt", 1:200)
pairs <- t(combn(1:12, 2))
pair_names <- paste(indicators[pairs[, 1]], indicators[pairs[, 2]], sep = "__")

# Preserve the exclusions in the audited input bundle. Changing this list
# changes the analysis sample and should be a separate scientific decision.
excluded_collections <- "neural_brain"
excluded_sets <- "GOBP_NEGATIVE_REGULATION_OF_INTRACELLULAR_LIPID_TRANSPORT"
collections <- list.dirs(jackknife_dir, recursive = FALSE, full.names = FALSE) |>
  setdiff(excluded_collections) |>
  sort()
stopifnot(length(collections) > 0, all(file.exists(block_files)))

# ==============================================================================
# Read and align MAGMA enrichment coefficients
# ==============================================================================

# Each collection uses sets measured in all 12 GWAS. Gene-score missingness in
# script 01 is handled pairwise; enrichment inputs here already have full rows.
read_collection <- function(paths, id_column) {
  results <- map(paths, \(path) {
    x <- fread(path, skip = "VARIABLE", select = c(id_column, "BETA_STD"),
               showProgress = FALSE) |>
      as_tibble() |>
      rename(gene_set = all_of(id_column)) |>
      filter(!gene_set %in% excluded_sets)
    stopifnot(!anyDuplicated(x$gene_set), !anyNA(x$gene_set), all(is.finite(x$BETA_STD)))
    x
  })
  common <- reduce(map(results, "gene_set"), intersect)
  stopifnot(length(common) >= 2)
  x <- vapply(results, \(x) x$BETA_STD[match(common, x$gene_set)], numeric(length(common)))
  dimnames(x) <- list(common, indicators)
  x
}

pool_sets <- function(collection_data) {
  x <- do.call(rbind, collection_data)
  ids <- rownames(x)
  duplicates <- unique(ids[duplicated(ids)])
  for (id in duplicates) {
    values <- x[ids == id, , drop = FALSE]
    if (any(apply(values, 2, \(z) diff(range(z))) > 1e-12)) {
      stop("Conflicting pooled gene set: ", id)
    }
  }
  x[!duplicated(ids), , drop = FALSE]
}

# The archive is authoritative; do not silently substitute an extracted copy.
extracted <- tempfile("full_jack_gene_sets_")
dir.create(extracted)
members <- untar(archive, list = TRUE)
status <- untar(archive, files = members[str_detect(members, "\\.gsa\\.out$")], exdir = extracted)
stopifnot(status == 0)
full_paths <- map(collections, \(collection) {
  file.path(extracted, "geneset_output", collection, paste0(trait_files, ".gsa.out"))
}) |> set_names(collections)
id_columns <- map_chr(full_paths, \(paths) {
  columns <- names(fread(paths[1], skip = "VARIABLE", nrows = 0))
  if ("FULL_NAME" %in% columns) "FULL_NAME" else "VARIABLE"
})
data <- map2(full_paths, id_columns, read_collection)
data$pooled <- pool_sets(data)
unlink(extracted, recursive = TRUE)
counts <- imap_dfr(data, \(x, collection) tibble(collection = collection, n_sets = nrow(x)))
print(counts)
correlations <- map(data, cor)
stopifnot(all(map_lgl(correlations, \(x) all(is.finite(x)))))

# ==============================================================================
# Recalculate correlations in all 200 MAGMA deletion samples
# ==============================================================================

# These MAGMA files must use the same exclusion lists as script 01. Their
# correspondence was audited on September 28; retain hashes of both sources.
paths <- crossing(collection = collections, trait_file = trait_files, block = 1:200) |>
  mutate(path = file.path(jackknife_dir, collection, paste0("jack_", trait_file),
                          sprintf("%s_jk%03d.gsa.out", trait_file, block)))
stopifnot(all(file.exists(paths$path)))
jackknife <- map(data, \(x) matrix(NA_real_, 200, 66, dimnames = list(NULL, pair_names)))

for (block in 1:200) {
  deleted <- map2(collections, id_columns, \(collection, id_column) {
    files <- file.path(jackknife_dir, collection, paste0("jack_", trait_files),
                       sprintf("%s_jk%03d.gsa.out", trait_files, block))
    read_collection(files, id_column)
  }) |> set_names(collections)
  deleted$pooled <- pool_sets(deleted)
  for (collection in names(deleted)) {
    jackknife[[collection]][block, ] <- cor(deleted[[collection]])[pairs]
  }
  if (block %% 10 == 0) message("Gene-set correlations: ", block, "/200 blocks")
}
stopifnot(all(map_lgl(jackknife, \(x) all(is.finite(x)))))
sampling_vcov <- map(jackknife, \(x) {
  centered <- sweep(x, 2, colMeans(x), "-")
  (199 / 200) * crossprod(centered)
})
sources <- tibble(path = c(archive, block_files, paths$path)) |>
  mutate(md5 = unname(tools::md5sum(path)))
inputs <- list(
  correlation = correlations$pooled, jackknife = jackknife$pooled, vcov = sampling_vcov$pooled,
  traits = traits, indicators = indicators, pairs = pairs,
  block_hashes = unname(tools::md5sum(block_files)), sources = sources
)
saveRDS(inputs, "output/gene_sets/inputs.rds")
saveRDS(list(correlations = correlations, jackknife = jackknife, vcov = sampling_vcov),
        "output/gene_sets/collection_inputs.rds")

# ==============================================================================
# Does reliability predict cross-trait biological overlap?
# ==============================================================================

same_trait <- ceiling(pairs[, 1] / 2) == ceiling(pairs[, 2] / 2)
cross_half <- pairs[, 1] %% 2 != pairs[, 2] %% 2
iq_pair <- pairs[, 2] > 10

# Each trait pair contributes its two cross-half correlations equally. Same-half
# cross-trait correlations are excluded because of overlapping Big Five samples.
summarize_collection <- function(x) {
  tibble(
    reliability = mean(x[same_trait]),
    all_pairs = mean(x[!same_trait & cross_half]),
    iq_big_five = mean(x[!same_trait & cross_half & iq_pair]),
    big_five = mean(x[!same_trait & cross_half & !iq_pair])
  )
}
collection_summary <- map_dfr(collections, \(collection) {
  summarize_collection(correlations[[collection]][pairs]) |>
    mutate(collection = collection, .before = 1)
})
deleted_summary <- map_dfr(1:200, \(block) {
  map_dfr(collections, \(collection) {
    summarize_collection(jackknife[[collection]][block, ]) |>
      mutate(collection = collection, block = block, .before = 1)
  })
})

# OLS across collections; its uncertainty comes from deleting genes and refitting,
# not from treating overlapping collections as independent observations.
models <- list(
  all_pairs = feols(all_pairs ~ reliability, data = collection_summary),
  iq_big_five = feols(iq_big_five ~ reliability, data = collection_summary),
  big_five = feols(big_five ~ reliability, data = collection_summary)
)
deleted_models <- map(1:200, \(block) {
  deleted_data <- filter(deleted_summary, .data$block == .env$block)
  list(
    all_pairs = feols(all_pairs ~ reliability, data = deleted_data),
    iq_big_five = feols(iq_big_five ~ reliability, data = deleted_data),
    big_five = feols(big_five ~ reliability, data = deleted_data)
  )
})
regression <- imap_dfr(models, \(model, series) {
  estimates <- do.call(rbind, map(deleted_models, \(x) coef(x[[series]])))
  centered <- sweep(estimates, 2, colMeans(estimates), "-")
  tibble(series = series, term = names(coef(model)), estimate = unname(coef(model)),
         standard_error = sqrt((199 / 200) * colSums(centered^2)))
}) |>
  mutate(lower = estimate - qnorm(.975) * standard_error,
         upper = estimate + qnorm(.975) * standard_error)

uncertainty <- deleted_summary |>
  pivot_longer(c(reliability, all_pairs, iq_big_five, big_five), names_to = "quantity") |>
  group_by(collection, quantity) |>
  summarise(standard_error = sqrt((199 / 200) * sum((value - mean(value))^2)), .groups = "drop")
collection_results <- collection_summary |>
  pivot_longer(-collection, names_to = "quantity", values_to = "estimate") |>
  left_join(uncertainty, by = c("collection", "quantity"))

collection_labels <- c(
  GO_cell_component = "GO cellular component", akingbuwa_gs = "Akingbuwa",
  brainspan = "BrainSpan", computational_perturbation_signatures = "Computational perturbations",
  gtex = "GTEx tissues", hallmark_processes = "Hallmark processes",
  kegg_medicus_pathways = "KEGG Medicus", microRNA_targets = "microRNA targets"
)
figure <- collection_summary |>
  pivot_longer(c(iq_big_five, big_five), names_to = "series", values_to = "correlation") |>
  arrange(match(series, c("iq_big_five", "big_five"))) |>
  mutate(series = factor(recode(series, iq_big_five = "IQ - Big Five",
                               big_five = "Big Five - Big Five"),
                         levels = c("Big Five - Big Five", "IQ - Big Five"))) |>
  ggplot(aes(reliability, correlation)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#A3A3A3", linewidth = .65) +
  geom_segment(data = collection_summary,
               aes(x = reliability, xend = reliability, y = iq_big_five, yend = big_five),
               inherit.aes = FALSE, colour = "#A3A3A3", linewidth = .65) +
  geom_point(aes(colour = series), size = 4.4) +
  geom_label_repel(
    data = collection_summary |>
      mutate(midpoint = (iq_big_five + big_five) / 2) |>
      pivot_longer(c(midpoint, iq_big_five, big_five), names_to = "position", values_to = "y") |>
      arrange(match(position, c("midpoint", "iq_big_five", "big_five"))) |>
      mutate(label = if_else(position == "midpoint" & collection %in% names(collection_labels),
                             collection_labels[collection], "")),
    aes(x = reliability, y = y, label = label),
    inherit.aes = FALSE, size = 4.4, seed = 20260729,
    box.padding = .7, label.padding = grid::unit(.18, "lines"), point.padding = 0,
    label.r = grid::unit(.2, "lines"), linewidth = .22, min.segment.length = 0,
    force = 4, force_pull = .25, max.time = 8, max.iter = 100000,
    max.overlaps = Inf, direction = "both", fill = scales::alpha("#EAF7FC", .92),
    colour = "#30414D", segment.color = scales::alpha("#71838F", .82),
    segment.size = .45, show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Big Five - Big Five" = "#2C7FB8", "IQ - Big Five" = "#D95F02"),
                      drop = FALSE) +
  labs(x = "Average same-trait gene set correlation",
       y = "Average cross-trait gene set correlation", colour = NULL) +
  scale_x_continuous(expand = expansion(mult = c(.10, .06))) +
  scale_y_continuous(expand = expansion(mult = c(.12, .06))) +
  coord_fixed(ratio = 1) +
  ggthemes::theme_foundation(base_size = 17, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = NA),
    axis.title = element_text(face = "bold"), axis.line = element_line(colour = "black"),
    panel.grid.major = element_line(colour = "#f0f0f0"), panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = scales::alpha("white", .88), colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.position = "inside", legend.position.inside = c(.03, .97),
    legend.justification.inside = c(0, 1), legend.direction = "vertical",
    legend.title.position = "top", legend.margin = margin(5, 7, 5, 7)
  )
ggsave("output/gene_sets/reliability_overlap.pdf", figure, width = 8.2, height = 7.4)
ggsave("output/gene_sets/reliability_overlap.png", figure, width = 8.2, height = 7.4, dpi = 240)

# ==============================================================================
# Save paper tables and the collection-level results
# ==============================================================================

correlation_table <- imap_dfr(correlations, \(x, collection) {
  se <- sqrt(diag(sampling_vcov[[collection]]))
  tibble(collection = collection, indicator_1 = indicators[pairs[, 1]],
         indicator_2 = indicators[pairs[, 2]], correlation = x[pairs],
         standard_error = se)
})
write_csv(counts, "output/gene_sets/set_counts.csv")
write_csv(correlation_table, "output/gene_sets/correlations.csv")
write_csv(collection_results, "output/gene_sets/collection_summary.csv")
write_csv(regression, "output/gene_sets/reliability_regression.csv")
write_csv(deleted_summary, "output/gene_sets/jackknife_collection_summary.csv")
write_csv(sources, "output/gene_sets/sources.csv")
write_xlsx(list(correlations = correlation_table, collection_summary = collection_results,
                regression = regression, set_counts = counts), "output/gene_sets/gene_sets.xlsx")
writeLines(capture.output(sessionInfo()), "output/gene_sets/session_info.txt")
message("Saved enrichment correlations, reliability regressions and the paper figure.")
