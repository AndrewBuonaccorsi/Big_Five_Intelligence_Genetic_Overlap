#!/usr/bin/env Rscript

# Closed-form latent correlations from saved gene-score and gene-set inputs.

# ==============================================================================
# Load the two saved correlation inputs
# ==============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, writexl)

dir.create("output/closed_form_correlations", recursive = TRUE, showWarnings = FALSE)
paths <- c(gene_scores = "output/gene_scores/inputs.rds", gene_sets = "output/gene_sets/inputs.rds")
inputs <- map(paths, readRDS)
stopifnot(identical(inputs$gene_scores$block_hashes, inputs$gene_sets$block_hashes),
          identical(inputs$gene_scores$indicators, inputs$gene_sets$indicators))
traits <- inputs$gene_scores$traits
trait_pairs <- t(combn(1:6, 2))
pair_names <- paste(traits[trait_pairs[, 1]], traits[trait_pairs[, 2]], sep = "__")
iq_pair <- trait_pairs[, 2] == 6

# ==============================================================================
# Closed-form latent correlations and group contrasts
# ==============================================================================

closed_form_correlations <- function(r) {
  a <- trait_pairs[, 1]
  b <- trait_pairs[, 2]
  reliability_a <- r[cbind(2 * a - 1, 2 * a)]
  reliability_b <- r[cbind(2 * b - 1, 2 * b)]
  cross_12 <- r[cbind(2 * a - 1, 2 * b)]
  cross_21 <- r[cbind(2 * a, 2 * b - 1)]
  valid <- reliability_a > 0 & reliability_b > 0 & cross_12 * cross_21 >= 0
  rho <- rep(NA_real_, 15)

  # Positive trait loadings and independent errors across halves identify rho
  # from this product, even when paired loadings differ. Same-half correlations
  # can contain sample overlap and must not enter the numerator.
  rho[valid] <- sign(cross_12[valid] + cross_21[valid]) * sqrt(
    cross_12[valid] * cross_21[valid] / (reliability_a[valid] * reliability_b[valid])
  )
  # Do not clip estimates or project a matrix: that would change the estimator.
  tibble(
    quantity = pair_names, estimate = rho,
    status = case_when(
      reliability_a <= 0 | reliability_b <= 0 ~ "nonpositive_reliability",
      cross_12 * cross_21 < 0 ~ "opposite_cross_half_signs",
      rho == 0 ~ "zero_cross_half_correlation",
      abs(rho) > 1 ~ "outside_correlation_bounds",
      TRUE ~ "ok"
    )
  )
}

group_correlations <- function(rho) {
  c(mean_big_five = mean(rho[!iq_pair]),
    mean_iq_big_five = mean(rho[iq_pair]),
    difference_iq_minus_big_five = mean(rho[iq_pair]) - mean(rho[!iq_pair]))
}

results <- imap(inputs, \(input, analysis) {
  stopifnot(nrow(input$jackknife) == 200, all(is.finite(input$jackknife)))
  point <- closed_form_correlations(input$correlation)
  deleted <- map_dfr(1:200, \(block) {
    r <- input$correlation
    r[input$pairs] <- r[input$pairs[, 2:1]] <- input$jackknife[block, ]
    pair_estimates <- closed_form_correlations(r)
    groups <- enframe(group_correlations(pair_estimates$estimate), "quantity", "estimate") |>
      mutate(status = if_else(is.finite(estimate), "ok", "undefined_group"))
    bind_rows(pair_estimates, groups) |> mutate(block = block, .before = 1)
  })
  point <- bind_rows(
    point,
    enframe(group_correlations(point$estimate), "quantity", "estimate") |>
      mutate(status = if_else(is.finite(estimate), "ok", "undefined_group"))
  )
  summary <- deleted |>
    group_by(quantity) |>
    summarise(
      valid_replicates = sum(is.finite(estimate)),
      outside_bounds = sum(status == "outside_correlation_bounds"),
      zero_replicates = sum(status == "zero_cross_half_correlation"),
      # A single undefined deletion invalidates this quantity's jackknife SE.
      standard_error = if (all(is.finite(estimate))) {
        sqrt((199 / 200) * sum((estimate - mean(estimate))^2))
      } else NA_real_,
      .groups = "drop"
    ) |>
    right_join(point, by = "quantity") |>
    mutate(analysis = analysis, .before = 1,
           standard_error = if_else(is.finite(estimate), standard_error, NA_real_),
           se_status = if_else(is.finite(standard_error), "complete_200_blocks", "unavailable"))
  list(summary = summary, jackknife = mutate(deleted, analysis = analysis, .before = 1))
})

# ==============================================================================
# Report pairs, group means and uncertainty
# ==============================================================================

summary <- map(results, "summary") |> list_rbind()
correlations <- summary |> filter(quantity %in% pair_names)
groups <- summary |>
  filter(!quantity %in% pair_names) |>
  mutate(lower = estimate - qnorm(.975) * standard_error,
         upper = estimate + qnorm(.975) * standard_error,
         p_value = if_else(quantity == "difference_iq_minus_big_five" & standard_error > 0,
                           2 * pnorm(-abs(estimate / standard_error)), NA_real_))
write_csv(correlations, "output/closed_form_correlations/correlations.csv")
write_csv(groups, "output/closed_form_correlations/group_comparison.csv")
write_csv(map(results, "jackknife") |> list_rbind(), "output/closed_form_correlations/jackknife.csv")
write_csv(tibble(path = unname(paths), md5 = unname(tools::md5sum(paths))),
          "output/closed_form_correlations/sources.csv")
write_xlsx(list(correlations = correlations, group_comparison = groups),
           "output/closed_form_correlations/closed_form_correlations.xlsx")
writeLines(capture.output(sessionInfo()), "output/closed_form_correlations/session_info.txt")
print(correlations)
print(groups)
