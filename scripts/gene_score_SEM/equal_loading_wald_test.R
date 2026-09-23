#!/usr/bin/env Rscript

# Test an observable implication of equal paired trait loadings when
# sample/method factors are orthogonal to the trait factors and each other.
#
# For every pair of traits a and b, the null implies
#
#   cor(a_1, b_2) - cor(a_2, b_1) = 0.
#
# The joint Wald test uses the full delete-one-block jackknife sampling
# covariance matrix of the 66 off-diagonal correlations.

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) == 1L) {
  script_path <- normalizePath(
    sub("^--file=", "", script_argument),
    mustWork = TRUE
  )
} else {
  script_path <- normalizePath(
    "scripts/equal_loading_wald_test.R",
    mustWork = TRUE
  )
}

project_directory <- dirname(dirname(script_path))
arguments <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(arguments) >= 1L) {
  path.expand(arguments[[1L]])
} else {
  file.path(project_directory, "output", "wald_equal_loadings")
}
if (!grepl("^/", output_directory)) {
  output_directory <- file.path(project_directory, output_directory)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

indicator_map <- data.frame(
  trait = c(
    "Agreeableness",
    "Conscientiousness",
    "Extraversion",
    "Neuroticism",
    "Openness",
    "IQ"
  ),
  indicator_1 = c(
    "ReGPC_agr_half_one_no23_dir",
    "ReGPC_con_half_one_no23_dir",
    "ReGPC_ext_half_one_no23_dir",
    "ReGPC_neu_half_one_no23_dir",
    "ReGPC_ope_half_one_no23_dir",
    "iq_female_dir"
  ),
  indicator_2 = c(
    "ReGPC_agr_half_two_no23_dir",
    "ReGPC_con_half_two_no23_dir",
    "ReGPC_ext_half_two_no23_dir",
    "ReGPC_neu_half_two_no23_dir",
    "ReGPC_ope_half_two_no23_dir",
    "iq_male_dir"
  ),
  indicator_1_label = c(
    rep("Sample 1", 5L),
    "Female"
  ),
  indicator_2_label = c(
    rep("Sample 2", 5L),
    "Male"
  ),
  stringsAsFactors = FALSE
)

indicator_order <- as.vector(rbind(
  indicator_map$indicator_1,
  indicator_map$indicator_2
))
pair_indices <- t(combn(seq_along(indicator_order), 2L))
pair_names <- paste(
  indicator_order[pair_indices[, 1L]],
  indicator_order[pair_indices[, 2L]],
  sep = "__"
)

input_specifications <- data.frame(
  analysis = c("Gene score", "Gene-set enrichment"),
  analysis_id = c("gene_score", "gene_set_enrichment"),
  correlation_path = file.path(
    project_directory,
    "output",
    "correlation_matrices",
    c(
      "gene_score_correlation_matrix.csv",
      "all_gene_sets_correlation_matrix.csv"
    )
  ),
  sampling_vcov_path = file.path(
    project_directory,
    "output",
    "sampling_covariances",
    c(
      "gene_score_correlation_sampling_vcov.csv",
      "all_gene_sets_correlation_sampling_vcov.csv"
    )
  ),
  stringsAsFactors = FALSE
)

validate_symmetric_matrix <- function(
  matrix_value,
  description,
  tolerance = 1e-10
) {
  if (!is.matrix(matrix_value) ||
      nrow(matrix_value) != ncol(matrix_value)) {
    stop(description, " is not a square matrix.", call. = FALSE)
  }
  if (!all(is.finite(matrix_value))) {
    stop(description, " contains non-finite values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(
    matrix_value,
    t(matrix_value),
    tolerance = tolerance
  ))) {
    stop(description, " is not symmetric.", call. = FALSE)
  }
  invisible(TRUE)
}

read_named_matrix <- function(path, description) {
  if (!file.exists(path)) {
    stop("Missing ", description, ": ", path, call. = FALSE)
  }
  matrix_value <- as.matrix(
    read.csv(
      path,
      row.names = 1L,
      check.names = FALSE
    )
  )
  storage.mode(matrix_value) <- "double"
  validate_symmetric_matrix(matrix_value, description)
  matrix_value
}

resolve_pair_name <- function(left, right, available_names) {
  forward <- paste(left, right, sep = "__")
  reverse <- paste(right, left, sep = "__")
  if (forward %in% available_names) {
    return(forward)
  }
  if (reverse %in% available_names) {
    return(reverse)
  }
  NA_character_
}

canonical_pair_name <- function(left, right) {
  positions <- match(c(left, right), indicator_order)
  if (anyNA(positions) || positions[[1L]] == positions[[2L]]) {
    stop(
      "Could not construct a valid pair for ",
      left,
      " and ",
      right,
      ".",
      call. = FALSE
    )
  }
  ordered_names <- c(left, right)[order(positions)]
  paste(ordered_names, collapse = "__")
}

write_named_matrix <- function(matrix_value, path) {
  write.csv(
    matrix_value,
    path,
    row.names = TRUE,
    quote = TRUE
  )
}

run_wald_test <- function(input_row) {
  analysis <- input_row$analysis[[1L]]
  analysis_id <- input_row$analysis_id[[1L]]

  correlation_matrix <- read_named_matrix(
    input_row$correlation_path[[1L]],
    paste0(analysis, " correlation matrix")
  )
  if (!setequal(rownames(correlation_matrix), indicator_order) ||
      !setequal(colnames(correlation_matrix), indicator_order)) {
    stop(
      analysis,
      " correlation matrix does not contain the expected 12 indicators.",
      call. = FALSE
    )
  }
  correlation_matrix <- correlation_matrix[
    indicator_order,
    indicator_order,
    drop = FALSE
  ]
  if (max(abs(diag(correlation_matrix) - 1)) > 1e-10) {
    stop(
      analysis,
      " correlation matrix does not have a unit diagonal.",
      call. = FALSE
    )
  }

  sampling_vcov <- read_named_matrix(
    input_row$sampling_vcov_path[[1L]],
    paste0(analysis, " correlation sampling covariance matrix")
  )
  source_pair_names <- vapply(
    seq_len(nrow(pair_indices)),
    function(index) {
      resolve_pair_name(
        indicator_order[pair_indices[index, 1L]],
        indicator_order[pair_indices[index, 2L]],
        rownames(sampling_vcov)
      )
    },
    character(1)
  )
  if (anyNA(source_pair_names) ||
      !all(source_pair_names %in% colnames(sampling_vcov))) {
    stop(
      analysis,
      " sampling covariance does not contain all 66 expected correlations.",
      call. = FALSE
    )
  }
  sampling_vcov <- sampling_vcov[
    source_pair_names,
    source_pair_names,
    drop = FALSE
  ]
  dimnames(sampling_vcov) <- list(pair_names, pair_names)

  observed_correlations <- correlation_matrix[pair_indices]
  names(observed_correlations) <- pair_names

  trait_pairs <- t(combn(seq_len(nrow(indicator_map)), 2L))
  contrast_matrix <- matrix(
    0,
    nrow = nrow(trait_pairs),
    ncol = length(pair_names),
    dimnames = list(
      rep("", nrow(trait_pairs)),
      pair_names
    )
  )
  contrast_rows <- vector("list", nrow(trait_pairs))

  for (index in seq_len(nrow(trait_pairs))) {
    trait_a <- indicator_map[trait_pairs[index, 1L], , drop = FALSE]
    trait_b <- indicator_map[trait_pairs[index, 2L], , drop = FALSE]

    a1_b2_name <- canonical_pair_name(
      trait_a$indicator_1[[1L]],
      trait_b$indicator_2[[1L]]
    )
    a2_b1_name <- canonical_pair_name(
      trait_a$indicator_2[[1L]],
      trait_b$indicator_1[[1L]]
    )
    contrast_matrix[index, a1_b2_name] <- 1
    contrast_matrix[index, a2_b1_name] <- -1

    contrast_id <- paste(
      trait_a$trait[[1L]],
      trait_b$trait[[1L]],
      sep = "__"
    )
    rownames(contrast_matrix)[[index]] <- contrast_id
    contrast_rows[[index]] <- data.frame(
      analysis = analysis,
      analysis_id = analysis_id,
      contrast_id = contrast_id,
      trait_a = trait_a$trait[[1L]],
      trait_b = trait_b$trait[[1L]],
      correlation_a1_b2_name = a1_b2_name,
      correlation_a2_b1_name = a2_b1_name,
      correlation_a1_b2 = observed_correlations[[a1_b2_name]],
      correlation_a2_b1 = observed_correlations[[a2_b1_name]],
      stringsAsFactors = FALSE
    )
  }

  contrast_table <- do.call(rbind, contrast_rows)
  contrast_estimates <- drop(contrast_matrix %*% observed_correlations)
  contrast_vcov <- contrast_matrix %*%
    sampling_vcov %*%
    t(contrast_matrix)
  contrast_vcov <- (contrast_vcov + t(contrast_vcov)) / 2
  contrast_eigenvalues <- eigen(
    contrast_vcov,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  contrast_rank <- qr(contrast_vcov, tol = 1e-10)$rank
  if (contrast_rank != nrow(contrast_matrix) ||
      min(contrast_eigenvalues) <= 0) {
    stop(
      analysis,
      " contrast covariance is not positive definite and full rank.",
      call. = FALSE
    )
  }

  contrast_precision <- chol2inv(chol(contrast_vcov))
  wald_statistic <- drop(
    crossprod(
      contrast_estimates,
      contrast_precision %*% contrast_estimates
    )
  )
  degrees_of_freedom <- contrast_rank
  omnibus_p_value <- pchisq(
    wald_statistic,
    df = degrees_of_freedom,
    lower.tail = FALSE
  )

  contrast_standard_errors <- sqrt(diag(contrast_vcov))
  contrast_z <- contrast_estimates / contrast_standard_errors
  contrast_p_values <- 2 * pnorm(abs(contrast_z), lower.tail = FALSE)
  normal_critical_value <- qnorm(0.975)

  contrast_table$difference <- contrast_estimates
  contrast_table$standard_error <- contrast_standard_errors
  contrast_table$ci_95_lower <- contrast_estimates -
    normal_critical_value * contrast_standard_errors
  contrast_table$ci_95_upper <- contrast_estimates +
    normal_critical_value * contrast_standard_errors
  contrast_table$z_value <- contrast_z
  contrast_table$p_value <- contrast_p_values
  contrast_table$p_value_holm <- p.adjust(contrast_p_values, method = "holm")
  contrast_table$significant_unadjusted_05 <-
    contrast_table$p_value < 0.05
  contrast_table$significant_holm_05 <-
    contrast_table$p_value_holm < 0.05

  omnibus_table <- data.frame(
    analysis = analysis,
    analysis_id = analysis_id,
    null_hypothesis = paste0(
      "cor(a1,b2) = cor(a2,b1) for all 15 trait pairs"
    ),
    wald_statistic = wald_statistic,
    degrees_of_freedom = degrees_of_freedom,
    p_value = omnibus_p_value,
    reject_at_05 = omnibus_p_value < 0.05,
    number_of_contrasts = nrow(contrast_matrix),
    contrast_covariance_rank = contrast_rank,
    contrast_covariance_min_eigenvalue = min(contrast_eigenvalues),
    contrast_covariance_condition_number =
      max(contrast_eigenvalues) / min(contrast_eigenvalues),
    correlation_source = normalizePath(
      input_row$correlation_path[[1L]],
      mustWork = TRUE
    ),
    sampling_vcov_source = normalizePath(
      input_row$sampling_vcov_path[[1L]],
      mustWork = TRUE
    ),
    stringsAsFactors = FALSE
  )

  subset_definitions <- list(
    personality_only = (
      contrast_table$trait_a != "IQ" &
        contrast_table$trait_b != "IQ"
    ),
    iq_involving = (
      contrast_table$trait_a == "IQ" |
        contrast_table$trait_b == "IQ"
    )
  )
  subset_labels <- c(
    personality_only = "Personality-trait pairs only",
    iq_involving = "Pairs involving IQ"
  )
  subset_tables <- lapply(
    names(subset_definitions),
    function(subset_id) {
      keep <- subset_definitions[[subset_id]]
      subset_estimates <- contrast_estimates[keep]
      subset_vcov <- contrast_vcov[keep, keep, drop = FALSE]
      subset_eigenvalues <- eigen(
        subset_vcov,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      subset_rank <- qr(subset_vcov, tol = 1e-10)$rank
      if (subset_rank != sum(keep) || min(subset_eigenvalues) <= 0) {
        stop(
          analysis,
          " ",
          subset_id,
          " contrast covariance is not positive definite and full rank.",
          call. = FALSE
        )
      }
      subset_precision <- chol2inv(chol(subset_vcov))
      subset_wald <- drop(
        crossprod(
          subset_estimates,
          subset_precision %*% subset_estimates
        )
      )
      subset_p_value <- pchisq(
        subset_wald,
        df = subset_rank,
        lower.tail = FALSE
      )
      data.frame(
        analysis = analysis,
        analysis_id = analysis_id,
        subset = subset_labels[[subset_id]],
        subset_id = subset_id,
        wald_statistic = subset_wald,
        degrees_of_freedom = subset_rank,
        p_value = subset_p_value,
        reject_at_05 = subset_p_value < 0.05,
        number_of_contrasts = sum(keep),
        contrast_covariance_rank = subset_rank,
        contrast_covariance_min_eigenvalue = min(subset_eigenvalues),
        contrast_covariance_condition_number =
          max(subset_eigenvalues) / min(subset_eigenvalues),
        stringsAsFactors = FALSE
      )
    }
  )
  subset_table <- do.call(rbind, subset_tables)

  write_named_matrix(
    contrast_vcov,
    file.path(
      output_directory,
      paste0(analysis_id, "_equal_loading_contrast_vcov.csv")
    )
  )
  write_named_matrix(
    contrast_matrix,
    file.path(
      output_directory,
      paste0(analysis_id, "_equal_loading_contrast_matrix.csv")
    )
  )

  list(
    omnibus = omnibus_table,
    subsets = subset_table,
    contrasts = contrast_table
  )
}

test_results <- lapply(
  seq_len(nrow(input_specifications)),
  function(index) {
    run_wald_test(input_specifications[index, , drop = FALSE])
  }
)

omnibus_results <- do.call(
  rbind,
  lapply(test_results, function(result) result$omnibus)
)
contrast_results <- do.call(
  rbind,
  lapply(test_results, function(result) result$contrasts)
)
subset_results <- do.call(
  rbind,
  lapply(test_results, function(result) result$subsets)
)

omnibus_path <- file.path(
  output_directory,
  "equal_loading_wald_omnibus.csv"
)
contrasts_path <- file.path(
  output_directory,
  "equal_loading_wald_contrasts.csv"
)
subsets_path <- file.path(
  output_directory,
  "equal_loading_wald_subsets.csv"
)
report_path <- file.path(
  output_directory,
  "equal_loading_wald_report.md"
)

write.csv(
  omnibus_results,
  omnibus_path,
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  contrast_results,
  contrasts_path,
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  subset_results,
  subsets_path,
  row.names = FALSE,
  quote = TRUE
)

format_probability <- function(value) {
  if (value < 0.001) {
    return(format(value, digits = 3L, scientific = TRUE))
  }
  sprintf("%.4f", value)
}

plot_contrasts <- function(path, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(
      filename = path,
      width = 3000,
      height = 1800,
      res = 300
    )
  } else {
    pdf(
      file = path,
      width = 12,
      height = 7.2,
      useDingbats = FALSE
    )
  }
  on.exit(dev.off(), add = TRUE)

  global_limit <- max(abs(c(
    contrast_results$ci_95_lower,
    contrast_results$ci_95_upper
  )))
  global_limit <- ceiling(global_limit * 100) / 100
  point_colors <- c(
    nonsignificant = "#425466",
    significant_raw = "#B33A3A"
  )

  old_parameters <- par(
    mfrow = c(1L, 2L),
    mar = c(4.6, 10.2, 3.2, 1.2),
    oma = c(2.4, 0.2, 1.5, 0.2),
    las = 1,
    xaxs = "i"
  )
  on.exit(par(old_parameters), add = TRUE)

  for (analysis_name in omnibus_results$analysis) {
    analysis_contrasts <- contrast_results[
      contrast_results$analysis == analysis_name,
      ,
      drop = FALSE
    ]
    omnibus_row <- omnibus_results[
      omnibus_results$analysis == analysis_name,
      ,
      drop = FALSE
    ]
    y_positions <- rev(seq_len(nrow(analysis_contrasts)))
    pair_labels <- paste(
      analysis_contrasts$trait_a,
      analysis_contrasts$trait_b,
      sep = " - "
    )
    colors <- ifelse(
      analysis_contrasts$significant_unadjusted_05,
      point_colors[["significant_raw"]],
      point_colors[["nonsignificant"]]
    )

    plot(
      NA_real_,
      NA_real_,
      xlim = c(-global_limit, global_limit),
      ylim = c(0.5, nrow(analysis_contrasts) + 0.5),
      xlab = expression(r(a[1], b[2]) - r(a[2], b[1])),
      ylab = "",
      yaxt = "n",
      bty = "n",
      main = paste0(
        analysis_name,
        "\nW(",
        omnibus_row$degrees_of_freedom[[1L]],
        ") = ",
        sprintf("%.2f", omnibus_row$wald_statistic[[1L]]),
        ", p = ",
        format_probability(omnibus_row$p_value[[1L]])
      )
    )
    abline(v = 0, col = "#9AA6B2", lty = 2L, lwd = 1.2)
    abline(
      h = seq_len(nrow(analysis_contrasts)),
      col = "#EEF1F4",
      lwd = 0.8
    )
    segments(
      x0 = analysis_contrasts$ci_95_lower,
      y0 = y_positions,
      x1 = analysis_contrasts$ci_95_upper,
      y1 = y_positions,
      col = colors,
      lwd = 2
    )
    points(
      analysis_contrasts$difference,
      y_positions,
      pch = 21L,
      bg = colors,
      col = "white",
      cex = 1.05
    )
    axis(
      side = 2L,
      at = y_positions,
      labels = pair_labels,
      tick = FALSE,
      cex.axis = 0.68
    )
  }

  mtext(
    paste0(
      "95% jackknife Wald intervals. Red points have raw p < .05; ",
      "none is Holm-significant. The joint omnibus tests are primary."
    ),
    side = 1L,
    outer = TRUE,
    line = 0.6,
    cex = 0.78,
    col = "#425466"
  )
  mtext(
    "Cross-sample symmetry contrasts",
    side = 3L,
    outer = TRUE,
    line = 0.2,
    cex = 1.08,
    font = 2L
  )
}

plot_png_path <- file.path(
  output_directory,
  "equal_loading_wald_contrasts.png"
)
plot_pdf_path <- file.path(
  output_directory,
  "equal_loading_wald_contrasts.pdf"
)
plot_contrasts(plot_png_path, device = "png")
plot_contrasts(plot_pdf_path, device = "pdf")

report_lines <- c(
  "# Equal paired-loading Wald tests",
  "",
  paste0(
    "This analysis tests the 15 observable symmetry restrictions ",
    "`cor(a1,b2) = cor(a2,b1)` implied jointly by equal paired trait ",
    "loadings and orthogonal sample/method factors."
  ),
  "",
  paste0(
    "This is a necessary but not sufficient test of exact loading equality. ",
    "Under orthogonality, the contrasts test whether paired-loading ratios ",
    "`lambda[a1]/lambda[a2]` are the same across traits. A common non-unit ",
    "ratio for every trait would satisfy all 15 restrictions."
  ),
  "",
  "The Wald statistic uses the full delete-one-block jackknife covariance:",
  "",
  "`W = d' [D Var_JK(r) D']^{-1} d`, where `d = D r`.",
  "",
  "## Omnibus tests",
  ""
)

for (index in seq_len(nrow(omnibus_results))) {
  omnibus_row <- omnibus_results[index, , drop = FALSE]
  interpretation <- if (omnibus_row$reject_at_05[[1L]]) {
    "Reject the joint symmetry restrictions at alpha = .05."
  } else {
    "Do not reject the joint symmetry restrictions at alpha = .05."
  }
  report_lines <- c(
    report_lines,
    paste0(
      "- ",
      omnibus_row$analysis[[1L]],
      ": W(",
      omnibus_row$degrees_of_freedom[[1L]],
      ") = ",
      sprintf("%.4f", omnibus_row$wald_statistic[[1L]]),
      ", p = ",
      format_probability(omnibus_row$p_value[[1L]]),
      ". ",
      interpretation
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "## Joint tests localized by contrast type",
  "",
  paste0(
    "These are secondary joint Wald tests of the ten personality-only ",
    "contrasts and the five contrasts involving IQ. They are not an ",
    "additive partition of the omnibus statistic because the contrast ",
    "groups are sampling-correlated."
  ),
  "",
  "| Analysis | Contrast subset | W | df | p | Reject at .05 |",
  "|---|---|---:|---:|---:|:---:|"
)

for (index in seq_len(nrow(subset_results))) {
  subset_row <- subset_results[index, , drop = FALSE]
  report_lines <- c(
    report_lines,
    paste0(
      "| ",
      subset_row$analysis[[1L]],
      " | ",
      subset_row$subset[[1L]],
      " | ",
      sprintf("%.4f", subset_row$wald_statistic[[1L]]),
      " | ",
      subset_row$degrees_of_freedom[[1L]],
      " | ",
      format_probability(subset_row$p_value[[1L]]),
      " | ",
      if (subset_row$reject_at_05[[1L]]) "Yes" else "No",
      " |"
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "## Largest individual asymmetries",
  "",
  paste0(
    "Differences are `cor(a1,b2) - cor(a2,b1)`. Rows are ordered by the ",
    "absolute z statistic. These individual tests are secondary to the ",
    "omnibus test."
  )
)

for (analysis_name in omnibus_results$analysis) {
  analysis_contrasts <- contrast_results[
    contrast_results$analysis == analysis_name,
    ,
    drop = FALSE
  ]
  analysis_contrasts <- analysis_contrasts[
    order(abs(analysis_contrasts$z_value), decreasing = TRUE),
    ,
    drop = FALSE
  ]
  displayed_contrasts <- head(analysis_contrasts, 5L)
  report_lines <- c(
    report_lines,
    "",
    paste0("### ", analysis_name),
    "",
    "| Trait pair | Difference | SE | z | Raw p | Holm p |",
    "|---|---:|---:|---:|---:|---:|"
  )
  for (index in seq_len(nrow(displayed_contrasts))) {
    contrast_row <- displayed_contrasts[index, , drop = FALSE]
    report_lines <- c(
      report_lines,
      paste0(
        "| ",
        contrast_row$trait_a[[1L]],
        "–",
        contrast_row$trait_b[[1L]],
        " | ",
        sprintf("%.4f", contrast_row$difference[[1L]]),
        " | ",
        sprintf("%.4f", contrast_row$standard_error[[1L]]),
        " | ",
        sprintf("%.3f", contrast_row$z_value[[1L]]),
        " | ",
        format_probability(contrast_row$p_value[[1L]]),
        " | ",
        format_probability(contrast_row$p_value_holm[[1L]]),
        " |"
      )
    )
  }
  report_lines <- c(
    report_lines,
    "",
    paste0(
      "Unadjusted p < .05: ",
      sum(analysis_contrasts$significant_unadjusted_05),
      " of 15. Holm-adjusted p < .05: ",
      sum(analysis_contrasts$significant_holm_05),
      " of 15."
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "## Interpretation boundary",
  "",
  paste0(
    "Rejection shows that the observed correlation matrix violates the ",
    "symmetry restrictions. With orthogonality maintained, it shows that ",
    "the paired-loading ratios are not uniform across traits; it does not ",
    "identify either member of a loading pair as the source. Without ",
    "maintaining orthogonality, it also cannot distinguish loading-ratio ",
    "differences from nonorthogonal method effects or another violated ",
    "assumption. Individual contrast p-values are reported both raw and ",
    "Holm-adjusted; the omnibus test is primary."
  ),
  "",
  "## Output files",
  "",
  "- `equal_loading_wald_omnibus.csv`: joint tests and numerical diagnostics.",
  "- `equal_loading_wald_subsets.csv`: personality-only and IQ-involving joint tests.",
  "- `equal_loading_wald_contrasts.csv`: all 15 contrasts per analysis.",
  "- `*_equal_loading_contrast_vcov.csv`: propagated contrast covariance.",
  "- `*_equal_loading_contrast_matrix.csv`: the exact contrast definitions.",
  "- `equal_loading_wald_contrasts.png` and `.pdf`: confidence-interval plot."
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Wrote omnibus tests: ", omnibus_path)
message("Wrote subset tests: ", subsets_path)
message("Wrote contrast details: ", contrasts_path)
message("Wrote contrast plots: ", plot_png_path, " and ", plot_pdf_path)
message("Wrote report: ", report_path)

print(
  omnibus_results[
    ,
    c(
      "analysis",
      "wald_statistic",
      "degrees_of_freedom",
      "p_value",
      "reject_at_05"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)
