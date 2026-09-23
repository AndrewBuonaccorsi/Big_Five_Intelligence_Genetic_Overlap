#!/usr/bin/env Rscript

# Sensitivity analysis for proportional split-half trait loadings.
#
# This script fits the repository's two SEM structures to the pooled gene-score
# and gene-set correlation matrices under three paired-loading specifications:
#
#   1. equal:  lambda_2 / lambda_1 = 1 for every trait (reference model),
#   2. common: one estimated ratio for the five personality traits, with the
#              female and male IQ loadings constrained equal, and
#   3. fixed_chisq: trait-specific ratios fixed from the supplied GWAS
#                   summary-statistic mean chi-square,
#   4. fixed_magma_chisq: trait-specific ratios fixed from mean MAGMA ZSTAT^2.
#
# The fixed ratio is defined on the correlation scale. If a standardized score
# has reliability R_h = (mean_chisq_h - 1) / mean_chisq_h and its loading is
# proportional to sqrt(R_h), then
#
#   kappa = lambda_2 / lambda_1 = sqrt(R_2 / R_1).
#
# The MAGMA calibration uses the actual gene-level values that form the
# gene-score indicators, rather than the upstream SNP-level GWAS statistics.
# This is deliberately standalone: it does not change or overwrite the
# canonical equal-loading SEM outputs.

suppressPackageStartupMessages({
  library(Matrix)
  library(numDeriv)
})

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
    "scripts/split_half_ratio_models.R",
    mustWork = TRUE
  )
}

project_directory <- dirname(dirname(script_path))
arguments <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(arguments) >= 1L) {
  path.expand(arguments[[1L]])
} else {
  file.path(project_directory, "output", "split_half_ratio_models")
}
if (!grepl("^/", output_directory)) {
  output_directory <- file.path(project_directory, output_directory)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

results_path <- file.path(
  output_directory,
  "split_half_ratio_model_results.csv"
)
comparisons_path <- file.path(
  output_directory,
  "split_half_ratio_model_comparisons.csv"
)
calibrations_path <- file.path(
  output_directory,
  "split_half_ratio_model_calibrations.csv"
)
mean_correlations_path <- file.path(
  output_directory,
  "split_half_ratio_model_mean_correlations.csv"
)
report_path <- file.path(
  output_directory,
  "split_half_ratio_model_report.md"
)

jackknife_blocks <- 200L
long_names <- c(
  "ReGPC_agr_half_one_no23_dir",
  "ReGPC_agr_half_two_no23_dir",
  "ReGPC_con_half_one_no23_dir",
  "ReGPC_con_half_two_no23_dir",
  "ReGPC_ext_half_one_no23_dir",
  "ReGPC_ext_half_two_no23_dir",
  "ReGPC_neu_half_one_no23_dir",
  "ReGPC_neu_half_two_no23_dir",
  "ReGPC_ope_half_one_no23_dir",
  "ReGPC_ope_half_two_no23_dir",
  "iq_female_dir",
  "iq_male_dir"
)
short_names <- c(
  "agr1", "agr2",
  "con1", "con2",
  "ext1", "ext2",
  "neu1", "neu2",
  "ope1", "ope2",
  "iq_female", "iq_male"
)
trait_names <- c(
  "Agreeableness",
  "Conscientiousness",
  "Extraversion",
  "Neuroticism",
  "Openness",
  "IQ"
)
trait_ids <- c("agree", "consc", "extra", "neurot", "open", "iq")
indicator_trait <- rep(trait_ids, each = 2L)
indicator_trait_index <- match(indicator_trait, trait_ids)
indicator_half <- rep(1:2, length(trait_ids))
personality_indicator_indices <- seq_len(10L)
indicator_sample <- rep(1:2, 5L)

indicator_map <- data.frame(
  trait = trait_names,
  trait_id = trait_ids,
  indicator_1 = long_names[seq(1L, 12L, by = 2L)],
  indicator_2 = long_names[seq(2L, 12L, by = 2L)],
  indicator_1_label = c(rep("Sample 1", 5L), "Female"),
  indicator_2_label = c(rep("Sample 2", 5L), "Male"),
  stringsAsFactors = FALSE
)

mean_chisq_1 <- c(
  Agreeableness = 1.120713,
  Conscientiousness = 1.172888,
  Extraversion = 1.201807,
  Neuroticism = 1.169212,
  Openness = 1.178110,
  IQ = 1.275822
)
mean_chisq_2 <- c(
  Agreeableness = 1.181906,
  Conscientiousness = 1.323071,
  Extraversion = 1.343517,
  Neuroticism = 1.340171,
  Openness = 1.362600,
  IQ = 1.211333
)
chisq_signal_1 <- (mean_chisq_1 - 1) / mean_chisq_1
chisq_signal_2 <- (mean_chisq_2 - 1) / mean_chisq_2
fixed_ratio <- sqrt(chisq_signal_2 / chisq_signal_1)
raw_excess_chisq_ratio <- sqrt(
  (mean_chisq_2 - 1) / (mean_chisq_1 - 1)
)

# The SEM gene-score indicators are MAGMA gene-level Z statistics.  The
# MAGMA .genes.raw column 9 is ZSTAT; squaring it gives a one-df chi-square
# analogue under the null.  For each paired indicator, use the same common
# gene intersection in both halves so that the two means are comparable.
magma_input_directory <- file.path(
  project_directory,
  "data",
  "gene_enrichments"
)
magma_score_paths <- file.path(
  magma_input_directory,
  paste0(long_names, ".genes.raw")
)
names(magma_score_paths) <- long_names

read_magma_z_scores <- function(path) {
  if (!file.exists(path)) {
    stop("Missing MAGMA gene-level input: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#", lines)]
  fields <- strsplit(trimws(lines), "[[:space:]]+")
  if (any(lengths(fields) < 9L)) {
    stop("MAGMA input has fewer than nine fields: ", path, call. = FALSE)
  }
  gene_ids <- vapply(fields, `[[`, character(1), 1L)
  z_scores <- as.numeric(vapply(fields, `[[`, character(1), 9L))
  if (anyNA(gene_ids) || anyNA(z_scores) || !all(is.finite(z_scores))) {
    stop("Invalid MAGMA gene-level Z statistics: ", path, call. = FALSE)
  }
  if (anyDuplicated(gene_ids)) {
    stop("Duplicated genes in MAGMA input: ", path, call. = FALSE)
  }
  names(z_scores) <- gene_ids
  z_scores
}

magma_z_scores <- lapply(magma_score_paths, read_magma_z_scores)

magma_calibration_rows <- lapply(seq_along(trait_names), function(index) {
  indicator_1 <- indicator_map$indicator_1[[index]]
  indicator_2 <- indicator_map$indicator_2[[index]]
  common_genes <- intersect(
    names(magma_z_scores[[indicator_1]]),
    names(magma_z_scores[[indicator_2]])
  )
  if (length(common_genes) < 2L) {
    stop(
      "Fewer than two common genes for MAGMA calibration of ",
      trait_names[[index]],
      call. = FALSE
    )
  }
  z_1 <- magma_z_scores[[indicator_1]][common_genes]
  z_2 <- magma_z_scores[[indicator_2]][common_genes]
  mean_1 <- mean(z_1^2)
  mean_2 <- mean(z_2^2)
  if (mean_1 <= 1 || mean_2 <= 1) {
    stop(
      "MAGMA mean Z-squared must exceed the null value of one for ",
      trait_names[[index]],
      ". Observed values: ",
      format(mean_1, digits = 6),
      " and ",
      format(mean_2, digits = 6),
      call. = FALSE
    )
  }
  data.frame(
    trait = trait_names[[index]],
    indicator_1 = indicator_1,
    indicator_2 = indicator_2,
    common_gene_count = length(common_genes),
    magma_mean_chisq_1 = mean_1,
    magma_mean_chisq_2 = mean_2,
    magma_chisq_signal_1 = (mean_1 - 1) / mean_1,
    magma_chisq_signal_2 = (mean_2 - 1) / mean_2,
    magma_kappa_raw_excess_chisq = sqrt(
      (mean_2 - 1) / (mean_1 - 1)
    ),
    magma_kappa_correlation_scale = sqrt(
      ((mean_2 - 1) / mean_2) /
        ((mean_1 - 1) / mean_1)
    ),
    stringsAsFactors = FALSE
  )
})
magma_calibration <- do.call(rbind, magma_calibration_rows)
rownames(magma_calibration) <- magma_calibration$trait
magma_mean_chisq_1 <- setNames(
  magma_calibration$magma_mean_chisq_1,
  magma_calibration$trait
)
magma_mean_chisq_2 <- setNames(
  magma_calibration$magma_mean_chisq_2,
  magma_calibration$trait
)
magma_chisq_signal_1 <- setNames(
  magma_calibration$magma_chisq_signal_1,
  magma_calibration$trait
)
magma_chisq_signal_2 <- setNames(
  magma_calibration$magma_chisq_signal_2,
  magma_calibration$trait
)
magma_raw_excess_chisq_ratio <- setNames(
  magma_calibration$magma_kappa_raw_excess_chisq,
  magma_calibration$trait
)
magma_fixed_ratio <- setNames(
  magma_calibration$magma_kappa_correlation_scale,
  magma_calibration$trait
)
write.csv(
  magma_calibration,
  calibrations_path,
  row.names = FALSE,
  quote = TRUE,
  na = ""
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

pair_indices <- t(combn(seq_along(short_names), 2L))
pair_names <- paste(
  short_names[pair_indices[, 1L]],
  short_names[pair_indices[, 2L]],
  sep = "__"
)
trait_pair_indices <- t(combn(seq_along(trait_ids), 2L))

validate_symmetric_matrix <- function(x, description, tolerance = 1e-10) {
  if (!is.matrix(x) || nrow(x) != ncol(x)) {
    stop(description, " is not square.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop(description, " contains non-finite values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = tolerance))) {
    stop(description, " is not symmetric.", call. = FALSE)
  }
  invisible(TRUE)
}

read_named_matrix <- function(path, description) {
  if (!file.exists(path)) {
    stop("Missing ", description, ": ", path, call. = FALSE)
  }
  value <- as.matrix(read.csv(path, row.names = 1L, check.names = FALSE))
  storage.mode(value) <- "double"
  validate_symmetric_matrix(value, description)
  value
}

resolve_pair_name <- function(left, right, available_names) {
  forward <- paste(left, right, sep = "__")
  reverse <- paste(right, left, sep = "__")
  if (forward %in% available_names) {
    forward
  } else if (reverse %in% available_names) {
    reverse
  } else {
    NA_character_
  }
}

correlation_matrix_from_parameters <- function(parameters, dimension) {
  lower <- diag(dimension)
  lower[lower.tri(lower)] <- parameters
  covariance <- tcrossprod(lower)
  covariance / sqrt(outer(diag(covariance), diag(covariance)))
}

load_analysis_inputs <- function(input_row) {
  analysis <- input_row$analysis[[1L]]
  correlation_matrix <- read_named_matrix(
    input_row$correlation_path[[1L]],
    paste0(analysis, " correlation matrix")
  )
  if (!setequal(rownames(correlation_matrix), long_names) ||
      !setequal(colnames(correlation_matrix), long_names)) {
    stop(
      analysis,
      " correlation matrix does not contain the expected indicators.",
      call. = FALSE
    )
  }
  correlation_matrix <- correlation_matrix[
    long_names,
    long_names,
    drop = FALSE
  ]
  dimnames(correlation_matrix) <- list(short_names, short_names)
  if (max(abs(diag(correlation_matrix) - 1)) > 1e-10) {
    stop(analysis, " correlation matrix does not have unit diagonal.")
  }

  sampling_vcov <- read_named_matrix(
    input_row$sampling_vcov_path[[1L]],
    paste0(analysis, " sampling covariance matrix")
  )
  source_pair_names <- vapply(
    seq_len(nrow(pair_indices)),
    function(index) {
      resolve_pair_name(
        long_names[pair_indices[index, 1L]],
        long_names[pair_indices[index, 2L]],
        rownames(sampling_vcov)
      )
    },
    character(1)
  )
  if (anyNA(source_pair_names) ||
      !all(source_pair_names %in% colnames(sampling_vcov))) {
    stop(
      analysis,
      " sampling covariance lacks one or more expected correlations.",
      call. = FALSE
    )
  }
  sampling_vcov <- sampling_vcov[
    source_pair_names,
    source_pair_names,
    drop = FALSE
  ]
  dimnames(sampling_vcov) <- list(pair_names, pair_names)

  correlation_eigenvalues <- eigen(
    correlation_matrix,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  sampling_eigenvalues <- eigen(
    sampling_vcov,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (min(correlation_eigenvalues) <= 0) {
    stop(analysis, " correlation matrix is not positive definite.")
  }
  if (min(sampling_eigenvalues) <= 0) {
    stop(analysis, " sampling covariance is not positive definite.")
  }

  list(
    correlation_matrix = correlation_matrix,
    sampling_vcov = sampling_vcov,
    sampling_precision = solve(sampling_vcov),
    objective_precision =
      ((jackknife_blocks - 1) / jackknife_blocks) * solve(sampling_vcov)
  )
}

parameter_layout <- function(structure, loading_model) {
  ratio_count <- if (loading_model == "common") 1L else 0L
  structure_count <- if (structure == "correlated") {
    choose(length(trait_ids), 2L)
  } else {
    length(trait_ids)
  }
  starts <- c(
    trait = 1L,
    method = 1L + length(trait_ids),
    structure = 1L + length(trait_ids) +
      length(personality_indicator_indices),
    ratio = 1L + length(trait_ids) +
      length(personality_indicator_indices) + structure_count
  )
  count <- length(trait_ids) +
    length(personality_indicator_indices) +
    structure_count + ratio_count
  list(
    trait = starts[["trait"]] + seq_len(length(trait_ids)) - 1L,
    method = starts[["method"]] +
      seq_len(length(personality_indicator_indices)) - 1L,
    structure = starts[["structure"]] + seq_len(structure_count) - 1L,
    ratio = if (ratio_count > 0L) {
      starts[["ratio"]] + seq_len(ratio_count) - 1L
    } else {
      integer()
    },
    count = count
  )
}

decode_parameters <- function(
  parameters,
  structure,
  loading_model,
  layout
) {
  if (loading_model == "equal") {
    ratios <- rep(1, length(trait_ids))
  } else if (loading_model == "fixed_chisq") {
    ratios <- unname(fixed_ratio[trait_names])
  } else if (loading_model == "fixed_magma_chisq") {
    ratios <- unname(magma_fixed_ratio[trait_names])
  } else {
    ratio_parameters <- parameters[layout$ratio]
    ratios <- c(
      rep(exp(ratio_parameters[[1L]]), 5L),
      1
    )
  }
  names(ratios) <- trait_ids

  # This cap guarantees that both loadings in every pair are below one in
  # absolute value while preserving lambda_2 = kappa * lambda_1 exactly.
  first_loading_cap <- 0.999 / pmax(1, ratios)
  first_loadings <- first_loading_cap * tanh(parameters[layout$trait])
  second_loadings <- ratios * first_loadings
  trait_loadings <- as.vector(rbind(first_loadings, second_loadings))
  names(first_loadings) <- names(second_loadings) <- trait_ids
  names(trait_loadings) <- short_names

  method_cap <- 0.999 * sqrt(pmax(
    1e-12,
    1 - trait_loadings[personality_indicator_indices]^2
  ))
  method_loadings <- method_cap * tanh(parameters[layout$method])
  names(method_loadings) <- short_names[personality_indicator_indices]

  if (structure == "correlated") {
    trait_correlation <- correlation_matrix_from_parameters(
      parameters[layout$structure],
      length(trait_ids)
    )
    general_loadings <- rep(NA_real_, length(trait_ids))
  } else {
    general_loadings <- 0.999 * tanh(parameters[layout$structure])
    trait_correlation <- tcrossprod(general_loadings)
    diag(trait_correlation) <- 1
  }
  dimnames(trait_correlation) <- list(trait_ids, trait_ids)
  names(general_loadings) <- trait_ids

  measurement_matrix <- matrix(
    0,
    nrow = length(short_names),
    ncol = length(trait_ids),
    dimnames = list(short_names, trait_ids)
  )
  measurement_matrix[
    cbind(seq_along(short_names), indicator_trait_index)
  ] <- trait_loadings

  method_matrix <- matrix(
    0,
    nrow = length(short_names),
    ncol = 2L,
    dimnames = list(short_names, c("sample_one", "sample_two"))
  )
  method_matrix[
    cbind(personality_indicator_indices, indicator_sample)
  ] <- method_loadings

  implied_correlation <-
    measurement_matrix %*% trait_correlation %*% t(measurement_matrix) +
    tcrossprod(method_matrix)
  residual_variances <- 1 - diag(implied_correlation)
  diag(implied_correlation) <- 1
  dimnames(implied_correlation) <- list(short_names, short_names)

  list(
    ratios = ratios,
    first_loadings = first_loadings,
    second_loadings = second_loadings,
    trait_loadings = trait_loadings,
    method_loadings = method_loadings,
    general_loadings = general_loadings,
    trait_correlation = trait_correlation,
    residual_variances = residual_variances,
    implied_correlation = implied_correlation
  )
}

make_start <- function(
  observed_matrix,
  structure,
  loading_model,
  layout,
  start_number
) {
  ratios <- if (loading_model == "fixed_chisq") {
    unname(fixed_ratio[trait_names])
  } else if (loading_model == "fixed_magma_chisq") {
    unname(magma_fixed_ratio[trait_names])
  } else if (loading_model == "common") {
    c(
      rep(exp(mean(log(fixed_ratio[trait_names[1:5]]))), 5L),
      1
    )
  } else {
    rep(1, length(trait_ids))
  }
  within_trait_correlations <- observed_matrix[
    cbind(seq(1L, 12L, by = 2L), seq(2L, 12L, by = 2L))
  ]
  first_loading_start <- sqrt(
    pmax(0.01, abs(within_trait_correlations)) / ratios
  )
  first_loading_cap <- 0.999 / pmax(1, ratios)
  first_loading_start <- pmin(
    0.90 * first_loading_cap,
    first_loading_start
  )
  trait_raw <- atanh(pmin(
    0.95,
    pmax(-0.95, first_loading_start / first_loading_cap)
  ))

  method_raw <- rep(0, length(personality_indicator_indices))
  if (start_number == 2L) {
    method_raw <- rep(c(0.10, -0.05), 5L)
  } else if (start_number >= 3L) {
    method_raw <- stats::rnorm(
      length(personality_indicator_indices),
      sd = 0.12
    )
  }

  if (structure == "correlated") {
    structure_raw <- if (start_number == 1L) {
      rep(0, length(layout$structure))
    } else {
      stats::rnorm(length(layout$structure), sd = 0.12)
    }
  } else {
    general_start <- if (start_number == 1L) 0.55 else 0.35
    structure_raw <- rep(
      atanh(general_start / 0.999),
      length(layout$structure)
    )
    if (start_number >= 3L) {
      structure_raw <- structure_raw +
        stats::rnorm(length(layout$structure), sd = 0.15)
    }
  }

  ratio_raw <- if (loading_model == "common") {
    base <- mean(log(fixed_ratio[trait_names[1:5]]))
    if (start_number == 2L) {
      log(1.05)
    } else if (start_number >= 3L) {
      base + stats::rnorm(1L, sd = 0.10)
    } else {
      base
    }
  } else {
    numeric()
  }

  c(trait_raw, method_raw, structure_raw, ratio_raw)
}

fit_ratio_model <- function(
  analysis,
  analysis_id,
  inputs,
  structure,
  loading_model
) {
  layout <- parameter_layout(structure, loading_model)
  observed_vector <- inputs$correlation_matrix[pair_indices]
  objective <- function(parameters) {
    decoded <- decode_parameters(
      parameters,
      structure,
      loading_model,
      layout
    )
    residual <- observed_vector - decoded$implied_correlation[pair_indices]
    as.numeric(crossprod(
      residual,
      inputs$objective_precision %*% residual
    ))
  }

  set.seed(
    1100L +
      match(analysis_id, input_specifications$analysis_id) * 100L +
      match(structure, c("correlated", "hierarchical")) * 10L +
      match(
        loading_model,
        c("equal", "common", "fixed_chisq", "fixed_magma_chisq")
      )
  )
  starts <- lapply(
    seq_len(4L),
    function(start_number) {
      make_start(
        inputs$correlation_matrix,
        structure,
        loading_model,
        layout,
        start_number
      )
    }
  )
  lower <- rep(-10, layout$count)
  upper <- rep(10, layout$count)
  if (length(layout$ratio) > 0L) {
    lower[layout$ratio] <- log(0.35)
    upper[layout$ratio] <- log(2.85)
  }
  optimizers <- lapply(
    starts,
    function(start) {
      nlminb(
        start = pmax(lower, pmin(upper, start)),
        objective = objective,
        lower = lower,
        upper = upper,
        control = list(iter.max = 6000L, eval.max = 15000L)
      )
    }
  )
  objectives <- vapply(optimizers, `[[`, numeric(1), "objective")
  optimizer <- optimizers[[which.min(objectives)]]
  refinement <- nlminb(
    start = optimizer$par,
    objective = objective,
    lower = lower,
    upper = upper,
    control = list(
      iter.max = 20000L,
      eval.max = 50000L,
      rel.tol = 1e-12,
      x.tol = 1e-10
    )
  )
  if (is.finite(refinement$objective) &&
      refinement$objective <= optimizer$objective + 1e-9) {
    optimizer <- refinement
  }

  parameters <- optimizer$par
  decoded <- decode_parameters(
    parameters,
    structure,
    loading_model,
    layout
  )
  objective_gradient <- numDeriv::grad(objective, parameters)
  maximum_gradient <- max(abs(objective_gradient))
  projected_gradient <- objective_gradient
  at_lower_bound <- parameters <= lower + 1e-6
  at_upper_bound <- parameters >= upper - 1e-6
  projected_gradient[at_lower_bound & objective_gradient > 0] <- 0
  projected_gradient[at_upper_bound & objective_gradient < 0] <- 0
  maximum_projected_gradient <- max(abs(projected_gradient))
  moment_function <- function(value) {
    decode_parameters(
      value,
      structure,
      loading_model,
      layout
    )$implied_correlation[pair_indices]
  }
  moment_jacobian <- numDeriv::jacobian(moment_function, parameters)
  information <- crossprod(
    moment_jacobian,
    inputs$sampling_precision %*% moment_jacobian
  )
  information_eigenvalues <- eigen(
    information,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  information_tolerance <- max(information_eigenvalues) *
    max(dim(information)) * .Machine$double.eps^0.75
  information_rank <- sum(information_eigenvalues > information_tolerance)
  full_rank <- information_rank == layout$count &&
    min(information_eigenvalues) > 0
  parameter_vcov <- if (full_rank) {
    solve(information)
  } else {
    matrix(NA_real_, nrow = layout$count, ncol = layout$count)
  }

  loading_function <- function(value) {
    decode_parameters(
      value,
      structure,
      loading_model,
      layout
    )$trait_loadings
  }
  loading_estimates <- loading_function(parameters)
  loading_jacobian <- numDeriv::jacobian(loading_function, parameters)
  loading_se <- if (full_rank) {
    sqrt(pmax(0, diag(
      loading_jacobian %*% parameter_vcov %*% t(loading_jacobian)
    )))
  } else {
    rep(NA_real_, length(loading_estimates))
  }

  ratio_function <- function(value) {
    decode_parameters(
      value,
      structure,
      loading_model,
      layout
    )$ratios
  }
  ratio_estimates <- ratio_function(parameters)
  ratio_se <- if (loading_model == "common" && full_rank) {
    ratio_jacobian <- numDeriv::jacobian(ratio_function, parameters)
    sqrt(pmax(0, diag(
      ratio_jacobian %*% parameter_vcov %*% t(ratio_jacobian)
    )))
  } else {
    rep(NA_real_, length(ratio_estimates))
  }
  if (loading_model == "common") {
    # IQ remains an equality constraint in this specification; it is not an
    # estimated zero-variance ratio parameter.
    ratio_se[[length(trait_ids)]] <- NA_real_
  }

  residual <- observed_vector - decoded$implied_correlation[pair_indices]
  degrees_of_freedom <- length(observed_vector) - layout$count
  ratio_z <- (ratio_estimates - 1) / ratio_se
  ratio_p <- 2 * pnorm(-abs(ratio_z))
  minimum_trait_eigenvalue <- min(eigen(
    decoded$trait_correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values)
  # nlminb code 1 can be a terminal constrained solution when its internal
  # Hessian is singular at a boundary. Retain it only when the projected KKT
  # gradient passes and the decoded covariance matrices remain admissible.
  terminal_optimizer_solution <- optimizer$convergence %in% c(0L, 1L)
  admissible <- terminal_optimizer_solution &&
    maximum_projected_gradient < 1e-4 &&
    min(decoded$residual_variances) > 0 &&
    minimum_trait_eigenvalue > 0
  near_boundary <- min(decoded$residual_variances) < 0.01 ||
    minimum_trait_eigenvalue < 0.005
  weak_identification <- !full_rank
  gradient_check_passed <- maximum_projected_gradient < 1e-4
  inference_reliable <- admissible &&
    !near_boundary &&
    !weak_identification &&
    gradient_check_passed
  diagnostic_parts <- character()
  if (!admissible) {
    diagnostic_parts <- c(diagnostic_parts, "inadmissible covariance solution")
  }
  if (near_boundary) {
    diagnostic_parts <- c(diagnostic_parts, "near boundary")
  }
  if (weak_identification) {
    diagnostic_parts <- c(
      diagnostic_parts,
      paste0(
        "information rank ",
        information_rank,
        "/",
        layout$count
      )
    )
  }
  if (!gradient_check_passed) {
    diagnostic_parts <- c(diagnostic_parts, "projected gradient above 1e-4")
  }
  diagnostic_status <- if (length(diagnostic_parts) == 0L) {
    "stable"
  } else {
    paste(diagnostic_parts, collapse = "; ")
  }

  rows <- lapply(seq_along(trait_ids), function(index) {
    data.frame(
      analysis = analysis,
      analysis_id = analysis_id,
      structure = structure,
      loading_model = loading_model,
      trait = trait_names[[index]],
      indicator_1 = indicator_map$indicator_1[[index]],
      indicator_2 = indicator_map$indicator_2[[index]],
      indicator_1_label = indicator_map$indicator_1_label[[index]],
      indicator_2_label = indicator_map$indicator_2_label[[index]],
      mean_chisq_1 = mean_chisq_1[[trait_names[[index]]]],
      mean_chisq_2 = mean_chisq_2[[trait_names[[index]]]],
      chisq_signal_1 = chisq_signal_1[[trait_names[[index]]]],
      chisq_signal_2 = chisq_signal_2[[trait_names[[index]]]],
      kappa_raw_excess_chisq =
        raw_excess_chisq_ratio[[trait_names[[index]]]],
      kappa_correlation_scale = fixed_ratio[[trait_names[[index]]]],
      magma_common_gene_count = magma_calibration$common_gene_count[[index]],
      magma_mean_chisq_1 = magma_mean_chisq_1[[trait_names[[index]]]],
      magma_mean_chisq_2 = magma_mean_chisq_2[[trait_names[[index]]]],
      magma_chisq_signal_1 =
        magma_chisq_signal_1[[trait_names[[index]]]],
      magma_chisq_signal_2 =
        magma_chisq_signal_2[[trait_names[[index]]]],
      magma_kappa_raw_excess_chisq =
        magma_raw_excess_chisq_ratio[[trait_names[[index]]]],
      magma_kappa_correlation_scale =
        magma_fixed_ratio[[trait_names[[index]]]],
      loading_ratio = ratio_estimates[[index]],
      loading_ratio_se = ratio_se[[index]],
      loading_ratio_ci_lower = ratio_estimates[[index]] -
        qnorm(0.975) * ratio_se[[index]],
      loading_ratio_ci_upper = ratio_estimates[[index]] +
        qnorm(0.975) * ratio_se[[index]],
      loading_ratio_z_vs_1 = ratio_z[[index]],
      loading_ratio_p_vs_1 = ratio_p[[index]],
      loading_1 = loading_estimates[[2L * index - 1L]],
      loading_1_se = loading_se[[2L * index - 1L]],
      loading_2 = loading_estimates[[2L * index]],
      loading_2_se = loading_se[[2L * index]],
      chi_square = optimizer$objective,
      degrees_of_freedom = degrees_of_freedom,
      fit_p_value = pchisq(
        optimizer$objective,
        degrees_of_freedom,
        lower.tail = FALSE
      ),
      off_diagonal_rmsr = sqrt(mean(residual^2)),
      free_parameters = layout$count,
      optimizer_convergence_code = optimizer$convergence,
      optimizer_message = optimizer$message,
      maximum_absolute_gradient = maximum_gradient,
      maximum_absolute_projected_gradient = maximum_projected_gradient,
      gradient_check_passed = gradient_check_passed,
      information_rank = information_rank,
      information_dimension = layout$count,
      minimum_information_eigenvalue = min(information_eigenvalues),
      minimum_residual_variance = min(decoded$residual_variances),
      minimum_trait_correlation_eigenvalue = minimum_trait_eigenvalue,
      near_boundary = near_boundary,
      weak_identification = weak_identification,
      admissible = admissible,
      inference_reliable = inference_reliable,
      diagnostic_status = diagnostic_status,
      stringsAsFactors = FALSE
    )
  })

  message(
    sprintf(
      "%s | %s | %s: chi-square = %.3f, df = %d, ratio(s) = %s",
      analysis,
      structure,
      loading_model,
      optimizer$objective,
      degrees_of_freedom,
      paste(format(round(unique(ratio_estimates), 4), nsmall = 4), collapse = ", ")
    )
  )

  list(
    rows = do.call(rbind, rows),
    objective = optimizer$objective,
    df = degrees_of_freedom,
    parameters = parameters,
    parameter_vcov = parameter_vcov,
    decoded = decoded,
    admissible = admissible
  )
}

all_fits <- list()
result_rows <- list()
result_index <- 0L
for (analysis_index in seq_len(nrow(input_specifications))) {
  input_row <- input_specifications[analysis_index, , drop = FALSE]
  inputs <- load_analysis_inputs(input_row)
  for (structure in c("correlated", "hierarchical")) {
    for (loading_model in c(
      "equal",
      "common",
      "fixed_chisq",
      "fixed_magma_chisq"
    )) {
      fit_key <- paste(
        input_row$analysis_id[[1L]],
        structure,
        loading_model,
        sep = "__"
      )
      result_index <- result_index + 1L
      all_fits[[fit_key]] <- fit_ratio_model(
        input_row$analysis[[1L]],
        input_row$analysis_id[[1L]],
        inputs,
        structure,
        loading_model
      )
      result_rows[[result_index]] <- all_fits[[fit_key]]$rows
    }
  }
}
results <- do.call(rbind, result_rows)
rownames(results) <- NULL
write.csv(results, results_path, row.names = FALSE, quote = TRUE, na = "")

summarize_oriented_trait_correlations <- function(
  parameters,
  structure,
  loading_model,
  orientation
) {
  layout <- parameter_layout(structure, loading_model)
  trait_correlation <- decode_parameters(
    parameters,
    structure,
    loading_model,
    layout
  )$trait_correlation
  trait_correlation <- trait_correlation * outer(orientation, orientation)
  big_five_pair_indices <- t(combn(seq_len(5L), 2L))
  mean_big_five <- mean(trait_correlation[big_five_pair_indices])
  mean_iq <- mean(trait_correlation[seq_len(5L), 6L])
  c(
    mean_big_five_big_five = mean_big_five,
    mean_iq_big_five = mean_iq,
    difference_iq_minus_big_five = mean_iq - mean_big_five
  )
}

mean_correlation_rows <- list()
mean_correlation_index <- 0L
for (analysis_index in seq_len(nrow(input_specifications))) {
  analysis_id <- input_specifications$analysis_id[[analysis_index]]
  analysis <- input_specifications$analysis[[analysis_index]]
  for (structure in c("correlated", "hierarchical")) {
    for (loading_model in c(
      "equal",
      "common",
      "fixed_chisq",
      "fixed_magma_chisq"
    )) {
      fit_key <- paste(
        analysis_id,
        structure,
        loading_model,
        sep = "__"
      )
      fitted_model <- all_fits[[fit_key]]
      orientation <- sign(
        fitted_model$decoded$first_loadings +
          fitted_model$decoded$second_loadings
      )
      orientation[orientation == 0] <- 1
      estimates <- summarize_oriented_trait_correlations(
        fitted_model$parameters,
        structure,
        loading_model,
        orientation
      )
      if (all(is.finite(fitted_model$parameter_vcov))) {
        summary_jacobian <- numDeriv::jacobian(
          function(parameters) {
            summarize_oriented_trait_correlations(
              parameters,
              structure,
              loading_model,
              orientation
            )
          },
          fitted_model$parameters
        )
        summary_vcov <-
          summary_jacobian %*%
            fitted_model$parameter_vcov %*%
            t(summary_jacobian)
        standard_errors <- sqrt(pmax(0, diag(summary_vcov)))
      } else {
        standard_errors <- rep(NA_real_, length(estimates))
      }
      difference_z <-
        estimates[["difference_iq_minus_big_five"]] /
          standard_errors[[3L]]
      mean_correlation_index <- mean_correlation_index + 1L
      mean_correlation_rows[[mean_correlation_index]] <- data.frame(
        analysis = analysis,
        analysis_id = analysis_id,
        structure = structure,
        loading_model = loading_model,
        mean_big_five_big_five =
          estimates[["mean_big_five_big_five"]],
        se_big_five_big_five = standard_errors[[1L]],
        ci_big_five_big_five_lower =
          estimates[["mean_big_five_big_five"]] -
            qnorm(0.975) * standard_errors[[1L]],
        ci_big_five_big_five_upper =
          estimates[["mean_big_five_big_five"]] +
            qnorm(0.975) * standard_errors[[1L]],
        mean_iq_big_five = estimates[["mean_iq_big_five"]],
        se_iq_big_five = standard_errors[[2L]],
        ci_iq_big_five_lower =
          estimates[["mean_iq_big_five"]] -
            qnorm(0.975) * standard_errors[[2L]],
        ci_iq_big_five_upper =
          estimates[["mean_iq_big_five"]] +
            qnorm(0.975) * standard_errors[[2L]],
        difference_iq_minus_big_five =
          estimates[["difference_iq_minus_big_five"]],
        se_difference = standard_errors[[3L]],
        ci_difference_lower =
          estimates[["difference_iq_minus_big_five"]] -
            qnorm(0.975) * standard_errors[[3L]],
        ci_difference_upper =
          estimates[["difference_iq_minus_big_five"]] +
            qnorm(0.975) * standard_errors[[3L]],
        z_equal_means = difference_z,
        p_equal_means = 2 * pnorm(-abs(difference_z)),
        diagnostic_status = fitted_model$rows$diagnostic_status[[1L]],
        stringsAsFactors = FALSE
      )
    }
  }
}
mean_correlations <- do.call(rbind, mean_correlation_rows)
rownames(mean_correlations) <- NULL
mean_correlations$change_big_five_from_equal <- ave(
  mean_correlations$mean_big_five_big_five,
  interaction(
    mean_correlations$analysis_id,
    mean_correlations$structure,
    drop = TRUE
  ),
  FUN = function(values) values - values[[1L]]
)
mean_correlations$change_iq_from_equal <- ave(
  mean_correlations$mean_iq_big_five,
  interaction(
    mean_correlations$analysis_id,
    mean_correlations$structure,
    drop = TRUE
  ),
  FUN = function(values) values - values[[1L]]
)
write.csv(
  mean_correlations,
  mean_correlations_path,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

comparison_rows <- list()
comparison_index <- 0L
fixed_model_labels <- c(
  fixed_chisq = "GWAS chi-square",
  fixed_magma_chisq = "MAGMA Z-squared"
)
for (analysis_index in seq_len(nrow(input_specifications))) {
  analysis_id <- input_specifications$analysis_id[[analysis_index]]
  analysis <- input_specifications$analysis[[analysis_index]]
  for (structure in c("correlated", "hierarchical")) {
    equal_fit <- all_fits[[paste(analysis_id, structure, "equal", sep = "__")]]
    common_fit <- all_fits[[paste(analysis_id, structure, "common", sep = "__")]]
    comparison_index <- comparison_index + 1L
    comparison_rows[[comparison_index]] <- data.frame(
      analysis = analysis,
      analysis_id = analysis_id,
      structure = structure,
      comparison = "Equal vs estimated common personality ratio",
      reference_model = "equal",
      alternative_model = "common",
      reference_chi_square = equal_fit$objective,
      alternative_chi_square = common_fit$objective,
      chi_square_difference = equal_fit$objective - common_fit$objective,
      df_difference = equal_fit$df - common_fit$df,
      p_value = pchisq(
        equal_fit$objective - common_fit$objective,
        equal_fit$df - common_fit$df,
        lower.tail = FALSE
      ),
      comparison_type = "Nested likelihood-style WLS difference test",
      stringsAsFactors = FALSE
    )
    for (fixed_model in names(fixed_model_labels)) {
      fixed_fit <- all_fits[[paste(
        analysis_id,
        structure,
        fixed_model,
        sep = "__"
      )]]
      calibration_label <- fixed_model_labels[[fixed_model]]
      comparison_index <- comparison_index + 1L
      comparison_rows[[comparison_index]] <- data.frame(
        analysis = analysis,
        analysis_id = analysis_id,
        structure = structure,
        comparison = paste0(
          calibration_label,
          " ratios vs equal ratios"
        ),
        reference_model = "equal",
        alternative_model = fixed_model,
        reference_chi_square = equal_fit$objective,
        alternative_chi_square = fixed_fit$objective,
        chi_square_difference = equal_fit$objective - fixed_fit$objective,
        df_difference = equal_fit$df - fixed_fit$df,
        p_value = NA_real_,
        comparison_type = paste0(
          "Descriptive only: same parameter count and non-nested fixed values"
        ),
        stringsAsFactors = FALSE
      )
      comparison_index <- comparison_index + 1L
      comparison_rows[[comparison_index]] <- data.frame(
        analysis = analysis,
        analysis_id = analysis_id,
        structure = structure,
        comparison = paste0(
          "Estimated common personality ratio vs ",
          tolower(calibration_label),
          " ratios"
        ),
        reference_model = fixed_model,
        alternative_model = "common",
        reference_chi_square = fixed_fit$objective,
        alternative_chi_square = common_fit$objective,
        chi_square_difference = fixed_fit$objective - common_fit$objective,
        df_difference = fixed_fit$df - common_fit$df,
        p_value = NA_real_,
        comparison_type = "Descriptive only: models impose different constraints",
        stringsAsFactors = FALSE
      )
    }
  }
}
comparisons <- do.call(rbind, comparison_rows)
rownames(comparisons) <- NULL
write.csv(
  comparisons,
  comparisons_path,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

format_number <- function(x, digits = 3L) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}
format_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "< .001", sub("^0", "", formatC(x, digits = 3L, format = "f")))
  )
}

loading_model_labels <- c(
  equal = "Equal",
  common = "Common ratio",
  fixed_chisq = "Fixed GWAS chi-square",
  fixed_magma_chisq = "Fixed MAGMA Z-squared"
)

calibration_table_lines <- c(
  "| Trait | GWAS mean chi-square 1 | GWAS mean chi-square 2 | GWAS kappa | MAGMA mean Z-squared 1 | MAGMA mean Z-squared 2 | MAGMA kappa | Common genes |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_along(trait_names), function(index) {
    trait <- trait_names[[index]]
    paste0(
      "| ", trait,
      " | ", format_number(mean_chisq_1[[trait]], 6L),
      " | ", format_number(mean_chisq_2[[trait]], 6L),
      " | ", format_number(fixed_ratio[[trait]], 3L),
      " | ", format_number(magma_mean_chisq_1[[trait]], 6L),
      " | ", format_number(magma_mean_chisq_2[[trait]], 6L),
      " | ", format_number(magma_fixed_ratio[[trait]], 3L),
      " | ", magma_calibration$common_gene_count[[index]],
      " |"
    )
  }, character(1))
)

fit_unique <- results[!duplicated(results[c(
  "analysis", "structure", "loading_model"
)]), ]
fit_table_lines <- c(
  "| Analysis | Structure | Loading model | Chi-square | df | p | RMSR | Diagnostic status |",
  "|---|---|---|---:|---:|---:|---:|---|",
  vapply(seq_len(nrow(fit_unique)), function(index) {
    row <- fit_unique[index, ]
    paste0(
      "| ", row$analysis,
      " | ", row$structure,
      " | ", unname(loading_model_labels[[row$loading_model]]),
      " | ", format_number(row$chi_square),
      " | ", row$degrees_of_freedom,
      " | ", format_p(row$fit_p_value),
      " | ", format_number(row$off_diagonal_rmsr),
      " | ", row$diagnostic_status,
      " |"
    )
  }, character(1))
)

common_rows <- results[
  results$loading_model == "common" &
    results$trait %in% c("Agreeableness", "IQ"),
]
common_table_lines <- c(
  "| Analysis | Structure | Ratio | Estimate | SE | 95% CI | p vs 1 |",
  "|---|---|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(common_rows)), function(index) {
    row <- common_rows[index, ]
    ratio_label <- if (row$trait == "IQ") {
      "Male / female IQ (fixed equal)"
    } else {
      "Personality sample 2 / sample 1"
    }
    paste0(
      "| ", row$analysis,
      " | ", row$structure,
      " | ", ratio_label,
      " | ", format_number(row$loading_ratio),
      " | ", format_number(row$loading_ratio_se),
      " | [", format_number(row$loading_ratio_ci_lower),
      ", ", format_number(row$loading_ratio_ci_upper), "]",
      " | ", format_p(row$loading_ratio_p_vs_1),
      " |"
    )
  }, character(1))
)

format_estimate_se <- function(estimate, standard_error) {
  if (is.na(standard_error)) {
    paste0(format_number(estimate), " (SE unavailable)")
  } else {
    paste0(
      format_number(estimate),
      " (",
      format_number(standard_error),
      ")"
    )
  }
}
mean_correlation_table_lines <- c(
  paste0(
    "| Analysis | Structure | Loading model | Mean Big Five-Big Five (SE) | ",
    "Mean IQ-Big Five (SE) | Difference: IQ minus Big Five | p, equal means |"
  ),
  "|---|---|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(mean_correlations)), function(index) {
    row <- mean_correlations[index, ]
    paste0(
      "| ", row$analysis,
      " | ", row$structure,
      " | ", unname(loading_model_labels[[row$loading_model]]),
      " | ", format_estimate_se(
        row$mean_big_five_big_five,
        row$se_big_five_big_five
      ),
      " | ", format_estimate_se(
        row$mean_iq_big_five,
        row$se_iq_big_five
      ),
      " | ", format_number(row$difference_iq_minus_big_five),
      " | ", format_p(row$p_equal_means),
      " |"
    )
  }, character(1))
)

nested_rows <- comparisons[
  comparisons$comparison == "Equal vs estimated common personality ratio",
]
nested_table_lines <- c(
  "| Analysis | Structure | Delta chi-square | Delta df | p |",
  "|---|---|---:|---:|---:|",
  vapply(seq_len(nrow(nested_rows)), function(index) {
    row <- nested_rows[index, ]
    paste0(
      "| ", row$analysis,
      " | ", row$structure,
      " | ", format_number(row$chi_square_difference),
      " | ", row$df_difference,
      " | ", format_p(row$p_value),
      " |"
    )
  }, character(1))
)

descriptive_rows <- comparisons[
  grepl("ratios vs equal ratios$", comparisons$comparison),
]
descriptive_table_lines <- c(
  "| Analysis | Structure | Calibration | Chi-square(equal) - Chi-square(fixed) |",
  "|---|---|---|---:|",
  vapply(seq_len(nrow(descriptive_rows)), function(index) {
    row <- descriptive_rows[index, ]
    paste0(
      "| ", row$analysis,
      " | ", row$structure,
      " | ", unname(loading_model_labels[[row$alternative_model]]),
      " | ", format_number(row$chi_square_difference),
      " |"
    )
  }, character(1))
)

report_lines <- c(
  "# Split-half loading-ratio sensitivity analysis",
  "",
  paste0("Generated by `scripts/split_half_ratio_models.R` on ", Sys.Date(), "."),
  "",
  "## Question and models",
  "",
  paste0(
    "This standalone analysis asks whether the paired trait loadings should be ",
    "proportional rather than exactly equal. It fits both the correlated-trait ",
    "and hierarchical general-factor structures to the gene-score and pooled ",
    "gene-set correlation inputs. The canonical repository models and outputs ",
    "are not altered."
  ),
  "",
  "Four loading specifications are reported:",
  "",
  "1. **Equal (reference):** every half-2/half-1 loading ratio is one.",
  paste0(
      "2. **Estimated common ratio:** one ratio is estimated jointly for the ",
    "five personality traits; the female and male IQ loadings remain equal."
  ),
  paste0(
    "3. **Fixed GWAS chi-square ratios:** each trait's ratio is fixed from the ",
    "supplied mean SNP-level GWAS chi-square values."
  ),
  paste0(
    "4. **Fixed MAGMA Z-squared ratios:** each trait's ratio is fixed from the ",
    "mean squared MAGMA gene-level `ZSTAT` values in the common genes for the ",
    "paired indicators. Absolute trait loadings, sample-factor loadings, and ",
    "structural parameters remain estimated in both fixed specifications."
  ),
  "",
  paste0(
    "The IQ equality constraint applies to the estimated-common-ratio model. ",
    "The fixed GWAS and fixed MAGMA specifications use their respective female-",
    "male calibrations (GWAS kappa = ",
    format_number(fixed_ratio[["IQ"]]),
    "; MAGMA kappa = ",
    format_number(magma_fixed_ratio[["IQ"]]),
    ")."
  ),
  "",
  "## Choosing kappa from mean chi-square",
  "",
  paste0(
    "Because the SEM is fitted to correlations, both fixed models use ",
    "`R_h = (mean_chisq_h - 1) / mean_chisq_h` and ",
    "`kappa = sqrt(R_2 / R_1)`, where `lambda_2 = kappa * lambda_1`. ",
    "The square root is required because reliability-like signal enters a ",
    "loading as a squared quantity. The proposed `(chisq_1 - 1) / ",
    "(chisq_2 - 1)` is reversed for a half-2/half-1 ratio, omits that square ",
    "root, and does not account for standardizing each score to variance one."
  ),
  "",
  paste0(
    "For the GWAS calibration, mean chi-square is supplied externally. For the ",
    "MAGMA calibration, `chisq` means `ZSTAT^2` from column 9 of the MAGMA ",
    "`.genes.raw` files; the raw MAGMA gene test statistic is not averaged because ",
    "its scale and effective degrees of freedom vary across genes. The raw-excess ",
    "factors are retained in the calibration CSV for transparency but are not ",
    "imposed in the correlation model."
  ),
  "",
  calibration_table_lines,
  "",
  "## Model fit",
  "",
  fit_table_lines,
  "",
  paste0(
    "RMSR is the root mean squared residual across the 66 off-diagonal ",
    "correlations. The WLS chi-square uses the full delete-one-block ",
    "jackknife sampling covariance and the repository's 199/200 test scaling."
  ),
  "",
  paste0(
    "`Stable` means the covariance solution is admissible, the projected ",
    "optimizer gradient is below 1e-4, the parameter information is full rank, ",
    "and neither the observed residual variances nor the trait-correlation ",
    "eigenvalues are within the stated boundary thresholds. A near-boundary ",
    "point estimate can still be inspected, but its conventional standard ",
    "errors and chi-square reference distribution require caution."
  ),
  "",
  "## Estimated common ratio",
  "",
  common_table_lines,
  "",
  paste0(
    "The ratio standard errors use the delta method with the full jackknife ",
    "sampling covariance. The five personality rows share a single parameter; ",
    "only Agreeableness is shown as its representative in this table. The IQ ",
    "ratio is fixed to one and therefore has no standard error or test."
  ),
  "",
  "## Model-implied mean trait correlations",
  "",
  paste0(
    "For each fitted model, the table averages the 10 latent correlations ",
    "among the Big Five and separately averages the 5 correlations between IQ ",
    "and the Big Five. Trait factors are oriented so their paired indicator ",
    "loadings are positive. The difference is `mean(IQ-Big Five) - ",
    "mean(Big Five-Big Five)`, so a positive value indicates greater average ",
    "similarity to IQ."
  ),
  "",
  mean_correlation_table_lines,
  "",
  paste0(
    "Standard errors and the equality test use the delta method with the full ",
    "jackknife sampling covariance. Fixed calibration ratios are treated as ",
    "known, so their external calibration uncertainty is not included in these ",
    "standard errors. Missing standard errors occur where the parameter ",
    "information is rank deficient."
  ),
  paste0(
    "The remaining correlated gene-set standard errors are finite, but those ",
    "models are also near a trait-correlation boundary and should be interpreted ",
    "as sensitivity diagnostics rather than definitive inference."
  ),
  "",
  paste0(
    "The qualitative conclusion is stable. In the gene-score correlated model, ",
    "the equal-loading specification gives a modestly larger IQ-Big Five mean ",
    "(.433 versus .361; p = .036), but the difference becomes smaller and ",
    "nonsignificant under the estimated-common, fixed-GWAS, and fixed-MAGMA ",
    "calibrations. In the hierarchical gene-score models the two means are ",
    "already very similar. The gene-set sensitivity models raise the absolute ",
    "latent correlations, especially in the hierarchical model, without ",
    "producing evidence that the two mean correlation categories differ. Thus, ",
    "the substantive conclusion is not dependent on the choice between the two ",
    "external chi-square calibrations."
  ),
  "",
  "## Comparisons",
  "",
  "The common-ratio model nests the equal model with one added parameter:",
  "",
  nested_table_lines,
  "",
  paste0(
    "The fixed and equal models have the same number of parameters and impose ",
    "different point constraints, so their chi-square difference is descriptive ",
    "rather than a standard nested chi-square test:"
  ),
  "",
  descriptive_table_lines,
  "",
  "## Interpretation and recommendation",
  "",
  paste0(
    "Across both structural models, estimating the common personality ratio ",
    "tests whether the five personality sample-2 loadings share a proportional ",
    "departure from their sample-1 loadings while female and male IQ remain ",
    "equal by construction. This relaxation improves fit over exact equality ",
    "by Delta chi-square = 10.75 to 15.53 on 1 df (all p <= .0011). The stable ",
    "common-ratio estimates are 1.45 and 1.66 for the gene-score correlated and ",
    "hierarchical models, respectively, and 1.57 for the gene-set hierarchical ",
    "model."
  ),
  "",
  paste0(
    "The gene-set correlated-traits solutions are near a singular ",
    "trait-correlation boundary, including the equal reference model. The ",
    "common-ratio version additionally has rank-deficient parameter ",
    "information, so its ratio standard errors are intentionally reported as ",
    "missing. The gene-set fixed-ratio hierarchical solution also places one ",
    "sample-factor loading near its residual-variance boundary. These are real ",
    "weak-identification warnings, not ordinary nonconvergence."
  ),
  "",
  paste0(
    "The estimated-common-ratio model is the most useful primary sensitivity ",
    "analysis: it relaxes exact equality parsimoniously and lets the observed ",
    "cross-trait covariance pattern determine the proportionality factor. The ",
    "fixed GWAS and fixed MAGMA models are calibration sensitivity analyses. The ",
    "MAGMA calibration is closer to the gene-level indicators used by the ",
    "gene-score SEM, but mean MAGMA Z-squared also reflects gene-level ",
    "polygenicity, gene mapping, LD, gene size, and outliers. Neither external ",
    "chi-square quantity should be interpreted as pure score reliability."
  ),
  "",
  paste0(
    "A fully trait-specific estimated-ratio model would be the freely estimated ",
    "paired-loading model discussed previously. Although it can be locally ",
    "identified in the full MTMM structure, it is more weakly identified and is ",
    "not needed to answer this narrower sensitivity question."
  ),
  "",
  "## Files",
  "",
  "- `split_half_ratio_model_results.csv`: trait-level ratios, loadings, fit, and diagnostics.",
  "- `split_half_ratio_model_comparisons.csv`: nested and descriptive comparisons.",
  "- `split_half_ratio_model_calibrations.csv`: GWAS and MAGMA calibration inputs and ratios.",
  "- `split_half_ratio_model_mean_correlations.csv`: category means, differences, uncertainty, and equality tests.",
  "- `split_half_ratio_model_report.md`: this report."
)
writeLines(report_lines, report_path, useBytes = TRUE)

if (!all(results$inference_reliable)) {
  warning(
    "One or more fitted models has an inference diagnostic warning; see results.",
    call. = FALSE
  )
}

message("Wrote results: ", results_path)
message("Wrote comparisons: ", comparisons_path)
message("Wrote mean correlations: ", mean_correlations_path)
message("Wrote report: ", report_path)
