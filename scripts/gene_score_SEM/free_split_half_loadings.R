#!/usr/bin/env Rscript

# Standalone sensitivity analysis for freely estimated paired trait loadings.
#
# The canonical SEM constrains the two standardized trait loadings for each
# trait to equality. This script retains the canonical two orthogonal
# sample/method factors and both structural specifications, but estimates all
# 12 trait-to-indicator loadings separately. It fits both the gene-score and
# pooled gene-set correlation matrices with the canonical full-jackknife WLS
# precision matrix. Equal-loading reference models provide exact nested
# comparisons and like-for-like mean latent-correlation benchmarks; the
# loading estimates themselves come from the free-paired-loading models.
#
# This script is intentionally standalone and does not overwrite canonical SEM
# outputs. By default it writes to output/free_split_half_loadings/.

suppressPackageStartupMessages({
  library(numDeriv)
})


# ==============================================================================
# 1. Paths and fixed design information
# ==============================================================================

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
    "scripts/free_split_half_loadings.R",
    mustWork = TRUE
  )
}

project_directory <- dirname(dirname(script_path))
arguments <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(arguments) >= 1L) {
  path.expand(arguments[[1L]])
} else {
  file.path(project_directory, "output", "free_split_half_loadings")
}
if (!grepl("^/", output_directory)) {
  output_directory <- file.path(project_directory, output_directory)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

requested_cores <- if (length(arguments) >= 2L) {
  as.integer(arguments[[2L]])
} else {
  min(4L, max(1L, parallel::detectCores() - 1L))
}
if (is.na(requested_cores) || requested_cores < 1L) {
  stop("The core count must be a positive integer.", call. = FALSE)
}
analysis_cores <- if (.Platform$OS.type == "windows") {
  1L
} else {
  requested_cores
}

jackknife_blocks <- 200L
trait_ids <- c("agree", "consc", "extra", "neurot", "open", "iq")
trait_labels <- c(
  agree = "Agreeableness",
  consc = "Conscientiousness",
  extra = "Extraversion",
  neurot = "Neuroticism",
  open = "Openness",
  iq = "IQ"
)
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
indicator_labels <- c(
  "Agreeableness (half 1)",
  "Agreeableness (half 2)",
  "Conscientiousness (half 1)",
  "Conscientiousness (half 2)",
  "Extraversion (half 1)",
  "Extraversion (half 2)",
  "Neuroticism (half 1)",
  "Neuroticism (half 2)",
  "Openness (half 1)",
  "Openness (half 2)",
  "IQ (female)",
  "IQ (male)"
)
indicator_trait_index <- rep(seq_along(trait_ids), each = 2L)
indicator_half <- rep(1:2, length(trait_ids))
personality_indicator_indices <- seq_len(10L)
indicator_sample <- rep(1:2, 5L)
pair_indices <- t(combn(seq_along(short_names), 2L))
pair_names <- paste(
  short_names[pair_indices[, 1L]],
  short_names[pair_indices[, 2L]],
  sep = "__"
)
trait_pair_indices <- t(combn(seq_along(trait_ids), 2L))

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
  jackknife_path = file.path(
    project_directory,
    "output",
    c("gene_scores", "gene_sets"),
    c(
      "gene_score_jackknife_correlations.csv",
      "all_gene_sets_jackknife_correlations.csv"
    )
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 2. Input helpers
# ==============================================================================

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
    stop(analysis, " correlation diagonal is not one.", call. = FALSE)
  }
  if (min(eigen(
    correlation_matrix,
    symmetric = TRUE,
    only.values = TRUE
  )$values) <= 0) {
    stop(analysis, " correlation matrix is not positive definite.", call. = FALSE)
  }

  source_sampling_vcov <- read_named_matrix(
    input_row$sampling_vcov_path[[1L]],
    paste0(analysis, " jackknife sampling covariance")
  )
  source_pair_names <- vapply(
    seq_len(nrow(pair_indices)),
    function(index) {
      resolve_pair_name(
        long_names[pair_indices[index, 1L]],
        long_names[pair_indices[index, 2L]],
        rownames(source_sampling_vcov)
      )
    },
    character(1)
  )
  if (anyNA(source_pair_names) ||
      !all(source_pair_names %in% colnames(source_sampling_vcov))) {
    stop(
      analysis,
      " sampling covariance lacks one or more expected pairs.",
      call. = FALSE
    )
  }
  sampling_vcov <- source_sampling_vcov[
    source_pair_names,
    source_pair_names,
    drop = FALSE
  ]
  dimnames(sampling_vcov) <- list(pair_names, pair_names)
  if (min(eigen(
    sampling_vcov,
    symmetric = TRUE,
    only.values = TRUE
  )$values) <= 0) {
    stop(
      analysis,
      " jackknife sampling covariance is not positive definite.",
      call. = FALSE
    )
  }

  jackknife_table <- read.csv(
    input_row$jackknife_path[[1L]],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!"block" %in% names(jackknife_table) ||
      nrow(jackknife_table) != jackknife_blocks ||
      anyNA(jackknife_table$block) ||
      anyDuplicated(jackknife_table$block) ||
      !setequal(as.integer(jackknife_table$block), seq_len(jackknife_blocks))) {
    stop(
      analysis,
      " jackknife table must contain blocks 1 through 200.",
      call. = FALSE
    )
  }
  jackknife_table <- jackknife_table[
    match(seq_len(jackknife_blocks), as.integer(jackknife_table$block)),
    ,
    drop = FALSE
  ]
  source_jackknife_pair_names <- vapply(
    seq_len(nrow(pair_indices)),
    function(index) {
      resolve_pair_name(
        long_names[pair_indices[index, 1L]],
        long_names[pair_indices[index, 2L]],
        names(jackknife_table)
      )
    },
    character(1)
  )
  if (anyNA(source_jackknife_pair_names)) {
    stop(
      analysis,
      " jackknife table lacks one or more expected pairs.",
      call. = FALSE
    )
  }
  jackknife_correlations <- as.matrix(
    jackknife_table[, source_jackknife_pair_names, drop = FALSE]
  )
  storage.mode(jackknife_correlations) <- "double"
  dimnames(jackknife_correlations) <- list(
    sprintf("jk%03d", seq_len(jackknife_blocks)),
    pair_names
  )
  if (!all(is.finite(jackknife_correlations))) {
    stop(analysis, " jackknife correlations are non-finite.", call. = FALSE)
  }

  centered <- sweep(
    jackknife_correlations,
    2L,
    colMeans(jackknife_correlations),
    "-"
  )
  recovered_sampling_vcov <- (
    (jackknife_blocks - 1) / jackknife_blocks
  ) * crossprod(centered)
  dimnames(recovered_sampling_vcov) <- dimnames(sampling_vcov)
  if (max(abs(recovered_sampling_vcov - sampling_vcov)) > 1e-10) {
    stop(
      analysis,
      " jackknife replicates do not reproduce the supplied sampling covariance.",
      call. = FALSE
    )
  }

  sampling_precision <- solve(sampling_vcov)
  list(
    correlation_matrix = correlation_matrix,
    observed_vector = correlation_matrix[pair_indices],
    sampling_vcov = sampling_vcov,
    sampling_precision = sampling_precision,
    objective_precision =
      ((jackknife_blocks - 1) / jackknife_blocks) * sampling_precision,
    jackknife_correlations = jackknife_correlations
  )
}


# ==============================================================================
# 3. Constrained covariance models
# ==============================================================================

correlation_matrix_from_parameters <- function(parameters, dimension) {
  lower <- diag(dimension)
  lower[lower.tri(lower)] <- parameters
  covariance <- tcrossprod(lower)
  covariance / sqrt(outer(diag(covariance), diag(covariance)))
}

matrix_to_correlation_parameters <- function(x) {
  eigen_decomposition <- eigen(x, symmetric = TRUE)
  eigenvalues <- pmax(eigen_decomposition$values, 1e-5)
  positive_definite <- eigen_decomposition$vectors %*%
    diag(eigenvalues) %*%
    t(eigen_decomposition$vectors)
  positive_definite <- positive_definite /
    sqrt(outer(diag(positive_definite), diag(positive_definite)))
  lower <- t(chol(positive_definite))
  unit_lower <- sweep(lower, 1L, diag(lower), "/")
  unit_lower[lower.tri(unit_lower)]
}

parameter_layout <- function(structure, loading_model) {
  trait_count <- if (loading_model == "equal") {
    length(trait_ids)
  } else {
    length(short_names)
  }
  structure_count <- if (structure == "correlated") {
    choose(length(trait_ids), 2L)
  } else {
    length(trait_ids)
  }
  starts <- c(
    trait = 1L,
    method = 1L + trait_count,
    structure = 1L + trait_count + length(personality_indicator_indices)
  )
  count <- trait_count +
    length(personality_indicator_indices) +
    structure_count
  list(
    trait = starts[["trait"]] + seq_len(trait_count) - 1L,
    method = starts[["method"]] +
      seq_len(length(personality_indicator_indices)) - 1L,
    structure = starts[["structure"]] + seq_len(structure_count) - 1L,
    count = count
  )
}

decode_parameters <- function(parameters, structure, loading_model, layout) {
  if (loading_model == "equal") {
    paired_loadings <- 0.999 * plogis(parameters[layout$trait])
    trait_loadings <- rep(paired_loadings, each = 2L)
  } else {
    trait_loadings <- 0.999 * plogis(parameters[layout$trait])
  }
  names(trait_loadings) <- short_names

  method_cap <- 0.999 * sqrt(pmax(
    1e-12,
    1 - trait_loadings[personality_indicator_indices]^2
  ))
  method_loadings <- method_cap * tanh(parameters[layout$method])
  names(method_loadings) <- short_names[personality_indicator_indices]
  # Each orthogonal sample factor has an arbitrary global sign. Orient its
  # five loadings to have a nonnegative sum so full-sample and jackknife
  # estimates are reported on a consistent scale.
  for (sample_index in seq_len(2L)) {
    sample_rows <- indicator_sample == sample_index
    if (sum(method_loadings[sample_rows]) < 0) {
      method_loadings[sample_rows] <- -method_loadings[sample_rows]
    }
  }

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

  trait_measurement <- matrix(
    0,
    nrow = length(short_names),
    ncol = length(trait_ids),
    dimnames = list(short_names, trait_ids)
  )
  trait_measurement[
    cbind(seq_along(short_names), indicator_trait_index)
  ] <- trait_loadings

  method_measurement <- matrix(
    0,
    nrow = length(short_names),
    ncol = 2L,
    dimnames = list(short_names, c("sample_one", "sample_two"))
  )
  method_measurement[
    cbind(personality_indicator_indices, indicator_sample)
  ] <- method_loadings

  implied_correlation <-
    trait_measurement %*%
      trait_correlation %*%
      t(trait_measurement) +
    tcrossprod(method_measurement)
  residual_variances <- 1 - diag(implied_correlation)
  diag(implied_correlation) <- 1
  dimnames(implied_correlation) <- list(short_names, short_names)

  list(
    trait_loadings = trait_loadings,
    loading_ratios = trait_loadings[seq(2L, 12L, by = 2L)] /
      trait_loadings[seq(1L, 11L, by = 2L)],
    method_loadings = method_loadings,
    general_loadings = general_loadings,
    trait_correlation = trait_correlation,
    residual_variances = residual_variances,
    implied_correlation = implied_correlation
  )
}

loading_to_raw <- function(x) {
  qlogis(pmin(0.995, pmax(0.005, x / 0.999)))
}

make_start <- function(
  observed_matrix,
  structure,
  loading_model,
  layout,
  start_number
) {
  within_trait <- observed_matrix[
    cbind(seq(1L, 12L, by = 2L), seq(2L, 12L, by = 2L))
  ]
  base_loading <- sqrt(pmax(0.01, abs(within_trait)))
  if (loading_model == "equal") {
    trait_start <- base_loading
  } else {
    ratio_pattern <- switch(
      as.character(start_number),
      `1` = rep(1, 6L),
      `2` = c(rep(1.35, 5L), 0.90),
      `3` = c(1.20, 1.30, 1.25, 1.32, 1.33, 0.90),
      `4` = c(rep(1.60, 5L), 0.85),
      exp(rnorm(6L, mean = log(1.25), sd = 0.20))
    )
    first <- sqrt(pmax(0.005, abs(within_trait)) / ratio_pattern)
    second <- ratio_pattern * first
    scale_down <- pmax(1, pmax(first, second) / 0.85)
    trait_start <- as.vector(rbind(first / scale_down, second / scale_down))
  }
  trait_raw <- loading_to_raw(trait_start)

  method_raw <- rep(0, length(personality_indicator_indices))
  if (start_number == 2L) {
    method_raw <- rep(c(0.15, 0.10), 5L)
  } else if (start_number == 3L) {
    method_raw <- rep(c(0.25, -0.10), 5L)
  } else if (start_number >= 4L) {
    method_raw <- rnorm(length(personality_indicator_indices), sd = 0.18)
  }

  if (structure == "correlated") {
    empirical_trait_correlation <- diag(length(trait_ids))
    for (index in seq_len(nrow(trait_pair_indices))) {
      left <- trait_pair_indices[index, 1L]
      right <- trait_pair_indices[index, 2L]
      cross_values <- c(
        observed_matrix[2L * left - 1L, 2L * right - 1L],
        observed_matrix[2L * left - 1L, 2L * right],
        observed_matrix[2L * left, 2L * right - 1L],
        observed_matrix[2L * left, 2L * right]
      )
      denominator <- base_loading[[left]] * base_loading[[right]]
      estimate <- mean(cross_values) / max(0.02, denominator)
      empirical_trait_correlation[left, right] <-
        empirical_trait_correlation[right, left] <-
          pmin(0.80, pmax(-0.80, estimate))
    }
    if (start_number == 1L) {
      structure_raw <- rep(0, length(layout$structure))
    } else {
      structure_raw <- matrix_to_correlation_parameters(
        empirical_trait_correlation
      )
      if (start_number >= 4L) {
        structure_raw <- structure_raw +
          rnorm(length(structure_raw), sd = 0.10)
      }
    }
  } else {
    general_start <- if (start_number == 1L) 0.55 else 0.70
    structure_raw <- rep(
      atanh(general_start / 0.999),
      length(layout$structure)
    )
    if (start_number >= 3L) {
      structure_raw <- structure_raw +
        rnorm(length(structure_raw), sd = 0.12)
    }
  }

  c(trait_raw, method_raw, structure_raw)
}

fit_covariance_model <- function(
  analysis,
  analysis_id,
  inputs,
  structure,
  loading_model
) {
  layout <- parameter_layout(structure, loading_model)
  objective <- function(parameters, observed_vector = inputs$observed_vector) {
    decoded <- decode_parameters(
      parameters,
      structure,
      loading_model,
      layout
    )
    residual <- observed_vector -
      decoded$implied_correlation[pair_indices]
    as.numeric(crossprod(
      residual,
      inputs$objective_precision %*% residual
    ))
  }

  set.seed(
    3700L +
      match(analysis_id, input_specifications$analysis_id) * 100L +
      match(structure, c("correlated", "hierarchical")) * 10L +
      match(loading_model, c("equal", "free"))
  )
  starts <- lapply(
    seq_len(8L),
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
  fits <- lapply(
    starts,
    function(start) {
      nlminb(
        start = pmin(10, pmax(-10, start)),
        objective = objective,
        lower = rep(-10, layout$count),
        upper = rep(10, layout$count),
        control = list(
          iter.max = 20000L,
          eval.max = 50000L,
          rel.tol = 1e-11,
          x.tol = 1e-9
        )
      )
    }
  )
  finite_objectives <- vapply(
    fits,
    function(fit) if (is.finite(fit$objective)) fit$objective else Inf,
    numeric(1)
  )
  best <- fits[[which.min(finite_objectives)]]
  refinement <- nlminb(
    start = best$par,
    objective = objective,
    lower = rep(-10, layout$count),
    upper = rep(10, layout$count),
    control = list(
      iter.max = 40000L,
      eval.max = 100000L,
      rel.tol = 1e-13,
      x.tol = 1e-11
    )
  )
  candidates <- c(fits, list(refinement))
  candidate_objectives <- vapply(
    candidates,
    function(fit) if (is.finite(fit$objective)) fit$objective else Inf,
    numeric(1)
  )
  best_index <- which.min(candidate_objectives)
  best <- candidates[[best_index]]

  parameters <- best$par
  decoded <- decode_parameters(
    parameters,
    structure,
    loading_model,
    layout
  )
  implied_function <- function(value) {
    decode_parameters(
      value,
      structure,
      loading_model,
      layout
    )$implied_correlation[pair_indices]
  }
  jacobian <- numDeriv::jacobian(implied_function, parameters)
  information <- crossprod(
    jacobian,
    inputs$sampling_precision %*% jacobian
  )
  information_eigenvalues <- eigen(
    information,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  information_tolerance <- max(information_eigenvalues) * 1e-9
  information_rank <- sum(information_eigenvalues > information_tolerance)
  full_rank <- information_rank == layout$count
  parameter_vcov <- if (full_rank) solve(information) else NULL
  raw_gradient <- numDeriv::grad(objective, parameters)
  active_lower <- parameters <= -10 + 1e-5
  active_upper <- parameters >= 10 - 1e-5
  projected_gradient <- raw_gradient
  projected_gradient[active_lower & raw_gradient > 0] <- 0
  projected_gradient[active_upper & raw_gradient < 0] <- 0
  maximum_projected_gradient <- max(abs(projected_gradient))
  trait_eigenvalues <- eigen(
    decoded$trait_correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  # PORT can return code 1 ("false convergence") after taking an extremely
  # small final step in the bounded transformed coordinates. Treat codes 0 and
  # 1 as terminal solutions, but require the explicit projected-gradient,
  # information-rank, and covariance-space checks below.
  optimizer_terminal <- best$convergence %in% c(0L, 1L)
  admissible <- (
    optimizer_terminal &&
      maximum_projected_gradient <= 1e-4 &&
      full_rank &&
      min(decoded$residual_variances) > 1e-8 &&
      min(trait_eigenvalues) > 1e-6
  )

  fitted_residual <- inputs$correlation_matrix -
    decoded$implied_correlation
  chi_square <- best$objective
  degrees_of_freedom <- length(inputs$observed_vector) - layout$count
  null_chi_square <- as.numeric(crossprod(
    inputs$observed_vector,
    inputs$objective_precision %*% inputs$observed_vector
  ))
  null_df <- length(inputs$observed_vector)
  cfi_denominator <- max(
    null_chi_square - null_df,
    chi_square - degrees_of_freedom,
    0
  )
  cfi <- if (cfi_denominator > 0) {
    1 - max(chi_square - degrees_of_freedom, 0) / cfi_denominator
  } else {
    1
  }
  tli_denominator <- null_chi_square / null_df - 1
  tli <- if (abs(tli_denominator) > 1e-12) {
    (null_chi_square / null_df - chi_square / degrees_of_freedom) /
      tli_denominator
  } else {
    NA_real_
  }
  rmsea <- sqrt(max(
    (chi_square - degrees_of_freedom) /
      (degrees_of_freedom * jackknife_blocks),
    0
  ))
  srmr <- sqrt(
    2 * sum(fitted_residual[upper.tri(fitted_residual)]^2) /
      (nrow(fitted_residual) * (nrow(fitted_residual) + 1))
  )

  list(
    analysis = analysis,
    analysis_id = analysis_id,
    structure = structure,
    loading_model = loading_model,
    layout = layout,
    objective = objective,
    parameters = parameters,
    decoded = decoded,
    optimizer = best,
    optimizer_terminal = optimizer_terminal,
    chi_square = chi_square,
    degrees_of_freedom = degrees_of_freedom,
    fit_p_value = pchisq(
      chi_square,
      df = degrees_of_freedom,
      lower.tail = FALSE
    ),
    cfi = cfi,
    tli = tli,
    rmsea = rmsea,
    srmr = srmr,
    residual_matrix = fitted_residual,
    raw_gradient = raw_gradient,
    maximum_projected_gradient = maximum_projected_gradient,
    information = information,
    information_eigenvalues = information_eigenvalues,
    information_rank = information_rank,
    full_rank = full_rank,
    parameter_vcov = parameter_vcov,
    minimum_residual_variance = min(decoded$residual_variances),
    minimum_trait_correlation_eigenvalue = min(trait_eigenvalues),
    admissible = admissible
  )
}


# ==============================================================================
# 4. Full-sample fits and nested comparisons
# ==============================================================================

analysis_inputs <- lapply(
  seq_len(nrow(input_specifications)),
  function(index) load_analysis_inputs(input_specifications[index, ])
)
names(analysis_inputs) <- input_specifications$analysis_id

all_fits <- list()
for (input_index in seq_len(nrow(input_specifications))) {
  specification <- input_specifications[input_index, ]
  inputs <- analysis_inputs[[specification$analysis_id]]
  for (structure in c("correlated", "hierarchical")) {
    for (loading_model in c("equal", "free")) {
      key <- paste(
        specification$analysis_id,
        structure,
        loading_model,
        sep = "__"
      )
      all_fits[[key]] <- fit_covariance_model(
        specification$analysis,
        specification$analysis_id,
        inputs,
        structure,
        loading_model
      )
      message(
        specification$analysis,
        " | ", structure,
        " | ", loading_model,
        ": chi-square = ",
        format(all_fits[[key]]$chi_square, digits = 7),
        ", df = ", all_fits[[key]]$degrees_of_freedom,
        ", admissible = ", all_fits[[key]]$admissible
      )
    }
  }
}

fit_rows <- lapply(
  all_fits,
  function(result) {
    data.frame(
      analysis = result$analysis,
      analysis_id = result$analysis_id,
      structure = result$structure,
      loading_model = result$loading_model,
      chi_square = result$chi_square,
      degrees_of_freedom = result$degrees_of_freedom,
      chi_square_p_value = result$fit_p_value,
      cfi = result$cfi,
      tli = result$tli,
      rmsea = result$rmsea,
      srmr = result$srmr,
      free_parameters = result$layout$count,
      optimizer_convergence_code = result$optimizer$convergence,
      optimizer_message = result$optimizer$message,
      maximum_absolute_projected_gradient =
        result$maximum_projected_gradient,
      information_rank = result$information_rank,
      information_dimension = result$layout$count,
      minimum_information_eigenvalue =
        min(result$information_eigenvalues),
      information_condition_number =
        if (min(result$information_eigenvalues) > 0) {
          max(result$information_eigenvalues) /
            min(result$information_eigenvalues)
        } else {
          NA_real_
        },
      minimum_residual_variance = result$minimum_residual_variance,
      minimum_trait_correlation_eigenvalue =
        result$minimum_trait_correlation_eigenvalue,
      admissible = result$admissible,
      inference_reliable = result$admissible,
      diagnostic_status = if (result$admissible) {
        "stable"
      } else if (!result$full_rank) {
        paste0(
          "weak identification: information rank ",
          result$information_rank,
          "/",
          result$layout$count
        )
      } else if (!result$optimizer_terminal) {
        paste0(
          "optimizer code ",
          result$optimizer$convergence,
          ": ",
          result$optimizer$message
        )
      } else if (result$minimum_residual_variance <= 1e-8 ||
          result$minimum_trait_correlation_eigenvalue <= 1e-6) {
        "inadmissible covariance boundary"
      } else {
        "failed numerical admissibility check"
      },
      stringsAsFactors = FALSE
    )
  }
)
fit_table <- do.call(rbind, fit_rows)
rownames(fit_table) <- NULL

comparison_rows <- list()
comparison_index <- 0L
for (input_index in seq_len(nrow(input_specifications))) {
  analysis_id <- input_specifications$analysis_id[[input_index]]
  analysis <- input_specifications$analysis[[input_index]]
  for (structure in c("correlated", "hierarchical")) {
    equal_result <- all_fits[[paste(
      analysis_id,
      structure,
      "equal",
      sep = "__"
    )]]
    free_result <- all_fits[[paste(
      analysis_id,
      structure,
      "free",
      sep = "__"
    )]]
    delta_chi_square <- equal_result$chi_square - free_result$chi_square
    delta_df <- equal_result$degrees_of_freedom -
      free_result$degrees_of_freedom
    comparison_index <- comparison_index + 1L
    comparison_rows[[comparison_index]] <- data.frame(
      analysis = analysis,
      analysis_id = analysis_id,
      structure = structure,
      comparison = "Equal paired loadings vs all paired loadings free",
      equal_chi_square = equal_result$chi_square,
      free_chi_square = free_result$chi_square,
      chi_square_difference = delta_chi_square,
      degrees_of_freedom_difference = delta_df,
      p_value = if (delta_chi_square >= 0 && delta_df > 0) {
        pchisq(delta_chi_square, df = delta_df, lower.tail = FALSE)
      } else {
        NA_real_
      },
      comparison_type = paste0(
        "Nested fixed-weight WLS difference test for the six paired-loading ",
        "equality restrictions"
      ),
      stringsAsFactors = FALSE
    )
  }
}
comparison_table <- do.call(rbind, comparison_rows)


# ==============================================================================
# 5. Fixed-weight delete-one-block refits of free and equal-loading models
# ==============================================================================

extract_estimands <- function(parameters, structure, loading_model, layout) {
  decoded <- decode_parameters(
    parameters,
    structure,
    loading_model,
    layout
  )
  trait_correlations <- decoded$trait_correlation[trait_pair_indices]
  names(trait_correlations) <- paste(
    trait_ids[trait_pair_indices[, 1L]],
    trait_ids[trait_pair_indices[, 2L]],
    sep = "__"
  )
  big_five_rows <-
    trait_pair_indices[, 1L] <= 5L &
      trait_pair_indices[, 2L] <= 5L
  iq_rows <-
    trait_pair_indices[, 1L] == 6L |
      trait_pair_indices[, 2L] == 6L
  category_means <- c(
    mean_big_five = mean(trait_correlations[big_five_rows]),
    mean_iq_big_five = mean(trait_correlations[iq_rows]),
    difference_iq_minus_big_five =
      mean(trait_correlations[iq_rows]) -
      mean(trait_correlations[big_five_rows])
  )
  general_estimands <- if (structure == "hierarchical") {
    setNames(
      decoded$general_loadings,
      paste0("general_loading__", trait_ids)
    )
  } else {
    numeric()
  }
  c(
    setNames(decoded$trait_loadings, paste0("loading__", short_names)),
    setNames(decoded$loading_ratios, paste0("ratio__", trait_ids)),
    setNames(
      decoded$method_loadings,
      paste0(
        "sample_loading__",
        short_names[personality_indicator_indices]
      )
    ),
    general_estimands,
    setNames(trait_correlations, paste0("correlation__", names(trait_correlations))),
    setNames(category_means, paste0("category__", names(category_means)))
  )
}

refit_jackknife_block <- function(block, full_result, inputs) {
  observed_vector <- inputs$jackknife_correlations[block, ]
  objective <- function(parameters) {
    decoded <- decode_parameters(
      parameters,
      full_result$structure,
      full_result$loading_model,
      full_result$layout
    )
    residual <- observed_vector -
      decoded$implied_correlation[pair_indices]
    as.numeric(crossprod(
      residual,
      inputs$objective_precision %*% residual
    ))
  }
  starts <- list(
    full_result$parameters,
    pmin(
      10,
      pmax(
        -10,
        full_result$parameters +
          sin(seq_along(full_result$parameters) * (block + 1L)) * 0.015
      )
    )
  )
  fits <- lapply(
    starts,
    function(start) {
      tryCatch(
        nlminb(
          start = start,
          objective = objective,
          lower = rep(-10, full_result$layout$count),
          upper = rep(10, full_result$layout$count),
          control = list(
            iter.max = 12000L,
            eval.max = 30000L,
            rel.tol = 1e-10,
            x.tol = 1e-8
          )
        ),
        error = function(condition) condition
      )
    }
  )
  valid <- vapply(
    fits,
    function(value) {
      !inherits(value, "error") && is.finite(value$objective)
    },
    logical(1)
  )
  if (!any(valid)) {
    return(list(
      block = block,
      estimates = rep(
        NA_real_,
        length(extract_estimands(
          full_result$parameters,
          full_result$structure,
          full_result$loading_model,
          full_result$layout
        ))
      ),
      converged = FALSE,
      post_check = FALSE,
      gradient_check_passed = FALSE,
      maximum_absolute_projected_gradient = NA_real_,
      objective = NA_real_,
      optimizer_message = paste(
        vapply(
          fits,
          function(value) {
            if (inherits(value, "error")) conditionMessage(value) else "invalid"
          },
          character(1)
        ),
        collapse = " | "
      )
    ))
  }
  valid_fits <- fits[valid]
  best <- valid_fits[[which.min(vapply(
    valid_fits,
    function(value) value$objective,
    numeric(1)
  ))]]
  refinement <- tryCatch(
    nlminb(
      start = best$par,
      objective = objective,
      lower = rep(-10, full_result$layout$count),
      upper = rep(10, full_result$layout$count),
      control = list(
        iter.max = 24000L,
        eval.max = 60000L,
        rel.tol = 1e-12,
        x.tol = 1e-10
      )
    ),
    error = function(condition) condition
  )
  candidates <- valid_fits
  if (!inherits(refinement, "error") && is.finite(refinement$objective)) {
    candidates <- c(candidates, list(refinement))
  }
  candidate_objectives <- vapply(
    candidates,
    function(value) value$objective,
    numeric(1)
  )
  best_index <- which.min(candidate_objectives)
  best <- candidates[[best_index]]
  decoded <- decode_parameters(
    best$par,
    full_result$structure,
    full_result$loading_model,
    full_result$layout
  )
  raw_gradient <- numDeriv::grad(objective, best$par)
  active_lower <- best$par <= -10 + 1e-5
  active_upper <- best$par >= 10 - 1e-5
  projected_gradient <- raw_gradient
  projected_gradient[active_lower & raw_gradient > 0] <- 0
  projected_gradient[active_upper & raw_gradient < 0] <- 0
  maximum_projected_gradient <- max(abs(projected_gradient))
  minimum_trait_eigenvalue <- min(eigen(
    decoded$trait_correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values)
  post_check <-
    min(decoded$residual_variances) > 1e-8 &&
      minimum_trait_eigenvalue > 1e-6
  gradient_check_passed <- maximum_projected_gradient <= 1e-4
  list(
    block = block,
    estimates = extract_estimands(
      best$par,
      full_result$structure,
      full_result$loading_model,
      full_result$layout
    ),
    converged = best$convergence %in% c(0L, 1L),
    post_check = post_check,
    gradient_check_passed = gradient_check_passed,
    maximum_absolute_projected_gradient = maximum_projected_gradient,
    objective = best$objective,
    optimizer_message = best$message
  )
}

jackknife_results <- list()
for (input_index in seq_len(nrow(input_specifications))) {
  analysis_id <- input_specifications$analysis_id[[input_index]]
  inputs <- analysis_inputs[[analysis_id]]
  for (structure in c("correlated", "hierarchical")) {
    for (loading_model in c("equal", "free")) {
      key <- paste(analysis_id, structure, loading_model, sep = "__")
      full_result <- all_fits[[key]]
      message(
        "Jackknife refits: ", full_result$analysis,
        " | ", structure,
        " | ", loading_model
      )
      block_results <- parallel::mclapply(
        seq_len(jackknife_blocks),
        refit_jackknife_block,
        full_result = full_result,
        inputs = inputs,
        mc.cores = analysis_cores,
        mc.preschedule = TRUE
      )
      estimate_names <- names(block_results[[1L]]$estimates)
      estimate_matrix <- do.call(
        rbind,
        lapply(block_results, function(value) value$estimates)
      )
      colnames(estimate_matrix) <- estimate_names
      diagnostics <- do.call(
        rbind,
        lapply(
          block_results,
          function(value) {
            data.frame(
              analysis = full_result$analysis,
              analysis_id = analysis_id,
              structure = structure,
              loading_model = loading_model,
              block = value$block,
              converged = value$converged,
              post_check = value$post_check,
              gradient_check_passed = value$gradient_check_passed,
              maximum_absolute_projected_gradient =
                value$maximum_absolute_projected_gradient,
              objective = value$objective,
              optimizer_message = value$optimizer_message,
              terminal_covariance_admissible =
                value$converged &&
                  value$post_check &&
                  value$gradient_check_passed,
              stringsAsFactors = FALSE
            )
          }
        )
      )
      if (any(!is.finite(estimate_matrix))) {
        stop(
          "At least one ", full_result$analysis, " ", structure, " ",
          loading_model,
          " jackknife refit failed to return finite estimates.",
          call. = FALSE
        )
      }
      jackknife_results[[key]] <- list(
        estimates = estimate_matrix,
        diagnostics = diagnostics
      )
    }
  }
}

jackknife_se <- function(values) {
  centered <- values - mean(values)
  sqrt(((length(values) - 1) / length(values)) * sum(centered^2))
}


# ==============================================================================
# 6. Estimand tables
# ==============================================================================

loading_rows <- list()
factor_loading_rows <- list()
correlation_rows <- list()
category_rows <- list()
jackknife_diagnostic_rows <- list()
row_index <- 0L

make_factor_loading_row <- function(
  result,
  estimate_name,
  loading_type,
  factor,
  indicator,
  full_estimands,
  replicate_estimates
) {
  estimate <- full_estimands[[estimate_name]]
  standard_error <- jackknife_se(replicate_estimates[, estimate_name])
  reliable <- isTRUE(result$admissible)
  z_value <- if (reliable && standard_error > 0) {
    estimate / standard_error
  } else {
    NA_real_
  }
  data.frame(
    analysis = result$analysis,
    analysis_id = result$analysis_id,
    structure = result$structure,
    loading_model = result$loading_model,
    loading_type = loading_type,
    factor = factor,
    indicator = indicator,
    standardized_loading = estimate,
    standard_error = standard_error,
    se_method = "200-block fixed-weight delete-one jackknife",
    ci_95_lower = estimate - qnorm(0.975) * standard_error,
    ci_95_upper = estimate + qnorm(0.975) * standard_error,
    z = z_value,
    p_value = if (is.finite(z_value)) {
      2 * pnorm(-abs(z_value))
    } else {
      NA_real_
    },
    equality_constraint = "",
    inference_reliable = reliable,
    diagnostic_status = if (reliable) {
      "stable"
    } else if (!result$full_rank) {
      paste0(
        "weak identification: information rank ",
        result$information_rank,
        "/",
        result$layout$count
      )
    } else {
      "failed numerical admissibility check"
    },
    stringsAsFactors = FALSE
  )
}

for (input_index in seq_len(nrow(input_specifications))) {
  analysis_id <- input_specifications$analysis_id[[input_index]]
  for (structure in c("correlated", "hierarchical")) {
    key <- paste(analysis_id, structure, "free", sep = "__")
    jk_key <- key
    result <- all_fits[[key]]
    replicate_estimates <- jackknife_results[[jk_key]]$estimates
    full_estimands <- extract_estimands(
      result$parameters,
      structure,
      "free",
      result$layout
    )

    for (trait_index in seq_along(trait_ids)) {
      first_indicator <- 2L * trait_index - 1L
      second_indicator <- 2L * trait_index
      first_name <- paste0("loading__", short_names[[first_indicator]])
      second_name <- paste0("loading__", short_names[[second_indicator]])
      ratio_name <- paste0("ratio__", trait_ids[[trait_index]])
      first_se <- jackknife_se(replicate_estimates[, first_name])
      second_se <- jackknife_se(replicate_estimates[, second_name])
      ratio_se <- jackknife_se(replicate_estimates[, ratio_name])
      row_index <- row_index + 1L
      loading_rows[[row_index]] <- data.frame(
        analysis = result$analysis,
        analysis_id = analysis_id,
        structure = structure,
        trait = trait_labels[[trait_ids[[trait_index]]]],
        trait_id = trait_ids[[trait_index]],
        indicator_1 = indicator_labels[[first_indicator]],
        indicator_2 = indicator_labels[[second_indicator]],
        loading_1 = full_estimands[[first_name]],
        loading_1_jackknife_se = first_se,
        loading_1_ci_95_lower = full_estimands[[first_name]] -
          qnorm(0.975) * first_se,
        loading_1_ci_95_upper = full_estimands[[first_name]] +
          qnorm(0.975) * first_se,
        loading_2 = full_estimands[[second_name]],
        loading_2_jackknife_se = second_se,
        loading_2_ci_95_lower = full_estimands[[second_name]] -
          qnorm(0.975) * second_se,
        loading_2_ci_95_upper = full_estimands[[second_name]] +
          qnorm(0.975) * second_se,
        loading_ratio_2_over_1 = full_estimands[[ratio_name]],
        loading_ratio_jackknife_se = ratio_se,
        loading_ratio_ci_95_lower = full_estimands[[ratio_name]] -
          qnorm(0.975) * ratio_se,
        loading_ratio_ci_95_upper = full_estimands[[ratio_name]] +
          qnorm(0.975) * ratio_se,
        stringsAsFactors = FALSE
      )
      for (indicator_index in c(first_indicator, second_indicator)) {
        estimate_name <- paste0(
          "loading__",
          short_names[[indicator_index]]
        )
        factor_loading_rows[[length(factor_loading_rows) + 1L]] <-
          make_factor_loading_row(
            result = result,
            estimate_name = estimate_name,
            loading_type = "Trait factor to split-half indicator",
            factor = trait_labels[[trait_ids[[trait_index]]]],
            indicator = indicator_labels[[indicator_index]],
            full_estimands = full_estimands,
            replicate_estimates = replicate_estimates
          )
      }
    }

    for (indicator_index in personality_indicator_indices) {
      sample_index <- indicator_sample[[indicator_index]]
      estimate_name <- paste0(
        "sample_loading__",
        short_names[[indicator_index]]
      )
      factor_loading_rows[[length(factor_loading_rows) + 1L]] <-
        make_factor_loading_row(
          result = result,
          estimate_name = estimate_name,
          loading_type = "Sample-specific factor to indicator",
          factor = paste("Sample-specific factor", sample_index),
          indicator = indicator_labels[[indicator_index]],
          full_estimands = full_estimands,
          replicate_estimates = replicate_estimates
        )
    }

    if (structure == "hierarchical") {
      general_factor_label <- paste(
        "General",
        tolower(result$analysis)
      )
      for (trait_index in seq_along(trait_ids)) {
        estimate_name <- paste0(
          "general_loading__",
          trait_ids[[trait_index]]
        )
        factor_loading_rows[[length(factor_loading_rows) + 1L]] <-
          make_factor_loading_row(
            result = result,
            estimate_name = estimate_name,
            loading_type = "General factor to trait factor",
            factor = general_factor_label,
            indicator = trait_labels[[trait_ids[[trait_index]]]],
            full_estimands = full_estimands,
            replicate_estimates = replicate_estimates
          )
      }
    }

    for (pair_index in seq_len(nrow(trait_pair_indices))) {
      left_id <- trait_ids[trait_pair_indices[pair_index, 1L]]
      right_id <- trait_ids[trait_pair_indices[pair_index, 2L]]
      estimate_name <- paste0(
        "correlation__", left_id, "__", right_id
      )
      estimate <- full_estimands[[estimate_name]]
      standard_error <- jackknife_se(
        replicate_estimates[, estimate_name]
      )
      correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
        analysis = result$analysis,
        analysis_id = analysis_id,
        structure = structure,
        loading_model = "free",
        trait_1 = trait_labels[[left_id]],
        trait_2 = trait_labels[[right_id]],
        correlation = estimate,
        jackknife_se = standard_error,
        ci_95_lower = estimate - qnorm(0.975) * standard_error,
        ci_95_upper = estimate + qnorm(0.975) * standard_error,
        z = estimate / standard_error,
        p_value = 2 * pnorm(-abs(estimate / standard_error)),
        stringsAsFactors = FALSE
      )
    }

    category_definitions <- data.frame(
      statistic = c(
        "Mean Big Five-Big Five correlation",
        "Mean IQ-Big Five correlation",
        "Difference: IQ minus Big Five"
      ),
      estimate_name = paste0(
        "category__",
        c(
          "mean_big_five",
          "mean_iq_big_five",
          "difference_iq_minus_big_five"
        )
      ),
      stringsAsFactors = FALSE
    )
    for (category_index in seq_len(nrow(category_definitions))) {
      estimate_name <- category_definitions$estimate_name[[category_index]]
      estimate <- full_estimands[[estimate_name]]
      standard_error <- jackknife_se(
        replicate_estimates[, estimate_name]
      )
      category_rows[[length(category_rows) + 1L]] <- data.frame(
        analysis = result$analysis,
        analysis_id = analysis_id,
        structure = structure,
        loading_model = "free",
        statistic = category_definitions$statistic[[category_index]],
        estimate = estimate,
        jackknife_se = standard_error,
        ci_95_lower = estimate - qnorm(0.975) * standard_error,
        ci_95_upper = estimate + qnorm(0.975) * standard_error,
        z = estimate / standard_error,
        p_value = 2 * pnorm(-abs(estimate / standard_error)),
        stringsAsFactors = FALSE
      )
    }

    jackknife_diagnostic_rows[[length(jackknife_diagnostic_rows) + 1L]] <-
      jackknife_results[[jk_key]]$diagnostics
  }
}

# Add like-for-like category summaries from the constrained equal-loading
# reference models. These use the same objective weights and the same 200
# delete-one-block refitting procedure as the focal free-loading summaries.
for (input_index in seq_len(nrow(input_specifications))) {
  analysis_id <- input_specifications$analysis_id[[input_index]]
  for (structure in c("correlated", "hierarchical")) {
    key <- paste(analysis_id, structure, "equal", sep = "__")
    result <- all_fits[[key]]
    replicate_estimates <- jackknife_results[[key]]$estimates
    full_estimands <- extract_estimands(
      result$parameters,
      structure,
      "equal",
      result$layout
    )
    for (category_index in seq_len(nrow(category_definitions))) {
      estimate_name <- category_definitions$estimate_name[[category_index]]
      estimate <- full_estimands[[estimate_name]]
      standard_error <- jackknife_se(
        replicate_estimates[, estimate_name]
      )
      category_rows[[length(category_rows) + 1L]] <- data.frame(
        analysis = result$analysis,
        analysis_id = analysis_id,
        structure = structure,
        loading_model = "equal",
        statistic = category_definitions$statistic[[category_index]],
        estimate = estimate,
        jackknife_se = standard_error,
        ci_95_lower = estimate - qnorm(0.975) * standard_error,
        ci_95_upper = estimate + qnorm(0.975) * standard_error,
        z = estimate / standard_error,
        p_value = 2 * pnorm(-abs(estimate / standard_error)),
        stringsAsFactors = FALSE
      )
    }
    jackknife_diagnostic_rows[[length(jackknife_diagnostic_rows) + 1L]] <-
      jackknife_results[[key]]$diagnostics
  }
}

loading_table <- do.call(rbind, loading_rows)
factor_loading_table <- do.call(rbind, factor_loading_rows)
correlation_table <- do.call(rbind, correlation_rows)
category_table <- do.call(rbind, category_rows)
jackknife_diagnostics <- do.call(rbind, jackknife_diagnostic_rows)
rownames(loading_table) <- NULL
rownames(factor_loading_table) <- NULL
rownames(correlation_table) <- NULL
rownames(category_table) <- NULL
rownames(jackknife_diagnostics) <- NULL
factor_loading_type_order <- c(
  "General factor to trait factor",
  "Trait factor to split-half indicator",
  "Sample-specific factor to indicator"
)
factor_loading_table <- factor_loading_table[order(
  match(factor_loading_table$analysis_id, input_specifications$analysis_id),
  match(factor_loading_table$structure, c("correlated", "hierarchical")),
  match(factor_loading_table$loading_type, factor_loading_type_order),
  factor_loading_table$factor,
  factor_loading_table$indicator
), ]
rownames(factor_loading_table) <- NULL
category_table <- category_table[order(
  match(category_table$analysis_id, input_specifications$analysis_id),
  match(category_table$structure, c("correlated", "hierarchical")),
  match(category_table$loading_model, c("equal", "free")),
  match(category_table$statistic, category_definitions$statistic)
), ]
rownames(category_table) <- NULL


# ==============================================================================
# 7. Write outputs and concise report
# ==============================================================================

fit_path <- file.path(output_directory, "free_loading_model_fit.csv")
comparison_path <- file.path(
  output_directory,
  "free_loading_model_comparisons.csv"
)
loading_path <- file.path(
  output_directory,
  "free_loading_estimates.csv"
)
factor_loading_path <- file.path(
  output_directory,
  "free_loading_factor_loadings.csv"
)
correlation_path <- file.path(
  output_directory,
  "free_loading_trait_correlations.csv"
)
category_path <- file.path(
  output_directory,
  "free_loading_correlation_means.csv"
)
jackknife_diagnostic_path <- file.path(
  output_directory,
  "free_loading_jackknife_diagnostics.csv"
)
report_path <- file.path(
  output_directory,
  "free_split_half_loading_report.md"
)

write.csv(fit_table, fit_path, row.names = FALSE, na = "")
write.csv(comparison_table, comparison_path, row.names = FALSE, na = "")
write.csv(loading_table, loading_path, row.names = FALSE, na = "")
write.csv(
  factor_loading_table,
  factor_loading_path,
  row.names = FALSE,
  na = ""
)
write.csv(correlation_table, correlation_path, row.names = FALSE, na = "")
write.csv(category_table, category_path, row.names = FALSE, na = "")
write.csv(
  jackknife_diagnostics,
  jackknife_diagnostic_path,
  row.names = FALSE,
  na = ""
)

format_number <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")
}
format_p <- function(x) {
  ifelse(
    !is.finite(x),
    "NA",
    ifelse(x < 0.001, "< .001", sub("^0", "", formatC(x, 3L, format = "f")))
  )
}
format_estimate_se <- function(estimate, standard_error) {
  paste0(
    format_number(estimate),
    " (", format_number(standard_error), ")"
  )
}

focal_fit <- fit_table[fit_table$loading_model == "free", , drop = FALSE]
fit_lines <- c(
  "| Analysis | Structure | Chi-square | df | p | CFI | TLI | RMSEA | SRMR | Admissible |",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|:---:|",
  vapply(
    seq_len(nrow(focal_fit)),
    function(index) {
      row <- focal_fit[index, ]
      paste0(
        "| ", row$analysis,
        " | ", row$structure,
        " | ", format_number(row$chi_square),
        " | ", row$degrees_of_freedom,
        " | ", format_p(row$chi_square_p_value),
        " | ", format_number(row$cfi),
        " | ", format_number(row$tli),
        " | ", format_number(row$rmsea),
        " | ", format_number(row$srmr),
        " | ", ifelse(row$admissible, "Yes", "No"), " |"
      )
    },
    character(1)
  )
)

comparison_lines <- c(
  "| Analysis | Structure | Delta chi-square | Delta df | p |",
  "|---|---|---:|---:|---:|",
  vapply(
    seq_len(nrow(comparison_table)),
    function(index) {
      row <- comparison_table[index, ]
      paste0(
        "| ", row$analysis,
        " | ", row$structure,
        " | ", format_number(row$chi_square_difference),
        " | ", row$degrees_of_freedom_difference,
        " | ", format_p(row$p_value), " |"
      )
    },
    character(1)
  )
)

loading_lines <- c(
  "| Analysis | Structure | Trait | Loading 1 | Loading 2 | Ratio 2/1 |",
  "|---|---|---|---:|---:|---:|",
  vapply(
    seq_len(nrow(loading_table)),
    function(index) {
      row <- loading_table[index, ]
      paste0(
        "| ", row$analysis,
        " | ", row$structure,
        " | ", row$trait,
        " | ", format_number(row$loading_1),
        " | ", format_number(row$loading_2),
        " | ", format_number(row$loading_ratio_2_over_1), " |"
      )
    },
    character(1)
  )
)

sample_factor_pairs <- loading_table[
  loading_table$trait_id != "iq",
  ,
  drop = FALSE
]
sample_factor_lines <- c(
  paste0(
    "| Analysis | Structure | Trait | Sample factor 1, half 1 (SE) | ",
    "Sample factor 2, half 2 (SE) | Inference |"
  ),
  "|---|---|---|---:|---:|---|",
  vapply(
    seq_len(nrow(sample_factor_pairs)),
    function(index) {
      pair <- sample_factor_pairs[index, ]
      first_row <- factor_loading_table[
        factor_loading_table$analysis_id == pair$analysis_id &
          factor_loading_table$structure == pair$structure &
          factor_loading_table$loading_type ==
            "Sample-specific factor to indicator" &
          factor_loading_table$factor == "Sample-specific factor 1" &
          factor_loading_table$indicator == pair$indicator_1,
        ,
        drop = FALSE
      ]
      second_row <- factor_loading_table[
        factor_loading_table$analysis_id == pair$analysis_id &
          factor_loading_table$structure == pair$structure &
          factor_loading_table$loading_type ==
            "Sample-specific factor to indicator" &
          factor_loading_table$factor == "Sample-specific factor 2" &
          factor_loading_table$indicator == pair$indicator_2,
        ,
        drop = FALSE
      ]
      if (nrow(first_row) != 1L || nrow(second_row) != 1L) {
        stop(
          "Could not match a unique pair of sample-factor loadings for ",
          pair$analysis,
          " ",
          pair$structure,
          " ",
          pair$trait,
          ".",
          call. = FALSE
        )
      }
      paste0(
        "| ", pair$analysis,
        " | ", pair$structure,
        " | ", pair$trait,
        " | ", format_estimate_se(
          first_row$standardized_loading,
          first_row$standard_error
        ),
        " | ", format_estimate_se(
          second_row$standardized_loading,
          second_row$standard_error
        ),
        " | ", ifelse(
          first_row$inference_reliable && second_row$inference_reliable,
          "Reliable",
          "Diagnostic"
        ),
        " |"
      )
    },
    character(1)
  )
)

iq_correlations <- correlation_table[
  correlation_table$trait_1 == "IQ" |
    correlation_table$trait_2 == "IQ",
  ,
  drop = FALSE
]
iq_lines <- c(
  "| Analysis | Structure | Pair | Correlation | Jackknife SE | 95% CI |",
  "|---|---|---|---:|---:|---:|",
  vapply(
    seq_len(nrow(iq_correlations)),
    function(index) {
      row <- iq_correlations[index, ]
      paste0(
        "| ", row$analysis,
        " | ", row$structure,
        " | ", row$trait_1, " - ", row$trait_2,
        " | ", format_number(row$correlation),
        " | ", format_number(row$jackknife_se),
        " | [", format_number(row$ci_95_lower),
        ", ", format_number(row$ci_95_upper), "] |"
      )
    },
    character(1)
  )
)

category_groups <- unique(category_table[, c(
  "analysis",
  "analysis_id",
  "structure",
  "loading_model"
)])
category_summary <- do.call(
  rbind,
  lapply(
    seq_len(nrow(category_groups)),
    function(index) {
      group <- category_groups[index, ]
      rows <- category_table[
        category_table$analysis_id == group$analysis_id &
          category_table$structure == group$structure &
          category_table$loading_model == group$loading_model,
        ,
        drop = FALSE
      ]
      big_five <- rows[
        rows$statistic == "Mean Big Five-Big Five correlation",
        ,
        drop = FALSE
      ]
      iq_big_five <- rows[
        rows$statistic == "Mean IQ-Big Five correlation",
        ,
        drop = FALSE
      ]
      difference <- rows[
        rows$statistic == "Difference: IQ minus Big Five",
        ,
        drop = FALSE
      ]
      fit_row <- fit_table[
        fit_table$analysis_id == group$analysis_id &
          fit_table$structure == group$structure &
          fit_table$loading_model == group$loading_model,
        ,
        drop = FALSE
      ]
      data.frame(
        group,
        mean_big_five = big_five$estimate,
        se_big_five = big_five$jackknife_se,
        mean_iq_big_five = iq_big_five$estimate,
        se_iq_big_five = iq_big_five$jackknife_se,
        difference = difference$estimate,
        se_difference = difference$jackknife_se,
        p_equal_means = difference$p_value,
        inference_reliable = fit_row$inference_reliable,
        stringsAsFactors = FALSE
      )
    }
  )
)

category_lines <- c(
  paste0(
    "| Analysis | Structure | Loadings | Mean Big Five-Big Five (SE) | ",
    "Mean IQ-Big Five (SE) | IQ-Big Five minus Big Five-Big Five (SE) | ",
    "p for difference | Inference |"
  ),
  "|---|---|---|---:|---:|---:|---:|---|",
  vapply(
    seq_len(nrow(category_summary)),
    function(index) {
      row <- category_summary[index, ]
      paste0(
        "| ", row$analysis,
        " | ", row$structure,
        " | ", ifelse(row$loading_model == "equal", "Equal", "Free"),
        " | ", format_estimate_se(row$mean_big_five, row$se_big_five),
        " | ", format_estimate_se(row$mean_iq_big_five, row$se_iq_big_five),
        " | ", format_estimate_se(row$difference, row$se_difference),
        " | ", format_p(row$p_equal_means),
        " | ", ifelse(row$inference_reliable, "Reliable", "Diagnostic"),
        " |"
      )
    },
    character(1)
  )
)

jackknife_summary <- aggregate(
  cbind(
    converged,
    post_check,
    gradient_check_passed,
    terminal_covariance_admissible
  ) ~
    analysis + structure + loading_model,
  data = jackknife_diagnostics,
  FUN = mean
)
jackknife_lines <- c(
  paste0(
    "| Analysis | Structure | Loadings | Converged | Post-check | Gradient | ",
    "Terminal + covariance |"
  ),
  "|---|---|---|---:|---:|---:|---:|",
  vapply(
    seq_len(nrow(jackknife_summary)),
    function(index) {
      row <- jackknife_summary[index, ]
      paste0(
        "| ", row$analysis,
        " | ", row$structure,
        " | ", ifelse(row$loading_model == "equal", "Equal", "Free"),
        " | ", format_number(row$converged),
        " | ", format_number(row$post_check),
        " | ", format_number(row$gradient_check_passed),
        " | ",
        format_number(row$terminal_covariance_admissible),
        " |"
      )
    },
    character(1)
  )
)

report_lines <- c(
  "# SEM sensitivity analysis with freely estimated paired loadings",
  "",
  paste0("Generated by `scripts/free_split_half_loadings.R` on ", Sys.Date(), "."),
  "",
  "## Specification",
  "",
  paste0(
    "All 12 trait-to-indicator loadings are estimated separately. The two ",
    "orthogonal Big Five sample/method factors, the correlated-traits and ",
    "hierarchical structures, the full jackknife WLS precision matrix, and the ",
    "200 fixed-weight delete-one-block refits are retained from the canonical ",
    "analysis. Equal-loading models are fitted as nested references and as ",
    "like-for-like benchmarks for the mean latent correlations."
  ),
  "",
  "The constrained parameterization keeps observed residual variances positive and the freely correlated six-trait matrix positive definite.",
  "",
  paste0(
    "A preliminary bare covariance-model fit can return a lower objective by ",
    "leaving the admissible covariance space—for example, by implying a ",
    "negative residual variance or a latent correlation matrix that is not ",
    "positive definite. Such a solution is an optimizer artifact, not an ",
    "interpretable SEM. The estimates in this report do not use that solution. ",
    "The remaining gene-set warning is different: even inside the admissible ",
    "space, the free-loading parameter information is rank deficient, so the ",
    "data do not distinguish all freed parameters locally."
  ),
  "",
  "## Free-loading model fit",
  "",
  fit_lines,
  "",
  paste0(
    "RMSEA uses 200 as the repository's WLS scaling convention, not as a ",
    "literal independent-observation count. AIC and BIC are omitted because ",
    "the supplied-jackknife WLS estimator has no likelihood."
  ),
  "",
  paste0(
    "Models marked inadmissible have rank-deficient parameter information or ",
    "another failed numerical check. Their point estimates remain diagnostic, ",
    "but their conventional tests, confidence intervals, and nested chi-square ",
    "reference distributions should not be treated as reliable inference."
  ),
  "",
  "## Test of the six paired-loading equality restrictions",
  "",
  comparison_lines,
  "",
  paste0(
    "The gene-set comparisons require particular caution whenever either ",
    "model has rank-deficient information or lies near a covariance boundary."
  ),
  "",
  "## Freely estimated standardized trait loadings",
  "",
  loading_lines,
  "",
  "Full loading and ratio jackknife standard errors and confidence intervals are in `free_loading_estimates.csv`.",
  "",
  "## Sample-specific factor loadings",
  "",
  paste0(
    "Sample-specific factor 1 loads the five half-1 personality indicators; ",
    "sample-specific factor 2 loads the five half-2 indicators. IQ does not ",
    "load on either sample-specific factor. Parentheses contain 200-block ",
    "jackknife standard errors. Each factor is oriented so the sum of its five ",
    "loadings is nonnegative."
  ),
  "",
  sample_factor_lines,
  "",
  paste0(
    "The gene-set free-loading models have rank-deficient parameter ",
    "information, so their loadings and jackknife intervals are diagnostic; ",
    "z tests and p-values are omitted for those models. The complete long-form ",
    "loading table, including general-factor loadings for hierarchical models, ",
    "is in `free_loading_factor_loadings.csv`."
  ),
  "",
  "## IQ-Big Five latent correlations in the free-loading models",
  "",
  iq_lines,
  "",
  "## Mean latent-correlation comparison",
  "",
  paste0(
    "For each loading specification, the table averages the 10 correlations ",
    "among the Big Five and the 5 correlations between IQ and the Big Five. ",
    "The requested difference is `mean(IQ-Big Five) - ",
    "mean(Big Five-Big Five)`. Parentheses contain 200-block jackknife standard ",
    "errors; the p-value tests whether that difference is zero."
  ),
  "",
  category_lines,
  "",
  "## Jackknife numerical stability rates",
  "",
  jackknife_lines,
  "",
  paste0(
    "These rates summarize optimizer, gradient, and covariance-space checks ",
    "for the delete-one refits. They do not override a rank-deficient ",
    "full-sample information matrix."
  ),
  "",
  "## Output files",
  "",
  "- `free_loading_model_fit.csv`",
  "- `free_loading_model_comparisons.csv`",
  "- `free_loading_estimates.csv`",
  "- `free_loading_factor_loadings.csv`",
  "- `free_loading_trait_correlations.csv`",
  "- `free_loading_correlation_means.csv`",
  "- `free_loading_jackknife_diagnostics.csv`",
  "- `free_split_half_loading_report.md`"
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Wrote free-loading sensitivity outputs to ", output_directory)
