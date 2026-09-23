#!/usr/bin/env Rscript

# Simulation study for the split-half measurement design.
#
# Data are generated with no sample/method factors. A larger Sample 2 mean
# chi-square changes only the trait loading through
#
#   lambda_h = sqrt((mean_chisq_h - 1) / mean_chisq_h).
#
# The script then fits the original equal-loading model with sample factors and
# three diagnostic comparators. It first fits exact population correlation
# vectors, then performs targeted Monte Carlo simulations using the repository's
# jackknife covariance matrices as realistic summary-statistic noise templates.

suppressPackageStartupMessages({
  library(Matrix)
  library(numDeriv)
  library(ggplot2)
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
    "scripts/split_half_design_simulation.R",
    mustWork = TRUE
  )
}

project_directory <- dirname(dirname(script_path))
arguments <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(arguments) >= 1L) {
  path.expand(arguments[[1L]])
} else {
  file.path(project_directory, "output", "split_half_design_simulation")
}
if (!grepl("^/", output_directory)) {
  output_directory <- file.path(project_directory, output_directory)
}
monte_carlo_replications <- if (length(arguments) >= 2L) {
  as.integer(arguments[[2L]])
} else {
  100L
}
if (is.na(monte_carlo_replications) || monte_carlo_replications < 1L) {
  stop("Monte Carlo replications must be a positive integer.", call. = FALSE)
}
requested_cores <- if (length(arguments) >= 3L) {
  as.integer(arguments[[3L]])
} else {
  min(4L, max(1L, parallel::detectCores() - 1L))
}
if (is.na(requested_cores) || requested_cores < 1L) {
  stop("The core count must be a positive integer.", call. = FALSE)
}
simulation_cores <- if (.Platform$OS.type == "windows") 1L else requested_cores
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

population_path <- file.path(output_directory, "population_results.csv")
monte_carlo_replicates_path <- file.path(
  output_directory,
  "monte_carlo_replicates.csv"
)
monte_carlo_summary_path <- file.path(
  output_directory,
  "monte_carlo_summary.csv"
)
truth_path <- file.path(output_directory, "truth_scenarios.csv")
grid_path <- file.path(output_directory, "chi_square_grid.csv")
population_figure_path <- file.path(
  output_directory,
  "population_category_bias.png"
)
monte_carlo_figure_path <- file.path(
  output_directory,
  "monte_carlo_gap_bias.png"
)
report_path <- file.path(output_directory, "split_half_design_simulation_report.md")

trait_names <- c(
  "Agreeableness",
  "Conscientiousness",
  "Extraversion",
  "Neuroticism",
  "Openness",
  "IQ"
)
trait_ids <- c("agree", "consc", "extra", "neurot", "open", "iq")
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
  "iq1", "iq2"
)
indicator_trait_index <- rep(seq_along(trait_ids), each = 2L)
indicator_half <- rep(1:2, length(trait_ids))
pair_indices <- t(combn(seq_along(short_names), 2L))
pair_names <- paste(
  short_names[pair_indices[, 1L]],
  short_names[pair_indices[, 2L]],
  sep = "__"
)
big_five_pair_indices <- t(combn(seq_len(5L), 2L))

mean_chisq_1 <- mean(c(1.120713, 1.172888, 1.201807, 1.169212, 1.178110))
observed_mean_chisq_2 <- mean(c(
  1.181906,
  1.323071,
  1.343517,
  1.340171,
  1.362600
))
population_chisq_2_grid <- c(
  mean_chisq_1,
  1.20,
  1.25,
  1.30,
  observed_mean_chisq_2,
  1.35,
  1.40,
  1.50,
  1.75,
  2.00
)
population_chisq_2_grid <- sort(unique(population_chisq_2_grid))
monte_carlo_chisq_2_grid <- c(
  mean_chisq_1,
  observed_mean_chisq_2,
  1.50
)

chisq_signal_fraction <- function(mean_chisq) {
  (mean_chisq - 1) / mean_chisq
}
loading_from_chisq <- function(mean_chisq) {
  sqrt(chisq_signal_fraction(mean_chisq))
}

chi_square_grid <- data.frame(
  mean_chisq_1 = mean_chisq_1,
  mean_chisq_2 = population_chisq_2_grid,
  lambda_1 = loading_from_chisq(mean_chisq_1),
  lambda_2 = loading_from_chisq(population_chisq_2_grid),
  kappa_2_over_1 =
    loading_from_chisq(population_chisq_2_grid) /
      loading_from_chisq(mean_chisq_1),
  used_in_monte_carlo = population_chisq_2_grid %in%
    monte_carlo_chisq_2_grid,
  stringsAsFactors = FALSE
)
write.csv(chi_square_grid, grid_path, row.names = FALSE, quote = TRUE)

design_specifications <- list(
  repository_faithful = list(
    label = "Repository-faithful",
    affected_traits = seq_len(5L),
    method_indicator_indices = seq_len(10L),
    description = paste0(
      "The Big Five receive the Sample 2 chi-square increase; IQ remains ",
      "equal and is excluded from the two fitted sample factors."
    )
  ),
  all_six_traits = list(
    label = "All six traits grouped",
    affected_traits = seq_len(6L),
    method_indicator_indices = seq_len(12L),
    description = paste0(
      "All six traits receive the Sample 2 chi-square increase and all six ",
      "pairs are included in the two fitted sample factors."
    )
  )
)

model_specifications <- list(
  original_equal_sample = list(
    label = "Original: equal + sample factors",
    ratio_kind = "equal",
    include_method = TRUE
  ),
  equal_no_sample = list(
    label = "Equal, no sample factors",
    ratio_kind = "equal",
    include_method = FALSE
  ),
  estimated_common_sample = list(
    label = "Estimated common ratio + sample factors",
    ratio_kind = "common",
    include_method = TRUE
  ),
  fixed_chisq_sample = list(
    label = "Oracle fixed chi-square + sample factors",
    ratio_kind = "fixed",
    include_method = TRUE
  )
)

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

load_sampling_vcov <- function(path, description) {
  source_vcov <- read_named_matrix(path, description)
  source_pair_names <- vapply(
    seq_len(nrow(pair_indices)),
    function(index) {
      resolve_pair_name(
        long_names[pair_indices[index, 1L]],
        long_names[pair_indices[index, 2L]],
        rownames(source_vcov)
      )
    },
    character(1)
  )
  if (anyNA(source_pair_names) ||
      !all(source_pair_names %in% colnames(source_vcov))) {
    stop(description, " lacks one or more expected pairs.", call. = FALSE)
  }
  value <- source_vcov[source_pair_names, source_pair_names, drop = FALSE]
  dimnames(value) <- list(pair_names, pair_names)
  eigenvalues <- eigen(value, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= 0) {
    stop(description, " is not positive definite.", call. = FALSE)
  }
  value
}

sampling_vcov_templates <- list(
  gene_score = load_sampling_vcov(
    file.path(
      project_directory,
      "output",
      "sampling_covariances",
      "gene_score_correlation_sampling_vcov.csv"
    ),
    "Gene-score sampling covariance"
  ),
  gene_set = load_sampling_vcov(
    file.path(
      project_directory,
      "output",
      "sampling_covariances",
      "all_gene_sets_correlation_sampling_vcov.csv"
    ),
    "Gene-set sampling covariance"
  )
)

build_empirical_trait_correlation <- function(analysis_label) {
  path <- file.path(
    project_directory,
    "output",
    "latent_factor_correlations",
    "implied_trait_correlations.csv"
  )
  table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  rows <- table[
    table$analysis == analysis_label &
      table$model == "Correlated-factors model",
    ,
    drop = FALSE
  ]
  if (nrow(rows) != choose(length(trait_names), 2L)) {
    stop(
      "Could not recover the empirical trait-correlation matrix for ",
      analysis_label,
      ".",
      call. = FALSE
    )
  }
  value <- diag(length(trait_names))
  dimnames(value) <- list(trait_names, trait_names)
  for (index in seq_len(nrow(rows))) {
    left <- rows$trait_1[[index]]
    right <- rows$trait_2[[index]]
    value[left, right] <- value[right, left] <- rows$correlation[[index]]
  }
  minimum_eigenvalue <- min(eigen(
    value,
    symmetric = TRUE,
    only.values = TRUE
  )$values)
  if (minimum_eigenvalue <= 1e-8) {
    value <- as.matrix(Matrix::nearPD(
      value,
      corr = TRUE,
      keepDiag = TRUE,
      eig.tol = 1e-8,
      posd.tol = 1e-8
    )$mat)
  }
  dimnames(value) <- list(trait_ids, trait_ids)
  value
}

equicorrelation_matrix <- function(correlation) {
  value <- matrix(
    correlation,
    nrow = length(trait_ids),
    ncol = length(trait_ids),
    dimnames = list(trait_ids, trait_ids)
  )
  diag(value) <- 1
  value
}

truth_scenarios <- list(
  equal_categories_gene_score = list(
    label = "Equal categories, gene-score strength/noise",
    description = paste0(
      "All 15 trait correlations equal .40; gene-score jackknife noise."
    ),
    trait_correlation = equicorrelation_matrix(0.40),
    noise_template = "gene_score"
  ),
  equal_categories_gene_set = list(
    label = "Equal categories, gene-set strength/noise",
    description = paste0(
      "All 15 trait correlations equal .65; gene-set jackknife noise."
    ),
    trait_correlation = equicorrelation_matrix(0.65),
    noise_template = "gene_set"
  ),
  empirical_gene_score = list(
    label = "Empirical gene-score correlations",
    description = paste0(
      "Canonical correlated-traits gene-score correlation matrix and ",
      "gene-score jackknife noise."
    ),
    trait_correlation = build_empirical_trait_correlation("Gene scores"),
    noise_template = "gene_score"
  ),
  empirical_gene_set = list(
    label = "Empirical gene-set correlations",
    description = paste0(
      "Canonical correlated-traits pooled gene-set correlation matrix and ",
      "gene-set jackknife noise."
    ),
    trait_correlation = build_empirical_trait_correlation(
      "Pooled gene-set enrichment"
    ),
    noise_template = "gene_set"
  )
)

summarize_trait_correlation <- function(value) {
  mean_big_five <- mean(value[big_five_pair_indices])
  mean_iq <- mean(value[seq_len(5L), 6L])
  c(
    mean_big_five_big_five = mean_big_five,
    mean_iq_big_five = mean_iq,
    difference_iq_minus_big_five = mean_iq - mean_big_five
  )
}

truth_table <- do.call(rbind, lapply(names(truth_scenarios), function(id) {
  scenario <- truth_scenarios[[id]]
  summary <- summarize_trait_correlation(scenario$trait_correlation)
  data.frame(
    truth_scenario = id,
    truth_label = scenario$label,
    description = scenario$description,
    noise_template = scenario$noise_template,
    true_mean_big_five_big_five = summary[[1L]],
    true_mean_iq_big_five = summary[[2L]],
    true_difference_iq_minus_big_five = summary[[3L]],
    minimum_trait_correlation_eigenvalue = min(eigen(
      scenario$trait_correlation,
      symmetric = TRUE,
      only.values = TRUE
    )$values),
    stringsAsFactors = FALSE
  )
}))
rownames(truth_table) <- NULL
write.csv(truth_table, truth_path, row.names = FALSE, quote = TRUE)

correlation_matrix_from_parameters <- function(parameters, dimension) {
  lower <- diag(dimension)
  lower[lower.tri(lower)] <- parameters
  covariance <- tcrossprod(lower)
  covariance / sqrt(outer(diag(covariance), diag(covariance)))
}

matrix_to_correlation_parameters <- function(value) {
  positive_definite <- as.matrix(Matrix::nearPD(
    value,
    corr = TRUE,
    keepDiag = TRUE,
    eig.tol = 1e-8,
    posd.tol = 1e-8
  )$mat)
  lower <- t(chol(positive_definite))
  unit_lower <- sweep(lower, 1L, diag(lower), "/")
  unit_lower[lower.tri(unit_lower)]
}

generate_population_moments <- function(
  trait_correlation,
  design,
  mean_chisq_2
) {
  lambda_1 <- rep(loading_from_chisq(mean_chisq_1), length(trait_ids))
  lambda_2 <- lambda_1
  lambda_2[design$affected_traits] <- loading_from_chisq(mean_chisq_2)
  trait_loadings <- as.vector(rbind(lambda_1, lambda_2))
  measurement <- matrix(
    0,
    nrow = length(short_names),
    ncol = length(trait_ids)
  )
  measurement[cbind(
    seq_along(short_names),
    indicator_trait_index
  )] <- trait_loadings
  implied <- measurement %*% trait_correlation %*% t(measurement)
  diag(implied) <- 1
  dimnames(implied) <- list(short_names, short_names)
  list(
    correlation = implied,
    pair_vector = implied[pair_indices],
    lambda_1 = lambda_1,
    lambda_2 = lambda_2,
    ratios = lambda_2 / lambda_1
  )
}

parameter_layout <- function(model, design) {
  method_count <- if (model$include_method) {
    length(design$method_indicator_indices)
  } else {
    0L
  }
  ratio_count <- if (model$ratio_kind == "common") 1L else 0L
  trait_count <- length(trait_ids)
  correlation_count <- choose(trait_count, 2L)
  trait_indices <- seq_len(trait_count)
  method_indices <- if (method_count > 0L) {
    max(trait_indices) + seq_len(method_count)
  } else {
    integer()
  }
  correlation_start <- trait_count + method_count
  correlation_indices <- correlation_start + seq_len(correlation_count)
  ratio_indices <- if (ratio_count > 0L) {
    max(correlation_indices) + seq_len(ratio_count)
  } else {
    integer()
  }
  list(
    trait = trait_indices,
    method = method_indices,
    correlation = correlation_indices,
    ratio = ratio_indices,
    count = trait_count + method_count + correlation_count + ratio_count
  )
}

decode_model_parameters <- function(
  parameters,
  model,
  design,
  fixed_ratios,
  layout
) {
  if (model$ratio_kind == "equal") {
    ratios <- rep(1, length(trait_ids))
  } else if (model$ratio_kind == "fixed") {
    ratios <- fixed_ratios
  } else {
    ratios <- rep(1, length(trait_ids))
    ratios[design$affected_traits] <- exp(parameters[layout$ratio])
  }
  first_loading_cap <- 0.999 / pmax(1, ratios)
  first_loadings <- first_loading_cap * tanh(parameters[layout$trait])
  second_loadings <- ratios * first_loadings
  trait_loadings <- as.vector(rbind(first_loadings, second_loadings))

  trait_correlation <- correlation_matrix_from_parameters(
    parameters[layout$correlation],
    length(trait_ids)
  )
  dimnames(trait_correlation) <- list(trait_ids, trait_ids)

  measurement <- matrix(
    0,
    nrow = length(short_names),
    ncol = length(trait_ids)
  )
  measurement[cbind(
    seq_along(short_names),
    indicator_trait_index
  )] <- trait_loadings

  method_matrix <- matrix(0, nrow = length(short_names), ncol = 2L)
  method_loadings <- numeric()
  if (length(layout$method) > 0L) {
    method_indices <- design$method_indicator_indices
    method_cap <- 0.999 * sqrt(pmax(
      1e-12,
      1 - trait_loadings[method_indices]^2
    ))
    method_loadings <- method_cap * tanh(parameters[layout$method])
    method_matrix[cbind(
      method_indices,
      indicator_half[method_indices]
    )] <- method_loadings
  }

  implied <-
    measurement %*% trait_correlation %*% t(measurement) +
      tcrossprod(method_matrix)
  residual_variances <- 1 - diag(implied)
  diag(implied) <- 1
  dimnames(implied) <- list(short_names, short_names)

  list(
    ratios = ratios,
    first_loadings = first_loadings,
    second_loadings = second_loadings,
    trait_loadings = trait_loadings,
    method_loadings = method_loadings,
    method_matrix = method_matrix,
    trait_correlation = trait_correlation,
    residual_variances = residual_variances,
    implied_correlation = implied
  )
}

build_start <- function(
  observed_vector,
  model,
  design,
  fixed_ratios,
  layout,
  start_number,
  seed
) {
  observed <- diag(length(short_names))
  observed[pair_indices] <- observed_vector
  observed[pair_indices[, 2:1, drop = FALSE]] <- observed_vector
  ratios <- if (model$ratio_kind == "fixed") {
    fixed_ratios
  } else if (model$ratio_kind == "common") {
    common_value <- exp(mean(log(fixed_ratios[design$affected_traits])))
    value <- rep(1, length(trait_ids))
    value[design$affected_traits] <- common_value
    value
  } else {
    rep(1, length(trait_ids))
  }
  within_trait <- observed[cbind(
    seq(1L, 12L, by = 2L),
    seq(2L, 12L, by = 2L)
  )]
  first_loading <- sqrt(pmax(0.01, abs(within_trait)) / ratios)
  first_cap <- 0.999 / pmax(1, ratios)
  first_loading <- pmin(0.90 * first_cap, first_loading)
  trait_raw <- atanh(pmin(0.95, first_loading / first_cap))
  second_loading <- ratios * first_loading

  estimated_trait_correlation <- diag(length(trait_ids))
  for (left in seq_len(length(trait_ids) - 1L)) {
    for (right in seq.int(left + 1L, length(trait_ids))) {
      left_indicators <- c(2L * left - 1L, 2L * left)
      right_indicators <- c(2L * right - 1L, 2L * right)
      numerator <- c(
        observed[left_indicators[[1L]], right_indicators[[1L]]],
        observed[left_indicators[[1L]], right_indicators[[2L]]],
        observed[left_indicators[[2L]], right_indicators[[1L]]],
        observed[left_indicators[[2L]], right_indicators[[2L]]]
      )
      denominator <- c(
        first_loading[[left]] * first_loading[[right]],
        first_loading[[left]] * second_loading[[right]],
        second_loading[[left]] * first_loading[[right]],
        second_loading[[left]] * second_loading[[right]]
      )
      estimate <- mean(numerator / denominator)
      estimate <- pmax(-0.90, pmin(0.90, estimate))
      estimated_trait_correlation[left, right] <-
        estimated_trait_correlation[right, left] <- estimate
    }
  }
  correlation_raw <- matrix_to_correlation_parameters(
    estimated_trait_correlation
  )

  method_raw <- numeric()
  if (length(layout$method) > 0L) {
    set.seed(seed + start_number * 1009L)
    if (start_number == 1L) {
      method_raw <- rep(0.05, length(layout$method))
    } else if (start_number == 2L) {
      method_raw <- rep(c(0.08, -0.04), length.out = length(layout$method))
    } else {
      method_raw <- rnorm(length(layout$method), sd = 0.10)
    }
  }
  ratio_raw <- if (model$ratio_kind == "common") {
    log(exp(mean(log(fixed_ratios[design$affected_traits]))))
  } else {
    numeric()
  }
  c(trait_raw, method_raw, correlation_raw, ratio_raw)
}

fit_simulation_model <- function(
  observed_vector,
  sampling_vcov,
  model_id,
  design_id,
  fixed_ratios,
  seed,
  starts = 2L,
  calculate_information = FALSE
) {
  model <- model_specifications[[model_id]]
  design <- design_specifications[[design_id]]
  layout <- parameter_layout(model, design)
  precision <- solve(sampling_vcov)
  objective <- function(parameters) {
    implied <- decode_model_parameters(
      parameters,
      model,
      design,
      fixed_ratios,
      layout
    )$implied_correlation[pair_indices]
    residual <- observed_vector - implied
    as.numeric(crossprod(residual, precision %*% residual))
  }
  lower <- rep(-8, layout$count)
  upper <- rep(8, layout$count)
  if (length(layout$ratio) > 0L) {
    lower[layout$ratio] <- log(0.30)
    upper[layout$ratio] <- log(3.00)
  }
  optimizers <- lapply(seq_len(starts), function(start_number) {
    start <- build_start(
      observed_vector,
      model,
      design,
      fixed_ratios,
      layout,
      start_number,
      seed
    )
    tryCatch(
      nlminb(
        start = pmax(lower, pmin(upper, start)),
        objective = objective,
        lower = lower,
        upper = upper,
        control = list(
          iter.max = if (calculate_information) 5000L else 1500L,
          eval.max = if (calculate_information) 12000L else 4000L,
          rel.tol = 1e-9
        )
      ),
      error = function(condition) NULL
    )
  })
  valid <- vapply(
    optimizers,
    function(value) !is.null(value) && is.finite(value$objective),
    logical(1)
  )
  if (!any(valid)) {
    return(list(success = FALSE, parameter_count = layout$count))
  }
  valid_optimizers <- optimizers[valid]
  objectives <- vapply(valid_optimizers, `[[`, numeric(1), "objective")
  optimizer <- valid_optimizers[[which.min(objectives)]]
  decoded <- decode_model_parameters(
    optimizer$par,
    model,
    design,
    fixed_ratios,
    layout
  )
  orientation <- sign(decoded$first_loadings + decoded$second_loadings)
  orientation[orientation == 0] <- 1
  oriented_trait_correlation <- decoded$trait_correlation *
    outer(orientation, orientation)
  category_summary <- summarize_trait_correlation(
    oriented_trait_correlation
  )
  affected_ratio <- mean(decoded$ratios[design$affected_traits])
  method_variance <- if (length(decoded$method_loadings) > 0L) {
    mean(decoded$method_loadings^2)
  } else {
    0
  }
  method_covariance <- tcrossprod(decoded$method_matrix)
  method_pair_values <- method_covariance[pair_indices]
  same_half <- indicator_half[pair_indices[, 1L]] ==
    indicator_half[pair_indices[, 2L]]
  included_method_pair <-
    pair_indices[, 1L] %in% design$method_indicator_indices &
      pair_indices[, 2L] %in% design$method_indicator_indices &
      same_half
  method_covariance_rms <- if (any(included_method_pair)) {
    sqrt(mean(method_pair_values[included_method_pair]^2))
  } else {
    0
  }
  residual <- observed_vector -
    decoded$implied_correlation[pair_indices]
  degrees_of_freedom <- length(observed_vector) - layout$count
  minimum_trait_eigenvalue <- min(eigen(
    oriented_trait_correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values)
  terminal <- optimizer$convergence %in% c(0L, 1L)
  admissible <- terminal &&
    min(decoded$residual_variances) > 0 &&
    minimum_trait_eigenvalue > 0

  information_rank <- NA_integer_
  information_dimension <- layout$count
  if (calculate_information) {
    moment_jacobian <- numDeriv::jacobian(
      function(parameters) {
        decode_model_parameters(
          parameters,
          model,
          design,
          fixed_ratios,
          layout
        )$implied_correlation[pair_indices]
      },
      optimizer$par
    )
    information <- crossprod(moment_jacobian, precision %*% moment_jacobian)
    eigenvalues <- eigen(
      information,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    tolerance <- max(eigenvalues) * max(dim(information)) *
      .Machine$double.eps^0.75
    information_rank <- sum(eigenvalues > tolerance)
  }

  list(
    success = TRUE,
    objective = optimizer$objective,
    degrees_of_freedom = degrees_of_freedom,
    fit_p_value = pchisq(
      optimizer$objective,
      degrees_of_freedom,
      lower.tail = FALSE
    ),
    off_diagonal_rmsr = sqrt(mean(residual^2)),
    mean_big_five_big_five = category_summary[[1L]],
    mean_iq_big_five = category_summary[[2L]],
    difference_iq_minus_big_five = category_summary[[3L]],
    estimated_affected_ratio = affected_ratio,
    mean_method_variance = method_variance,
    method_covariance_rms = method_covariance_rms,
    minimum_residual_variance = min(decoded$residual_variances),
    minimum_trait_correlation_eigenvalue = minimum_trait_eigenvalue,
    optimizer_convergence_code = optimizer$convergence,
    admissible = admissible,
    information_rank = information_rank,
    information_dimension = information_dimension,
    parameter_count = layout$count
  )
}

result_row_from_fit <- function(
  fit,
  design_id,
  truth_id,
  mean_chisq_2,
  model_id,
  replication = NA_integer_
) {
  truth <- truth_scenarios[[truth_id]]
  true_summary <- summarize_trait_correlation(truth$trait_correlation)
  generated <- generate_population_moments(
    truth$trait_correlation,
    design_specifications[[design_id]],
    mean_chisq_2
  )
  true_ratio <- mean(
    generated$ratios[design_specifications[[design_id]]$affected_traits]
  )
  empty <- !isTRUE(fit$success)
  extract <- function(name, default = NA_real_) {
    if (empty || is.null(fit[[name]])) default else fit[[name]]
  }
  data.frame(
    replication = replication,
    design = design_id,
    design_label = design_specifications[[design_id]]$label,
    truth_scenario = truth_id,
    truth_label = truth$label,
    noise_template = truth$noise_template,
    mean_chisq_1 = mean_chisq_1,
    mean_chisq_2 = mean_chisq_2,
    true_affected_ratio = true_ratio,
    model = model_id,
    model_label = model_specifications[[model_id]]$label,
    true_mean_big_five_big_five = true_summary[[1L]],
    estimated_mean_big_five_big_five = extract(
      "mean_big_five_big_five"
    ),
    bias_big_five_big_five =
      extract("mean_big_five_big_five") - true_summary[[1L]],
    true_mean_iq_big_five = true_summary[[2L]],
    estimated_mean_iq_big_five = extract("mean_iq_big_five"),
    bias_iq_big_five = extract("mean_iq_big_five") - true_summary[[2L]],
    true_difference_iq_minus_big_five = true_summary[[3L]],
    estimated_difference_iq_minus_big_five = extract(
      "difference_iq_minus_big_five"
    ),
    bias_difference_iq_minus_big_five =
      extract("difference_iq_minus_big_five") - true_summary[[3L]],
    estimated_affected_ratio = extract("estimated_affected_ratio"),
    bias_affected_ratio = extract("estimated_affected_ratio") - true_ratio,
    mean_method_variance = extract("mean_method_variance"),
    method_covariance_rms = extract("method_covariance_rms"),
    chi_square = extract("objective"),
    degrees_of_freedom = extract("degrees_of_freedom"),
    fit_p_value = extract("fit_p_value"),
    off_diagonal_rmsr = extract("off_diagonal_rmsr"),
    optimizer_convergence_code = extract("optimizer_convergence_code"),
    admissible = if (empty) FALSE else isTRUE(fit$admissible),
    minimum_residual_variance = extract("minimum_residual_variance"),
    minimum_trait_correlation_eigenvalue = extract(
      "minimum_trait_correlation_eigenvalue"
    ),
    information_rank = extract("information_rank", NA_integer_),
    information_dimension = extract("information_dimension", NA_integer_),
    parameter_count = extract("parameter_count", NA_integer_),
    stringsAsFactors = FALSE
  )
}

message("Running exact-population simulation grid...")
population_rows <- list()
population_index <- 0L
for (design_id in names(design_specifications)) {
  for (truth_id in names(truth_scenarios)) {
    truth <- truth_scenarios[[truth_id]]
    sampling_vcov <- sampling_vcov_templates[[truth$noise_template]]
    for (chisq_index in seq_along(population_chisq_2_grid)) {
      mean_chisq_2 <- population_chisq_2_grid[[chisq_index]]
      generated <- generate_population_moments(
        truth$trait_correlation,
        design_specifications[[design_id]],
        mean_chisq_2
      )
      for (model_index in seq_along(model_specifications)) {
        model_id <- names(model_specifications)[[model_index]]
        seed <- 100000L +
          match(design_id, names(design_specifications)) * 10000L +
          match(truth_id, names(truth_scenarios)) * 1000L +
          chisq_index * 10L + model_index
        fit <- fit_simulation_model(
          generated$pair_vector,
          sampling_vcov,
          model_id,
          design_id,
          generated$ratios,
          seed = seed,
          starts = 3L,
          calculate_information = TRUE
        )
        population_index <- population_index + 1L
        population_rows[[population_index]] <- result_row_from_fit(
          fit,
          design_id,
          truth_id,
          mean_chisq_2,
          model_id
        )
      }
    }
    message(
      "  completed ",
      design_specifications[[design_id]]$label,
      " | ",
      truth$label
    )
  }
}
population_results <- do.call(rbind, population_rows)
rownames(population_results) <- NULL
write.csv(
  population_results,
  population_path,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

monte_carlo_conditions <- expand.grid(
  design = names(design_specifications),
  truth_scenario = names(truth_scenarios),
  mean_chisq_2 = monte_carlo_chisq_2_grid,
  replication = seq_len(monte_carlo_replications),
  stringsAsFactors = FALSE
)
monte_carlo_conditions$condition_index <- seq_len(nrow(monte_carlo_conditions))

run_monte_carlo_condition <- function(condition_index) {
  condition <- monte_carlo_conditions[condition_index, , drop = FALSE]
  design_id <- condition$design[[1L]]
  truth_id <- condition$truth_scenario[[1L]]
  mean_chisq_2 <- condition$mean_chisq_2[[1L]]
  replication <- condition$replication[[1L]]
  truth <- truth_scenarios[[truth_id]]
  sampling_vcov <- sampling_vcov_templates[[truth$noise_template]]
  generated <- generate_population_moments(
    truth$trait_correlation,
    design_specifications[[design_id]],
    mean_chisq_2
  )
  seed <- 5000000L +
    match(design_id, names(design_specifications)) * 1000000L +
    match(truth_id, names(truth_scenarios)) * 100000L +
    match(mean_chisq_2, monte_carlo_chisq_2_grid) * 10000L +
    replication
  set.seed(seed)
  noise <- as.numeric(
    t(chol(sampling_vcov)) %*% rnorm(nrow(sampling_vcov))
  )
  observed_vector <- generated$pair_vector + noise
  rows <- lapply(seq_along(model_specifications), function(model_index) {
    model_id <- names(model_specifications)[[model_index]]
    fit <- fit_simulation_model(
      observed_vector,
      sampling_vcov,
      model_id,
      design_id,
      generated$ratios,
      seed = seed + model_index * 101L,
      starts = 2L,
      calculate_information = FALSE
    )
    result_row_from_fit(
      fit,
      design_id,
      truth_id,
      mean_chisq_2,
      model_id,
      replication
    )
  })
  do.call(rbind, rows)
}

message(
  "Running ",
  nrow(monte_carlo_conditions),
  " Monte Carlo datasets (",
  monte_carlo_replications,
  " replications per condition; ",
  simulation_cores,
  " core(s))..."
)
if (simulation_cores > 1L) {
  monte_carlo_rows <- parallel::mclapply(
    seq_len(nrow(monte_carlo_conditions)),
    run_monte_carlo_condition,
    mc.cores = simulation_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  monte_carlo_rows <- lapply(
    seq_len(nrow(monte_carlo_conditions)),
    run_monte_carlo_condition
  )
}
monte_carlo_replicates <- do.call(rbind, monte_carlo_rows)
rownames(monte_carlo_replicates) <- NULL
write.csv(
  monte_carlo_replicates,
  monte_carlo_replicates_path,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

safe_mean <- function(value) {
  if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE)
}
safe_sd <- function(value) {
  if (sum(!is.na(value)) < 2L) NA_real_ else sd(value, na.rm = TRUE)
}
safe_rmse <- function(value) {
  if (all(is.na(value))) NA_real_ else sqrt(mean(value^2, na.rm = TRUE))
}
safe_quantile <- function(value, probability) {
  if (all(is.na(value))) {
    NA_real_
  } else {
    unname(quantile(value, probability, na.rm = TRUE, type = 8L))
  }
}

summary_groups <- interaction(
  monte_carlo_replicates$design,
  monte_carlo_replicates$truth_scenario,
  monte_carlo_replicates$mean_chisq_2,
  monte_carlo_replicates$model,
  drop = TRUE
)
monte_carlo_summary_rows <- lapply(split(
  monte_carlo_replicates,
  summary_groups
), function(rows) {
  data.frame(
    design = rows$design[[1L]],
    design_label = rows$design_label[[1L]],
    truth_scenario = rows$truth_scenario[[1L]],
    truth_label = rows$truth_label[[1L]],
    noise_template = rows$noise_template[[1L]],
    mean_chisq_1 = rows$mean_chisq_1[[1L]],
    mean_chisq_2 = rows$mean_chisq_2[[1L]],
    true_affected_ratio = rows$true_affected_ratio[[1L]],
    model = rows$model[[1L]],
    model_label = rows$model_label[[1L]],
    replications = nrow(rows),
    admissible_replications = sum(rows$admissible),
    admissible_rate = mean(rows$admissible),
    true_mean_big_five_big_five =
      rows$true_mean_big_five_big_five[[1L]],
    mean_estimated_big_five_big_five = safe_mean(
      rows$estimated_mean_big_five_big_five
    ),
    bias_big_five_big_five = safe_mean(rows$bias_big_five_big_five),
    rmse_big_five_big_five = safe_rmse(rows$bias_big_five_big_five),
    mc_sd_big_five_big_five = safe_sd(
      rows$estimated_mean_big_five_big_five
    ),
    true_mean_iq_big_five = rows$true_mean_iq_big_five[[1L]],
    mean_estimated_iq_big_five = safe_mean(rows$estimated_mean_iq_big_five),
    bias_iq_big_five = safe_mean(rows$bias_iq_big_five),
    rmse_iq_big_five = safe_rmse(rows$bias_iq_big_five),
    mc_sd_iq_big_five = safe_sd(rows$estimated_mean_iq_big_five),
    true_difference_iq_minus_big_five =
      rows$true_difference_iq_minus_big_five[[1L]],
    mean_estimated_difference_iq_minus_big_five = safe_mean(
      rows$estimated_difference_iq_minus_big_five
    ),
    bias_difference_iq_minus_big_five = safe_mean(
      rows$bias_difference_iq_minus_big_five
    ),
    rmse_difference_iq_minus_big_five = safe_rmse(
      rows$bias_difference_iq_minus_big_five
    ),
    mc_sd_difference_iq_minus_big_five = safe_sd(
      rows$estimated_difference_iq_minus_big_five
    ),
    q025_difference = safe_quantile(
      rows$estimated_difference_iq_minus_big_five,
      0.025
    ),
    q975_difference = safe_quantile(
      rows$estimated_difference_iq_minus_big_five,
      0.975
    ),
    mean_estimated_affected_ratio = safe_mean(rows$estimated_affected_ratio),
    bias_affected_ratio = safe_mean(rows$bias_affected_ratio),
    mean_method_variance = safe_mean(rows$mean_method_variance),
    mean_method_covariance_rms = safe_mean(rows$method_covariance_rms),
    mean_chi_square = safe_mean(rows$chi_square),
    fit_rejection_rate_05 = safe_mean(rows$fit_p_value < 0.05),
    mean_off_diagonal_rmsr = safe_mean(rows$off_diagonal_rmsr),
    stringsAsFactors = FALSE
  )
})
monte_carlo_summary <- do.call(rbind, monte_carlo_summary_rows)
rownames(monte_carlo_summary) <- NULL
monte_carlo_summary <- monte_carlo_summary[order(
  match(monte_carlo_summary$design, names(design_specifications)),
  match(monte_carlo_summary$truth_scenario, names(truth_scenarios)),
  monte_carlo_summary$mean_chisq_2,
  match(monte_carlo_summary$model, names(model_specifications))
), ]
write.csv(
  monte_carlo_summary,
  monte_carlo_summary_path,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

model_colors <- c(
  original_equal_sample = "#D55E00",
  equal_no_sample = "#999999",
  estimated_common_sample = "#0072B2",
  fixed_chisq_sample = "#009E73"
)
population_long <- rbind(
  data.frame(
    population_results[c(
      "design", "design_label", "truth_scenario", "truth_label",
      "mean_chisq_2", "model", "model_label"
    )],
    category = "Big Five-Big Five",
    bias = population_results$bias_big_five_big_five
  ),
  data.frame(
    population_results[c(
      "design", "design_label", "truth_scenario", "truth_label",
      "mean_chisq_2", "model", "model_label"
    )],
    category = "IQ-Big Five",
    bias = population_results$bias_iq_big_five
  )
)
population_plot <- ggplot(
  population_long,
  aes(
    x = mean_chisq_2,
    y = bias,
    color = model,
    linetype = category,
    group = interaction(model, category)
  )
) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
  geom_line(linewidth = 0.65) +
  facet_grid(design_label ~ truth_label, scales = "free_y") +
  scale_color_manual(
    values = model_colors,
    labels = vapply(model_specifications, `[[`, character(1), "label")
  ) +
  labs(
    x = "Sample 2 mean chi-square",
    y = "Population bias in mean latent correlation",
    color = "Fitted model",
    linetype = "Correlation category"
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )
ggsave(
  population_figure_path,
  population_plot,
  width = 15,
  height = 8.5,
  dpi = 180
)

monte_carlo_plot <- ggplot(
  monte_carlo_summary,
  aes(
    x = mean_chisq_2,
    y = bias_difference_iq_minus_big_five,
    color = model,
    group = model
  )
) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.5) +
  facet_grid(design_label ~ truth_label, scales = "free_y") +
  scale_color_manual(
    values = model_colors,
    labels = vapply(model_specifications, `[[`, character(1), "label")
  ) +
  labs(
    x = "Sample 2 mean chi-square",
    y = "Monte Carlo bias in IQ-minus-Big-Five mean difference",
    color = "Fitted model"
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 7),
    legend.text = element_text(size = 7)
  )
ggsave(
  monte_carlo_figure_path,
  monte_carlo_plot,
  width = 15,
  height = 8.5,
  dpi = 180
)

format_number <- function(value, digits = 3L) {
  ifelse(is.na(value), "NA", formatC(value, digits = digits, format = "f"))
}
format_p <- function(value) {
  ifelse(
    is.na(value),
    "NA",
    ifelse(
      value < 0.001,
      "< .001",
      sub("^0", "", formatC(value, digits = 3L, format = "f"))
    )
  )
}

grid_table_lines <- c(
  "| Sample 2 mean chi-square | Loading 1 | Loading 2 | True kappa | Monte Carlo |",
  "|---:|---:|---:|---:|:---:|",
  vapply(seq_len(nrow(chi_square_grid)), function(index) {
    row <- chi_square_grid[index, ]
    paste0(
      "| ", format_number(row$mean_chisq_2, 3L),
      " | ", format_number(row$lambda_1),
      " | ", format_number(row$lambda_2),
      " | ", format_number(row$kappa_2_over_1),
      " | ", ifelse(row$used_in_monte_carlo, "Yes", "No"),
      " |"
    )
  }, character(1))
)

truth_table_lines <- c(
  "| Truth scenario | Big Five-Big Five | IQ-Big Five | True difference | Noise template |",
  "|---|---:|---:|---:|---|",
  vapply(seq_len(nrow(truth_table)), function(index) {
    row <- truth_table[index, ]
    paste0(
      "| ", row$truth_label,
      " | ", format_number(row$true_mean_big_five_big_five),
      " | ", format_number(row$true_mean_iq_big_five),
      " | ", format_number(row$true_difference_iq_minus_big_five),
      " | ", row$noise_template,
      " |"
    )
  }, character(1))
)

observed_population <- population_results[
  abs(population_results$mean_chisq_2 - observed_mean_chisq_2) < 1e-10,
  ,
  drop = FALSE
]
population_key <- do.call(rbind, lapply(split(
  observed_population,
  interaction(observed_population$design, observed_population$model, drop = TRUE)
), function(rows) {
  data.frame(
    design = rows$design[[1L]],
    design_label = rows$design_label[[1L]],
    model = rows$model[[1L]],
    model_label = rows$model_label[[1L]],
    mean_absolute_big_five_bias = mean(abs(rows$bias_big_five_big_five)),
    mean_absolute_iq_bias = mean(abs(rows$bias_iq_big_five)),
    mean_absolute_difference_bias = mean(abs(
      rows$bias_difference_iq_minus_big_five
    )),
    maximum_absolute_difference_bias = max(abs(
      rows$bias_difference_iq_minus_big_five
    )),
    mean_method_variance = mean(rows$mean_method_variance),
    mean_rmsr = mean(rows$off_diagonal_rmsr),
    stringsAsFactors = FALSE
  )
}))
population_key <- population_key[order(
  match(population_key$design, names(design_specifications)),
  match(population_key$model, names(model_specifications))
), ]
population_key_lines <- c(
  paste0(
    "| Design | Fitted model | Mean absolute Big Five bias | Mean absolute IQ bias | ",
    "Mean absolute difference bias | Maximum difference bias | Mean method variance |"
  ),
  "|---|---|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(population_key)), function(index) {
    row <- population_key[index, ]
    paste0(
      "| ", row$design_label,
      " | ", row$model_label,
      " | ", format_number(row$mean_absolute_big_five_bias),
      " | ", format_number(row$mean_absolute_iq_bias),
      " | ", format_number(row$mean_absolute_difference_bias),
      " | ", format_number(row$maximum_absolute_difference_bias),
      " | ", format_number(row$mean_method_variance),
      " |"
    )
  }, character(1))
)

identification_key <- do.call(rbind, lapply(split(
  observed_population,
  interaction(observed_population$design, observed_population$model, drop = TRUE)
), function(rows) {
  data.frame(
    design = rows$design[[1L]],
    design_label = rows$design_label[[1L]],
    model = rows$model[[1L]],
    model_label = rows$model_label[[1L]],
    information_dimension = rows$information_dimension[[1L]],
    minimum_information_rank = min(rows$information_rank),
    maximum_information_rank = max(rows$information_rank),
    full_rank_scenarios = sum(
      rows$information_rank == rows$information_dimension
    ),
    scenarios = nrow(rows),
    stringsAsFactors = FALSE
  )
}))
identification_key <- identification_key[order(
  match(identification_key$design, names(design_specifications)),
  match(identification_key$model, names(model_specifications))
), ]
identification_key_lines <- c(
  "| Design | Fitted model | Parameters | Information rank | Full-rank scenarios |",
  "|---|---|---:|---:|---:|",
  vapply(seq_len(nrow(identification_key)), function(index) {
    row <- identification_key[index, ]
    rank_text <- if (row$minimum_information_rank == row$maximum_information_rank) {
      as.character(row$minimum_information_rank)
    } else {
      paste0(row$minimum_information_rank, "-", row$maximum_information_rank)
    }
    paste0(
      "| ", row$design_label,
      " | ", row$model_label,
      " | ", row$information_dimension,
      " | ", rank_text,
      " | ", row$full_rank_scenarios, "/", row$scenarios,
      " |"
    )
  }, character(1))
)

observed_monte_carlo <- monte_carlo_summary[
  abs(monte_carlo_summary$mean_chisq_2 - observed_mean_chisq_2) < 1e-10,
  ,
  drop = FALSE
]
monte_carlo_key <- do.call(rbind, lapply(split(
  observed_monte_carlo,
  interaction(observed_monte_carlo$design, observed_monte_carlo$model, drop = TRUE)
), function(rows) {
  data.frame(
    design = rows$design[[1L]],
    design_label = rows$design_label[[1L]],
    model = rows$model[[1L]],
    model_label = rows$model_label[[1L]],
    mean_absolute_big_five_bias = mean(abs(rows$bias_big_five_big_five)),
    mean_absolute_iq_bias = mean(abs(rows$bias_iq_big_five)),
    mean_absolute_difference_bias = mean(abs(
      rows$bias_difference_iq_minus_big_five
    )),
    mean_difference_rmse = mean(rows$rmse_difference_iq_minus_big_five),
    mean_method_variance = mean(rows$mean_method_variance),
    mean_admissible_rate = mean(rows$admissible_rate),
    mean_fit_rejection_rate = mean(rows$fit_rejection_rate_05),
    stringsAsFactors = FALSE
  )
}))
monte_carlo_key <- monte_carlo_key[order(
  match(monte_carlo_key$design, names(design_specifications)),
  match(monte_carlo_key$model, names(model_specifications))
), ]
monte_carlo_key_lines <- c(
  paste0(
    "| Design | Fitted model | Mean absolute Big Five bias | Mean absolute IQ bias | ",
    "Mean absolute difference bias | Difference RMSE | Admissible rate | Fit rejection rate |"
  ),
  "|---|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(monte_carlo_key)), function(index) {
    row <- monte_carlo_key[index, ]
    paste0(
      "| ", row$design_label,
      " | ", row$model_label,
      " | ", format_number(row$mean_absolute_big_five_bias),
      " | ", format_number(row$mean_absolute_iq_bias),
      " | ", format_number(row$mean_absolute_difference_bias),
      " | ", format_number(row$mean_difference_rmse),
      " | ", format_number(row$mean_admissible_rate),
      " | ", format_number(row$mean_fit_rejection_rate),
      " |"
    )
  }, character(1))
)

stress_monte_carlo <- monte_carlo_summary[
  abs(monte_carlo_summary$mean_chisq_2 - 1.50) < 1e-10,
  ,
  drop = FALSE
]
stress_key <- do.call(rbind, lapply(split(
  stress_monte_carlo,
  interaction(stress_monte_carlo$design, stress_monte_carlo$model, drop = TRUE)
), function(rows) {
  data.frame(
    design = rows$design[[1L]],
    design_label = rows$design_label[[1L]],
    model = rows$model[[1L]],
    model_label = rows$model_label[[1L]],
    mean_absolute_big_five_bias = mean(abs(rows$bias_big_five_big_five)),
    mean_absolute_iq_bias = mean(abs(rows$bias_iq_big_five)),
    mean_absolute_difference_bias = mean(abs(
      rows$bias_difference_iq_minus_big_five
    )),
    mean_difference_rmse = mean(rows$rmse_difference_iq_minus_big_five),
    stringsAsFactors = FALSE
  )
}))
stress_key <- stress_key[order(
  match(stress_key$design, names(design_specifications)),
  match(stress_key$model, names(model_specifications))
), ]
stress_key_lines <- c(
  paste0(
    "| Design | Fitted model | Mean absolute Big Five bias | Mean absolute IQ bias | ",
    "Mean absolute difference bias | Difference RMSE |"
  ),
  "|---|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(stress_key)), function(index) {
    row <- stress_key[index, ]
    paste0(
      "| ", row$design_label,
      " | ", row$model_label,
      " | ", format_number(row$mean_absolute_big_five_bias),
      " | ", format_number(row$mean_absolute_iq_bias),
      " | ", format_number(row$mean_absolute_difference_bias),
      " | ", format_number(row$mean_difference_rmse),
      " |"
    )
  }, character(1))
)

equal_truth_population <- observed_population[
  grepl("^equal_categories", observed_population$truth_scenario),
  ,
  drop = FALSE
]
equal_truth_lines <- c(
  "| Design | Truth scenario | Fitted model | Estimated true-zero difference | Big Five bias | IQ bias |",
  "|---|---|---|---:|---:|---:|",
  vapply(seq_len(nrow(equal_truth_population)), function(index) {
    row <- equal_truth_population[index, ]
    paste0(
      "| ", row$design_label,
      " | ", row$truth_label,
      " | ", row$model_label,
      " | ", format_number(row$estimated_difference_iq_minus_big_five),
      " | ", format_number(row$bias_big_five_big_five),
      " | ", format_number(row$bias_iq_big_five),
      " |"
    )
  }, character(1))
)

empirical_truth_population <- observed_population[
  grepl("^empirical", observed_population$truth_scenario),
  ,
  drop = FALSE
]
empirical_truth_lines <- c(
  paste0(
    "| Design | Truth scenario | Fitted model | True Big Five | Estimated Big Five | ",
    "True IQ | Estimated IQ | True difference | Estimated difference |"
  ),
  "|---|---|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(empirical_truth_population)), function(index) {
    row <- empirical_truth_population[index, ]
    paste0(
      "| ", row$design_label,
      " | ", row$truth_label,
      " | ", row$model_label,
      " | ", format_number(row$true_mean_big_five_big_five),
      " | ", format_number(row$estimated_mean_big_five_big_five),
      " | ", format_number(row$true_mean_iq_big_five),
      " | ", format_number(row$estimated_mean_iq_big_five),
      " | ", format_number(row$true_difference_iq_minus_big_five),
      " | ", format_number(row$estimated_difference_iq_minus_big_five),
      " |"
    )
  }, character(1))
)

original_population_key <- population_key[
  population_key$model == "original_equal_sample",
  ,
  drop = FALSE
]
oracle_population_key <- population_key[
  population_key$model == "fixed_chisq_sample",
  ,
  drop = FALSE
]
original_mc_key <- monte_carlo_key[
  monte_carlo_key$model == "original_equal_sample",
  ,
  drop = FALSE
]
oracle_mc_key <- monte_carlo_key[
  monte_carlo_key$model == "fixed_chisq_sample",
  ,
  drop = FALSE
]

report_lines <- c(
  "# Split-half design simulation",
  "",
  paste0(
    "Generated by `scripts/split_half_design_simulation.R` on ",
    Sys.Date(),
    "."
  ),
  "",
  "## Purpose",
  "",
  paste0(
    "This simulation tests whether the original equal-loading design can ",
    "manufacture or distort the difference between mean IQ-Big Five and mean ",
    "Big Five-Big Five latent correlations when there are no true sample ",
    "factors but Sample 2 has a larger mean chi-square."
  ),
  "",
  "## Data-generating process",
  "",
  paste0(
    "Each standardized indicator is generated as `Y_th = lambda_h * eta_t + ",
    "sqrt(1 - lambda_h^2) * error_th`. All errors are mutually independent, so ",
    "both true sample factors are exactly zero. Loadings are determined by ",
    "`lambda_h = sqrt((mean_chisq_h - 1) / mean_chisq_h)`. The primary ",
    "repository-faithful design applies the imbalance to the Big Five while IQ ",
    "remains equal; the all-six design applies it to every trait."
  ),
  "",
  paste0(
    "The exact-population stage uses every chi-square value below. The Monte ",
    "Carlo stage uses the null, observed-average, and stress conditions, with ",
    monte_carlo_replications,
    " replications per design, truth scenario, and condition."
  ),
  "",
  grid_table_lines,
  "",
  "## Truth scenarios",
  "",
  truth_table_lines,
  "",
  "## Fitted models",
  "",
  "1. **Original:** equal paired loadings plus two sample factors.",
  "2. **Equal, no sample factors:** isolates the loading-equality restriction.",
  "3. **Estimated common ratio:** estimates the shared ratio for affected traits while retaining sample factors.",
  "4. **Oracle fixed chi-square:** fixes the ratio to its generating value while retaining sample factors.",
  "",
  paste0(
    "All fitted models use a freely correlated six-trait matrix. This isolates ",
    "measurement-design bias from any additional hierarchical structural ",
    "restriction. Population fits use the matching jackknife precision matrix; ",
    "Monte Carlo correlation vectors are generated by adding multivariate-normal ",
    "noise with that jackknife covariance."
  ),
  "",
  "### Local identification diagnostic",
  "",
  paste0(
    "The table reports the numerical rank of the expected-information matrix at ",
    "the observed-average population solutions across the four truth scenarios. ",
    "A rank smaller than the parameter count reflects the nonregular zero-method-",
    "factor solution; it does not prevent exact recovery of the trait correlations ",
    "in the correctly ratio-specified models."
  ),
  "",
  identification_key_lines,
  "",
  "## Exact-population results at the observed average imbalance",
  "",
  paste0(
    "The observed Big Five averages are mean chi-square 1 = ",
    format_number(mean_chisq_1),
    " and mean chi-square 2 = ",
    format_number(observed_mean_chisq_2),
    ", corresponding to kappa = ",
    format_number(
      loading_from_chisq(observed_mean_chisq_2) /
        loading_from_chisq(mean_chisq_1)
    ),
    ". The table averages absolute bias across the four truth scenarios."
  ),
  "",
  population_key_lines,
  "",
  "### Equal-category truth",
  "",
  paste0(
    "In these scenarios the true IQ-minus-Big-Five difference is exactly zero. ",
    "Any nonzero fitted value is therefore manufactured by the measurement ",
    "model."
  ),
  "",
  equal_truth_lines,
  "",
  "### Empirical-correlation truth",
  "",
  paste0(
    "These rows show whether the fitted model preserves the category means and ",
    "their difference when the generating trait correlations match the canonical ",
    "gene-score or pooled gene-set estimates."
  ),
  "",
  empirical_truth_lines,
  "",
  paste0("![Population bias](", basename(population_figure_path), ")"),
  "",
  "## Monte Carlo results at the observed average imbalance",
  "",
  paste0(
    "The following values average each performance measure across the four truth ",
    "scenarios. Bias is computed from admissible and inadmissible terminal point ",
    "estimates when available; the admissible rate is reported separately."
  ),
  "",
  monte_carlo_key_lines,
  "",
  "### Stress condition: Sample 2 mean chi-square = 1.50",
  "",
  paste0(
    "The stress table averages absolute bias and RMSE across the four truth ",
    "scenarios. It shows how the models behave beyond the observed-average ",
    "imbalance."
  ),
  "",
  stress_key_lines,
  "",
  paste0("![Monte Carlo gap bias](", basename(monte_carlo_figure_path), ")"),
  "",
  "## Interpretation",
  "",
  paste0(
    "At the observed imbalance, the original model's mean absolute population ",
    "bias in the IQ-minus-Big-Five difference is ",
    paste(
      paste0(
        original_population_key$design_label,
        ": ",
        format_number(
          original_population_key$mean_absolute_difference_bias
        )
      ),
      collapse = "; "
    ),
    ". The oracle model's corresponding values are ",
    paste(
      paste0(
        oracle_population_key$design_label,
        ": ",
        format_number(oracle_population_key$mean_absolute_difference_bias)
      ),
      collapse = "; "
    ),
    "."
  ),
  "",
  paste0(
    "In Monte Carlo sampling, the original model's mean absolute difference bias ",
    "is ",
    paste(
      paste0(
        original_mc_key$design_label,
        ": ",
        format_number(original_mc_key$mean_absolute_difference_bias)
      ),
      collapse = "; "
    ),
    ", compared with ",
    paste(
      paste0(
        oracle_mc_key$design_label,
        ": ",
        format_number(oracle_mc_key$mean_absolute_difference_bias)
      ),
      collapse = "; "
    ),
    " for the oracle model."
  ),
  "",
  paste0(
    "The repository-faithful design is the relevant primary result because the ",
    "actual split-half model leaves IQ outside the Big Five sample factors. At the ",
    "observed imbalance, the original population model inflates the empirical ",
    "gene-score IQ-minus-Big-Five difference from 0.072 to 0.088 and changes the ",
    "empirical gene-set difference from -0.008 to 0.015. Both ratio models recover ",
    "the generating differences essentially exactly. When all six traits receive ",
    "the same imbalance, the original model shifts both correlation categories in ",
    "a more similar way and its difference bias is much smaller."
  ),
  "",
  paste0(
    "The estimated common-ratio and oracle fixed-ratio models have nearly the same ",
    "Monte Carlo bias and RMSE. This favors the estimated common ratio for the real ",
    "sensitivity analysis: it removes the systematic population bias without ",
    "requiring the chi-square transformation to be exactly correct. The fixed ",
    "chi-square model remains useful as a sharper, assumption-dependent benchmark."
  ),
  "",
  paste0(
    "Because method covariance is a product of method loadings, individual ",
    "sample-factor loadings are nonregular when their true values are all zero. ",
    "The simulation therefore reports their mean variance contribution and ",
    "empirical convergence rather than relying on ordinary Wald tests of zero ",
    "method loadings. The nominal chi-square rejection rate should likewise be ",
    "read as a diagnostic under this nonregular null."
  ),
  "",
  "## Output files",
  "",
  "- `population_results.csv`: all exact-population fits and bias measures.",
  "- `monte_carlo_replicates.csv`: every fitted Monte Carlo replication.",
  "- `monte_carlo_summary.csv`: bias, RMSE, dispersion, convergence, and fit summaries.",
  "- `truth_scenarios.csv`: the four latent-correlation truth scenarios.",
  "- `chi_square_grid.csv`: chi-square values and their implied loadings and ratios.",
  "- `population_category_bias.png` and `monte_carlo_gap_bias.png`: diagnostic figures.",
  "- `split_half_design_simulation_report.md`: this report."
)
writeLines(report_lines, report_path, useBytes = TRUE)

message("Wrote population results: ", population_path)
message("Wrote Monte Carlo replicates: ", monte_carlo_replicates_path)
message("Wrote Monte Carlo summary: ", monte_carlo_summary_path)
message("Wrote report: ", report_path)
