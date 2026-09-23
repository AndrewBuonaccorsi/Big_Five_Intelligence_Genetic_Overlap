#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(lavaan)
})

if (!identical(
  Sys.getenv("FULL_JACK_PAPER_RESULTS_MODE", unset = "0"),
  "1"
)) {
  stop(
    paste0(
      "This script is an internal helper for scripts/analysis.R. ",
      "Run the canonical analysis instead."
    ),
    call. = FALSE
  )
}

# ==============================================================================
# 1. Paths and configuration
# ==============================================================================

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) == 1L) {
  script_path <- normalizePath(
    sub("^--file=", "", script_argument),
    mustWork = TRUE
  )
} else {
  script_path <- normalizePath(
    "scripts/gene_set_enrichment_sem.R",
    mustWork = TRUE
  )
}

project_directory <- dirname(dirname(script_path))
arguments <- commandArgs(trailingOnly = TRUE)
model_variant <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  "both_samples"
}
if (!identical(model_variant, "both_samples")) {
  stop(
    paste0(
      "Only the primary 'both_samples' model is supported. The obsolete ",
      "'sample_one_only' result set has been retired."
    ),
    call. = FALSE
  )
}
include_sample_two_factor <- TRUE
analysis_kind <- Sys.getenv(
  "FULL_JACK_SEM_ANALYSIS_KIND",
  unset = "gene_sets"
)
if (!analysis_kind %in% c("gene_sets", "gene_scores")) {
  stop(
    "FULL_JACK_SEM_ANALYSIS_KIND must be 'gene_sets' or 'gene_scores'.",
    call. = FALSE
  )
}
is_gene_score_analysis <- identical(analysis_kind, "gene_scores")
artifact_base_name <- if (is_gene_score_analysis) {
  "gene_score_sem"
} else {
  "gene_set_enrichment_sem"
}
analysis_title <- if (is_gene_score_analysis) {
  "Gene-Score"
} else {
  "Gene-Set Enrichment"
}
analysis_description <- if (is_gene_score_analysis) {
  "gene-score"
} else {
  "gene-set enrichment"
}
general_factor_display_label <- if (is_gene_score_analysis) {
  "General gene score"
} else {
  "General gene-set enrichment"
}
analysis_output_label <- if (is_gene_score_analysis) {
  "Gene scores"
} else {
  "Pooled gene-set enrichment"
}
output_directory <- if (length(arguments) >= 1L) {
  path.expand(arguments[[1L]])
} else {
  file.path(project_directory, "output", "sem", artifact_base_name)
}
if (!grepl("^/", output_directory)) {
  output_directory <- file.path(project_directory, output_directory)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

default_correlation_path <- file.path(
  project_directory,
  "output",
  "correlation_matrices",
  if (is_gene_score_analysis) {
    "gene_score_correlation_matrix.csv"
  } else {
    "all_gene_sets_correlation_matrix.csv"
  }
)
default_sampling_vcov_path <- file.path(
  project_directory,
  "output",
  "sampling_covariances",
  if (is_gene_score_analysis) {
    "gene_score_correlation_sampling_vcov.csv"
  } else {
    "all_gene_sets_correlation_sampling_vcov.csv"
  }
)
default_jackknife_correlation_path <- file.path(
  project_directory,
  "output",
  if (is_gene_score_analysis) {
    "gene_scores"
  } else {
    "gene_sets"
  },
  if (is_gene_score_analysis) {
    "gene_score_jackknife_correlations.csv"
  } else {
    "all_gene_sets_jackknife_correlations.csv"
  }
)
correlation_path <- Sys.getenv(
  "FULL_JACK_SEM_CORRELATION_PATH",
  unset = default_correlation_path
)
sampling_vcov_path <- Sys.getenv(
  "FULL_JACK_SEM_SAMPLING_VCOV_PATH",
  unset = default_sampling_vcov_path
)
jackknife_correlation_path <- Sys.getenv(
  "FULL_JACK_SEM_JACKKNIFE_CORRELATION_PATH",
  unset = default_jackknife_correlation_path
)
workbook_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_results.xlsx")
)
tikz_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_model.tex")
)
diagram_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_model.pdf")
)
simplified_tikz_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_model_simplified.tex")
)
simplified_diagram_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_model_simplified.pdf")
)
model_fit_csv_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_model_fit.csv")
)
model_comparison_csv_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_model_comparison.csv")
)
structural_comparison_csv_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_structural_comparison.csv")
)
model_fit_tex_path <- file.path(
  output_directory,
  paste0(artifact_base_name, "_model_fit.tex")
)

required_files <- c(
  correlation_path,
  sampling_vcov_path,
  jackknife_correlation_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing required input file(s): ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

jackknife_blocks <- 200L
excluded_collection <- "neural_brain"
excluded_gene_set <- paste0(
  "GOBP_NEGATIVE_REGULATION_OF_INTRACELLULAR_",
  "LIPID_TRANSPORT"
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
trait_factor_order <- c(
  "iq", "extra", "agree", "consc", "neurot", "open"
)
trait_factor_labels <- c(
  iq = "IQ",
  extra = "Extraversion",
  agree = "Agreeableness",
  consc = "Conscientiousness",
  neurot = "Neuroticism",
  open = "Openness"
)

model_syntax <- '
    agree =~ equal_agr*agr1 + equal_agr*agr2
    consc =~ equal_con*con1 + equal_con*con2
    extra =~ equal_ext*ext1 + equal_ext*ext2
    neurot =~ equal_neu*neu1 + equal_neu*neu2
    open =~ equal_ope*ope1 + equal_ope*ope2
    iq =~ equal_iq*iq_female + equal_iq*iq_male

    gene_set_enrichment =~
      general_iq*iq +
      general_agree*agree +
      general_consc*consc +
      general_extra*extra +
      general_neurot*neurot +
      general_open*open

    sample_one =~ agr1 + con1 + ext1 + neu1 + ope1
    sample_two =~ agr2 + con2 + ext2 + neu2 + ope2

    gene_set_enrichment ~~ 0*sample_one
    gene_set_enrichment ~~ 0*sample_two
    agree + consc + extra + neurot + open + iq ~~ 0*sample_one
    agree + consc + extra + neurot + open + iq ~~ 0*sample_two
    sample_one ~~ 0*sample_two

    corr_agree_consc := general_agree * general_consc
    corr_agree_extra := general_agree * general_extra
    corr_agree_neurot := general_agree * general_neurot
    corr_agree_open := general_agree * general_open
    corr_agree_iq := general_agree * general_iq
    corr_consc_extra := general_consc * general_extra
    corr_consc_neurot := general_consc * general_neurot
    corr_consc_open := general_consc * general_open
    corr_consc_iq := general_consc * general_iq
    corr_extra_neurot := general_extra * general_neurot
    corr_extra_open := general_extra * general_open
    corr_extra_iq := general_extra * general_iq
    corr_neurot_open := general_neurot * general_open
    corr_neurot_iq := general_neurot * general_iq
    corr_open_iq := general_open * general_iq
  '
  correlated_model_syntax <- '
    agree =~ equal_agr*agr1 + equal_agr*agr2
    consc =~ equal_con*con1 + equal_con*con2
    extra =~ equal_ext*ext1 + equal_ext*ext2
    neurot =~ equal_neu*neu1 + equal_neu*neu2
    open =~ equal_ope*ope1 + equal_ope*ope2
    iq =~ equal_iq*iq_female + equal_iq*iq_male

    iq ~~ agree + consc + extra + neurot + open
    agree ~~ consc + extra + neurot + open
    consc ~~ extra + neurot + open
    extra ~~ neurot + open
    neurot ~~ open

    sample_one =~ agr1 + con1 + ext1 + neu1 + ope1
    sample_two =~ agr2 + con2 + ext2 + neu2 + ope2

    agree + consc + extra + neurot + open + iq ~~ 0*sample_one
    agree + consc + extra + neurot + open + iq ~~ 0*sample_two
    sample_one ~~ 0*sample_two
  '
  latent_factor_names <- c(
    "gene_set_enrichment",
    "iq",
    "agree",
    "consc",
    "extra",
    "neurot",
    "open",
    "sample_one",
    "sample_two"
  )
  model_title <- paste0(
    "Second-Order ",
    analysis_title,
    " SEM: Both Sample Factors"
  )
  model_description <- paste0(
    "Second-order six-trait model with equal paired loadings and orthogonal ",
    "Sample 1 and Sample 2 method factors"
  )
  model_subtitle <- paste0(
    "Within-trait indicator loadings are equal across samples; sample factors ",
    "capture Big Five sample-specific covariance; WLS inference uses the ",
    "full delete-one jackknife covariance."
  )
  factor_correlation_note <- paste0(
    "Sample factors are mutually orthogonal and orthogonal to the general ",
    "and six trait factors"
  )
  fitted_correlation_note <- paste0(
    "Correlations implied by the fitted second-order trait-factor model."
  )
  r_squared_note <- paste0(
    "Variance in each observed ",
    analysis_description,
    " profile explained by the ",
    "trait, general, and sample-factor structure."
  )
  factor_definitions <- data.frame(
    factor = c(
      "gene_set_enrichment",
      "agree",
      "consc",
      "extra",
      "neurot",
      "open",
      "iq",
      "sample_one",
      "sample_two",
      "Equal indicator loadings",
      "Orthogonality"
    ),
    definition = c(
      "Second-order factor loading on the six trait factors.",
      "Loads equally on agreeableness indicators from Samples 1 and 2.",
      "Loads equally on conscientiousness indicators from Samples 1 and 2.",
      "Loads equally on extraversion indicators from Samples 1 and 2.",
      "Loads equally on neuroticism indicators from Samples 1 and 2.",
      "Loads equally on openness indicators from Samples 1 and 2.",
      "Loads equally on female and male IQ indicators.",
      "Loads on the five sample-one personality indicators only.",
      "Loads on the five sample-two personality indicators only.",
      "The two loadings within every trait factor share one parameter.",
      paste0(
        "Sample factors are mutually orthogonal and orthogonal to the ",
        "general and trait factors."
      )
    ),
    stringsAsFactors = FALSE
  )


# ==============================================================================
# 2. Helpers
# ==============================================================================

validate_symmetric_matrix <- function(x, context, tolerance = 1e-10) {
  if (!is.matrix(x) || nrow(x) != ncol(x)) {
    stop(context, " is not a square matrix.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop(context, " contains non-finite values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = tolerance))) {
    stop(context, " is not symmetric.", call. = FALSE)
  }
  invisible(TRUE)
}

matrix_payload <- function(x) {
  list(
    row_names = rownames(x),
    column_names = colnames(x),
    values = unname(x)
  )
}

records_from_data_frame <- function(x) {
  rownames(x) <- NULL
  lapply(
    seq_len(nrow(x)),
    function(index) {
      as.list(x[index, , drop = FALSE])
    }
  )
}

format_loading <- function(loadings, factor_name, indicator_name) {
  value <- loadings$std_all[
    loadings$factor == factor_name &
      loadings$indicator == indicator_name
  ]
  if (length(value) != 1L || !is.finite(value)) {
    stop(
      "Could not find a finite standardized loading for ",
      factor_name,
      " -> ",
      indicator_name,
      call. = FALSE
    )
  }
  sprintf("%.2f", value)
}

format_loading_with_se <- function(
  standardized_loadings,
  factor_name,
  indicator_name
) {
  loading_row <- standardized_loadings[
    standardized_loadings$factor == factor_name &
      standardized_loadings$indicator == indicator_name,
    ,
    drop = FALSE
  ]
  if (
    nrow(loading_row) != 1L ||
      !is.finite(loading_row$estimate[[1L]])
  ) {
    stop(
      "Could not find a finite standardized loading for ",
      factor_name,
      " -> ",
      indicator_name,
      call. = FALSE
    )
  }
  standard_error <- loading_row$se[[1L]]
  standard_error_label <- if (is.finite(standard_error)) {
    sprintf("%.2f", standard_error)
  } else {
    "--"
  }
  sprintf(
    "\\shortstack{%.2f\\\\(%s)}",
    loading_row$estimate[[1L]],
    standard_error_label
  )
}

tex_indicator_label <- c(
  agr1 = "Agree.\\\\S1",
  con1 = "Consc.\\\\S1",
  ext1 = "Extra.\\\\S1",
  neu1 = "Neurot.\\\\S1",
  ope1 = "Open.\\\\S1",
  iq_female = "IQ\\\\female",
  iq_male = "IQ\\\\male",
  agr2 = "Agree.\\\\S2",
  con2 = "Consc.\\\\S2",
  ext2 = "Extra.\\\\S2",
  neu2 = "Neurot.\\\\S2",
  ope2 = "Open.\\\\S2"
)


# ==============================================================================
# 3. Read, align, and validate the input matrices
# ==============================================================================

correlation_matrix <- as.matrix(
  read.csv(
    correlation_path,
    row.names = 1L,
    check.names = FALSE
  )
)
storage.mode(correlation_matrix) <- "double"

if (!setequal(rownames(correlation_matrix), long_names) ||
    !setequal(colnames(correlation_matrix), long_names)) {
  stop(
    "The pooled correlation matrix does not contain the expected 12 indicators.",
    call. = FALSE
  )
}
correlation_matrix <- correlation_matrix[long_names, long_names, drop = FALSE]
dimnames(correlation_matrix) <- list(short_names, short_names)

validate_symmetric_matrix(
  correlation_matrix,
  paste0(analysis_title, " correlation matrix")
)
if (max(abs(diag(correlation_matrix) - 1)) > 1e-10) {
  stop("The correlation-matrix diagonal is not one.", call. = FALSE)
}

pair_indices <- t(combn(seq_along(short_names), 2L))
short_pair_names <- paste(
  short_names[pair_indices[, 1L]],
  short_names[pair_indices[, 2L]],
  sep = "__"
)

sampling_vcov <- as.matrix(
  read.csv(
    sampling_vcov_path,
    row.names = 1L,
    check.names = FALSE
  )
)
storage.mode(sampling_vcov) <- "double"

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
    "The sampling covariance matrix lacks one or more expected correlations.",
    call. = FALSE
  )
}
sampling_vcov <- sampling_vcov[
  source_pair_names,
  source_pair_names,
  drop = FALSE
]
dimnames(sampling_vcov) <- list(short_pair_names, short_pair_names)

jackknife_correlation_table <- read.csv(
  jackknife_correlation_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (
  !"block" %in% names(jackknife_correlation_table) ||
    nrow(jackknife_correlation_table) != jackknife_blocks ||
    anyNA(jackknife_correlation_table$block) ||
    anyDuplicated(jackknife_correlation_table$block) ||
    !setequal(
      as.integer(jackknife_correlation_table$block),
      seq_len(jackknife_blocks)
    )
) {
  stop(
    "The jackknife correlation table must contain blocks 1 through 200.",
    call. = FALSE
  )
}
jackknife_correlation_table <- jackknife_correlation_table[
  match(
    seq_len(jackknife_blocks),
    as.integer(jackknife_correlation_table$block)
  ),
  ,
  drop = FALSE
]
source_jackknife_pair_names <- vapply(
  seq_len(nrow(pair_indices)),
  function(index) {
    resolve_pair_name(
      long_names[pair_indices[index, 1L]],
      long_names[pair_indices[index, 2L]],
      names(jackknife_correlation_table)
    )
  },
  character(1)
)
if (anyNA(source_jackknife_pair_names)) {
  stop(
    "The jackknife correlation table lacks one or more expected pairs.",
    call. = FALSE
  )
}
jackknife_correlation_estimates <- as.matrix(
  jackknife_correlation_table[
    ,
    source_jackknife_pair_names,
    drop = FALSE
  ]
)
storage.mode(jackknife_correlation_estimates) <- "double"
dimnames(jackknife_correlation_estimates) <- list(
  sprintf("jk%03d", seq_len(jackknife_blocks)),
  short_pair_names
)
if (!all(is.finite(jackknife_correlation_estimates))) {
  stop(
    "The jackknife correlation table contains non-finite estimates.",
    call. = FALSE
  )
}
jackknife_correlation_means <- colMeans(
  jackknife_correlation_estimates
)
centered_jackknife_correlations <- sweep(
  jackknife_correlation_estimates,
  MARGIN = 2L,
  STATS = jackknife_correlation_means,
  FUN = "-"
)
recovered_sampling_vcov <- (
  (jackknife_blocks - 1) / jackknife_blocks
) * crossprod(centered_jackknife_correlations)
dimnames(recovered_sampling_vcov) <- dimnames(sampling_vcov)
if (max(abs(recovered_sampling_vcov - sampling_vcov)) > 1e-10) {
  stop(
    paste0(
      "The jackknife replicates do not reproduce the supplied sampling ",
      "covariance matrix."
    ),
    call. = FALSE
  )
}

validate_symmetric_matrix(
  sampling_vcov,
  "Jackknife correlation sampling covariance matrix"
)

correlation_eigenvalues <- eigen(
  correlation_matrix,
  symmetric = TRUE,
  only.values = TRUE
)$values
sampling_vcov_eigenvalues <- eigen(
  sampling_vcov,
  symmetric = TRUE,
  only.values = TRUE
)$values

if (min(correlation_eigenvalues) <= 0) {
  stop("The pooled correlation matrix is not positive definite.", call. = FALSE)
}
if (min(sampling_vcov_eigenvalues) <= 0) {
  stop(
    "The jackknife sampling covariance matrix is not positive definite.",
    call. = FALSE
  )
}


# ==============================================================================
# 4. Fit the second-order trait model with orthogonal sample factors
# ==============================================================================

# lavaan defines NACOV as N times the covariance of the sample statistics.
# Multiplying the direct delete-one jackknife covariance by the same N supplied
# to lavaan makes the WLS discrepancy, test statistic, and standard errors use
# the jackknife covariance itself. The choice N = 200 is therefore a scaling
# convention, not a claim that there are 200 independent observations.
nacov <- sampling_vcov * jackknife_blocks
wls_weight <- solve(nacov)

fit_sem_model <- function(
  syntax,
  model_label,
  optim_bounds = NULL
) {
  captured_warnings <- character()
  fitted_model <- withCallingHandlers(
    lavaan::cfa(
      model = syntax,
      sample.cov = correlation_matrix,
      sample.nobs = jackknife_blocks,
      sample.cov.rescale = FALSE,
      std.lv = TRUE,
      estimator = "WLS",
      correlation = TRUE,
      WLS.V = wls_weight,
      NACOV = nacov,
      se = "standard",
      test = "standard",
      optim.bounds = optim_bounds,
      # Retain the terminal solution so every diagnostic remains auditable.
      check.gradient = FALSE
    ),
    warning = function(warning_condition) {
      captured_warnings <<- unique(c(
        captured_warnings,
        conditionMessage(warning_condition)
      ))
      invokeRestart("muffleWarning")
    }
  )

  model_converged <- isTRUE(
    lavaan::lavInspect(fitted_model, "converged")
  )
  model_post_check <- withCallingHandlers(
    isTRUE(lavaan::lavInspect(fitted_model, "post.check")),
    warning = function(warning_condition) {
      captured_warnings <<- unique(c(
        captured_warnings,
        conditionMessage(warning_condition)
      ))
      invokeRestart("muffleWarning")
    }
  )
  if (!model_converged) {
    stop(
      model_label,
      " did not return a terminal solution that can be reported.",
      call. = FALSE
    )
  }

  model_maximum_gradient <- max(abs(fitted_model@optim$dx))
  model_gradient_tolerance <- lavaan::lavOptions()$optim.dx.tol
  model_gradient_passed <-
    model_maximum_gradient <= model_gradient_tolerance

  list(
    fit = fitted_model,
    warnings = captured_warnings,
    converged = model_converged,
    post_check = model_post_check,
    maximum_absolute_gradient = model_maximum_gradient,
    gradient_tolerance = model_gradient_tolerance,
    gradient_check_passed = model_gradient_passed,
    admissible = (
      model_converged &&
        model_post_check &&
        model_gradient_passed
    )
  )
}

hierarchical_result <- fit_sem_model(
  model_syntax,
  "Hierarchical general-factor model"
)
lavaan_correlated_result <- fit_sem_model(
  correlated_model_syntax,
  "Unconstrained correlated-factors starting model"
)

fit <- hierarchical_result$fit
fit_warnings <- hierarchical_result$warnings
converged <- hierarchical_result$converged
post_check <- hierarchical_result$post_check
maximum_absolute_gradient <-
  hierarchical_result$maximum_absolute_gradient
gradient_tolerance <- hierarchical_result$gradient_tolerance
gradient_check_passed <- hierarchical_result$gradient_check_passed
solution_admissible <- hierarchical_result$admissible

lavaan_correlated_fit <- lavaan_correlated_result$fit
analysis_status <- if (solution_admissible) {
  paste0(
    "ADMISSIBLE - the model converged, passed the gradient check, and ",
    "contains no negative observed residual variances."
  )
} else if (include_sample_two_factor) {
  paste0(
    "INADMISSIBLE - report for diagnostic purposes only. ",
    "The second-order model with both sample factors failed an estimation ",
    "or post-estimation admissibility check."
  )
} else {
  paste0(
    "INADMISSIBLE - report for diagnostic purposes only. ",
    "The second-order sample-one-only structure failed an estimation or ",
    "post-estimation admissibility check."
  )
}
analysis_status_short <- if (solution_admissible) {
  "ADMISSIBLE"
} else if (include_sample_two_factor) {
  "INADMISSIBLE - both-sample model failed admissibility checks"
} else {
  "INADMISSIBLE - sample-one model failed admissibility checks"
}

fit_measure_names <- c(
  "chisq", "df", "pvalue",
  "cfi", "tli",
  "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue",
  "srmr"
)
all_fit_measures <- lavaan::fitMeasures(fit)
fit_measures <- all_fit_measures[fit_measure_names]

extract_parameter_estimates <- function(fitted_model) {
  estimates <- lavaan::parameterEstimates(
    fitted_model,
    standardized = TRUE,
    ci = TRUE,
    level = 0.95
  )
  estimates <- estimates[
    ,
    c(
      "lhs", "op", "rhs",
      "est", "se", "z", "pvalue",
      "ci.lower", "ci.upper",
      "std.lv", "std.all"
    )
  ]
  names(estimates) <- c(
    "lhs", "operator", "rhs",
    "estimate", "se", "z", "p_value",
    "ci_lower", "ci_upper",
    "std_lv", "std_all"
  )
  estimates$significant_05 <- ifelse(
    is.na(estimates$p_value),
    NA,
    estimates$p_value < 0.05
  )
  estimates
}

parameter_estimates <- extract_parameter_estimates(fit)
lavaan_correlated_parameter_estimates <-
  extract_parameter_estimates(lavaan_correlated_fit)

# The unconstrained lavaan correlated-factor solution is useful for starting
# values, but lavaan does not constrain the latent correlation matrix to remain
# positive definite. Fit the same 31-parameter covariance structure directly
# with a Cholesky-parameterized trait correlation matrix and admissible observed
# residual variances.
correlation_matrix_from_parameters <- function(parameters, dimension) {
  lower <- diag(dimension)
  lower[lower.tri(lower)] <- parameters
  covariance <- tcrossprod(lower)
  covariance / sqrt(outer(diag(covariance), diag(covariance)))
}

matrix_to_correlation_parameters <- function(x) {
  positive_definite <- as.matrix(
    Matrix::nearPD(
      x,
      corr = TRUE,
      keepDiag = TRUE,
      eig.tol = 1e-10,
      posd.tol = 1e-10
    )$mat
  )
  lower <- t(chol(positive_definite))
  unit_lower <- sweep(lower, 1L, diag(lower), "/")
  unit_lower[lower.tri(unit_lower)]
}

indicator_trait <- c(
  "agree", "agree",
  "consc", "consc",
  "extra", "extra",
  "neurot", "neurot",
  "open", "open",
  "iq", "iq"
)
indicator_sample <- c(
  1L, 2L,
  1L, 2L,
  1L, 2L,
  1L, 2L,
  1L, 2L,
  NA_integer_, NA_integer_
)
personality_indicator_indices <- which(!is.na(indicator_sample))
correlated_parameter_count <-
  length(trait_factor_order) +
    length(personality_indicator_indices) +
    choose(length(trait_factor_order), 2L)

decode_correlated_parameters <- function(parameters) {
  trait_count <- length(trait_factor_order)
  method_count <- length(personality_indicator_indices)
  trait_raw <- parameters[seq_len(trait_count)]
  method_raw <- parameters[
    trait_count + seq_len(method_count)
  ]
  correlation_raw <- parameters[
    trait_count + method_count +
      seq_len(choose(trait_count, 2L))
  ]

  trait_loadings <- 0.999 * tanh(trait_raw)
  names(trait_loadings) <- trait_factor_order
  method_cap <- 0.999 * sqrt(
    pmax(
      1e-10,
      1 - trait_loadings[
        indicator_trait[personality_indicator_indices]
      ]^2
    )
  )
  method_loadings <- method_cap * tanh(method_raw)
  names(method_loadings) <-
    short_names[personality_indicator_indices]
  trait_correlation <- correlation_matrix_from_parameters(
    correlation_raw,
    trait_count
  )
  dimnames(trait_correlation) <- list(
    trait_factor_order,
    trait_factor_order
  )

  measurement_matrix <- matrix(
    0,
    nrow = length(short_names),
    ncol = trait_count,
    dimnames = list(short_names, trait_factor_order)
  )
  measurement_matrix[
    cbind(
      seq_along(short_names),
      match(indicator_trait, trait_factor_order)
    )
  ] <- trait_loadings[indicator_trait]

  method_matrix <- matrix(
    0,
    nrow = length(short_names),
    ncol = 2L,
    dimnames = list(short_names, c("sample_one", "sample_two"))
  )
  method_matrix[
    cbind(
      personality_indicator_indices,
      indicator_sample[personality_indicator_indices]
    )
  ] <- method_loadings

  implied_correlation <-
    measurement_matrix %*%
      trait_correlation %*%
      t(measurement_matrix) +
      tcrossprod(method_matrix)
  residual_variances <- 1 - diag(implied_correlation)
  diag(implied_correlation) <- 1
  dimnames(implied_correlation) <- list(short_names, short_names)

  list(
    trait_loadings = trait_loadings,
    method_loadings = method_loadings,
    trait_correlation = trait_correlation,
    measurement_matrix = measurement_matrix,
    method_matrix = method_matrix,
    residual_variances = residual_variances,
    implied_correlation = implied_correlation
  )
}

sampling_precision <- solve(sampling_vcov)
wls_test_scale <- (jackknife_blocks - 1) / jackknife_blocks
wls_objective_precision <- wls_test_scale * sampling_precision
observed_pair_vector <- correlation_matrix[pair_indices]
correlated_objective <- function(parameters) {
  implied_pair_vector <-
    decode_correlated_parameters(parameters)$implied_correlation[
      pair_indices
    ]
  residual <- observed_pair_vector - implied_pair_vector
  as.numeric(crossprod(
    residual,
    wls_objective_precision %*% residual
  ))
}

starting_trait_loadings <- vapply(
  trait_factor_order,
  function(trait) {
    rows <- lavaan_correlated_parameter_estimates[
      lavaan_correlated_parameter_estimates$operator == "=~" &
        lavaan_correlated_parameter_estimates$lhs == trait,
      ,
      drop = FALSE
    ]
    rows$estimate[[1L]]
  },
  numeric(1)
)
starting_trait_raw <- atanh(
  pmin(0.995, pmax(-0.995, starting_trait_loadings / 0.999))
)
starting_method_loadings <- vapply(
  personality_indicator_indices,
  function(indicator_index) {
    factor_name <- paste0(
      "sample_",
      if (indicator_sample[[indicator_index]] == 1L) "one" else "two"
    )
    rows <- lavaan_correlated_parameter_estimates[
      lavaan_correlated_parameter_estimates$operator == "=~" &
        lavaan_correlated_parameter_estimates$lhs == factor_name &
        lavaan_correlated_parameter_estimates$rhs ==
          short_names[[indicator_index]],
      ,
      drop = FALSE
    ]
    rows$estimate[[1L]]
  },
  numeric(1)
)
starting_method_cap <- 0.999 * sqrt(
  pmax(
    1e-10,
    1 - starting_trait_loadings[
      indicator_trait[personality_indicator_indices]
    ]^2
  )
)
starting_method_raw <- atanh(
  pmin(
    0.995,
    pmax(-0.995, starting_method_loadings / starting_method_cap)
  )
)
lavaan_starting_correlation <-
  lavaan::lavInspect(lavaan_correlated_fit, "cor.lv")[
    trait_factor_order,
    trait_factor_order,
    drop = FALSE
  ]
hierarchical_starting_correlation <-
  lavaan::lavInspect(fit, "cor.lv")[
    trait_factor_order,
    trait_factor_order,
    drop = FALSE
  ]
correlated_starts <- list(
  c(
    starting_trait_raw,
    starting_method_raw,
    matrix_to_correlation_parameters(lavaan_starting_correlation)
  ),
  c(
    starting_trait_raw,
    starting_method_raw,
    matrix_to_correlation_parameters(hierarchical_starting_correlation)
  ),
  c(
    starting_trait_raw,
    starting_method_raw,
    rep(0, choose(length(trait_factor_order), 2L))
  )
)
correlated_fits <- lapply(
  correlated_starts,
  function(start) {
    nlminb(
      start = pmin(10, pmax(-10, start)),
      objective = correlated_objective,
      lower = rep(-10, correlated_parameter_count),
      upper = rep(10, correlated_parameter_count),
      control = list(iter.max = 6000L, eval.max = 15000L)
    )
  }
)
correlated_objectives <- vapply(
  correlated_fits,
  `[[`,
  numeric(1),
  "objective"
)
correlated_optimizer <-
  correlated_fits[[which.min(correlated_objectives)]]
# Refine the best multistart solution under stricter stopping criteria. The
# correlated model is intentionally parameterized on an unconstrained scale,
# so this final pass is useful when the trait-correlation matrix is close to
# singular without actually reaching the positive-definiteness boundary.
correlated_refinement <- nlminb(
  start = correlated_optimizer$par,
  objective = correlated_objective,
  lower = rep(-10, correlated_parameter_count),
  upper = rep(10, correlated_parameter_count),
  control = list(
    iter.max = 20000L,
    eval.max = 50000L,
    rel.tol = 1e-12,
    x.tol = 1e-10
  )
)
if (
  is.finite(correlated_refinement$objective) &&
    correlated_refinement$objective <= correlated_optimizer$objective + 1e-10
) {
  correlated_optimizer <- correlated_refinement
}
correlated_parameters <- correlated_optimizer$par
correlated_decoded <-
  decode_correlated_parameters(correlated_parameters)
correlated_maximum_gradient <- max(
  abs(numDeriv::grad(correlated_objective, correlated_parameters))
)
correlated_jacobian <- numDeriv::jacobian(
  function(parameters) {
    decode_correlated_parameters(parameters)$implied_correlation[
      pair_indices
    ]
  },
  correlated_parameters
)
correlated_information <- crossprod(
  correlated_jacobian,
  sampling_precision %*% correlated_jacobian
)
correlated_information_eigenvalues <- eigen(
  correlated_information,
  symmetric = TRUE,
  only.values = TRUE
)$values
if (min(correlated_information_eigenvalues) <= 0) {
  stop(
    "The admissible correlated-factor information matrix is not ",
    "positive definite.",
    call. = FALSE
  )
}
correlated_parameter_vcov <- solve(correlated_information)

correlated_loading_function <- function(parameters) {
  decoded <- decode_correlated_parameters(parameters)
  c(
    decoded$trait_loadings[indicator_trait],
    decoded$method_loadings
  )
}
correlated_loading_estimates <-
  correlated_loading_function(correlated_parameters)
correlated_loading_jacobian <- numDeriv::jacobian(
  correlated_loading_function,
  correlated_parameters
)
correlated_loading_se <- sqrt(diag(
  correlated_loading_jacobian %*%
    correlated_parameter_vcov %*%
    t(correlated_loading_jacobian)
))

trait_pair_indices <- t(combn(seq_along(trait_factor_order), 2L))
correlated_correlation_function <- function(parameters) {
  decode_correlated_parameters(parameters)$trait_correlation[
    trait_pair_indices
  ]
}
correlated_correlation_estimates <-
  correlated_correlation_function(correlated_parameters)
correlated_correlation_jacobian <- numDeriv::jacobian(
  correlated_correlation_function,
  correlated_parameters
)
correlated_correlation_se <- sqrt(diag(
  correlated_correlation_jacobian %*%
    correlated_parameter_vcov %*%
    t(correlated_correlation_jacobian)
))

correlated_residual_function <- function(parameters) {
  decode_correlated_parameters(parameters)$residual_variances
}
correlated_residual_estimates <-
  correlated_residual_function(correlated_parameters)
correlated_residual_jacobian <- numDeriv::jacobian(
  correlated_residual_function,
  correlated_parameters
)
correlated_residual_se <- sqrt(diag(
  correlated_residual_jacobian %*%
    correlated_parameter_vcov %*%
    t(correlated_residual_jacobian)
))

parameter_statistics <- function(estimates, standard_errors) {
  z_values <- estimates / standard_errors
  p_values <- 2 * pnorm(-abs(z_values))
  data.frame(
    estimate = estimates,
    se = standard_errors,
    z = z_values,
    p_value = p_values,
    ci_lower = estimates - qnorm(0.975) * standard_errors,
    ci_upper = estimates + qnorm(0.975) * standard_errors,
    std_lv = estimates,
    std_all = estimates,
    significant_05 = p_values < 0.05,
    stringsAsFactors = FALSE
  )
}

trait_loading_statistics <- parameter_statistics(
  correlated_loading_estimates[seq_along(short_names)],
  correlated_loading_se[seq_along(short_names)]
)
method_loading_positions <-
  length(short_names) + seq_along(personality_indicator_indices)
method_loading_statistics <- parameter_statistics(
  correlated_loading_estimates[method_loading_positions],
  correlated_loading_se[method_loading_positions]
)
correlation_statistics <- parameter_statistics(
  correlated_correlation_estimates,
  correlated_correlation_se
)
residual_statistics <- parameter_statistics(
  correlated_residual_estimates,
  correlated_residual_se
)

correlated_parameter_estimates <- rbind(
  data.frame(
    lhs = indicator_trait,
    operator = "=~",
    rhs = short_names,
    trait_loading_statistics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    lhs = ifelse(
      indicator_sample[personality_indicator_indices] == 1L,
      "sample_one",
      "sample_two"
    ),
    operator = "=~",
    rhs = short_names[personality_indicator_indices],
    method_loading_statistics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    lhs = trait_factor_order[trait_pair_indices[, 1L]],
    operator = "~~",
    rhs = trait_factor_order[trait_pair_indices[, 2L]],
    correlation_statistics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    lhs = short_names,
    operator = "~~",
    rhs = short_names,
    residual_statistics,
    stringsAsFactors = FALSE
  )
)

correlated_fitted_matrix <- correlated_decoded$implied_correlation
correlated_residual_matrix <-
  correlation_matrix - correlated_fitted_matrix
correlated_chisq <- correlated_optimizer$objective
correlated_df <- length(observed_pair_vector) -
  correlated_parameter_count
correlated_srmr <- sqrt(
  2 * sum(
    correlated_residual_matrix[upper.tri(correlated_residual_matrix)]^2
  ) /
    (nrow(correlated_residual_matrix) *
      (nrow(correlated_residual_matrix) + 1))
)
correlated_fit_measures <- c(
  chisq = correlated_chisq,
  df = correlated_df,
  pvalue = pchisq(
    correlated_chisq,
    df = correlated_df,
    lower.tail = FALSE
  ),
  cfi = unname(lavaan:::lav_fit_cfi(
    X2 = correlated_chisq,
    df = correlated_df,
    X2.null = all_fit_measures[["baseline.chisq"]],
    df.null = all_fit_measures[["baseline.df"]]
  )),
  tli = unname(lavaan:::lav_fit_tli(
    X2 = correlated_chisq,
    df = correlated_df,
    X2.null = all_fit_measures[["baseline.chisq"]],
    df.null = all_fit_measures[["baseline.df"]]
  )),
  rmsea = unname(lavaan:::lav_fit_rmsea(
    X2 = correlated_chisq,
    df = correlated_df,
    N = jackknife_blocks
  )),
  rmsea.ci.lower = unname(lavaan:::lav_fit_rmsea_ci(
    X2 = correlated_chisq,
    df = correlated_df,
    N = jackknife_blocks,
    level = 0.90
  )[["rmsea.ci.lower"]]),
  rmsea.ci.upper = unname(lavaan:::lav_fit_rmsea_ci(
    X2 = correlated_chisq,
    df = correlated_df,
    N = jackknife_blocks,
    level = 0.90
  )[["rmsea.ci.upper"]]),
  rmsea.pvalue = unname(lavaan:::lav_fit_rmsea_closefit(
    X2 = correlated_chisq,
    df = correlated_df,
    N = jackknife_blocks,
    rmsea.h0 = 0.05
  )),
  srmr = correlated_srmr
)

correlated_trait_minimum_eigenvalue <- min(eigen(
  correlated_decoded$trait_correlation,
  symmetric = TRUE,
  only.values = TRUE
)$values)
# This gradient is evaluated in the unconstrained transformed coordinates.
# Near the positive-definiteness boundary those coordinates are poorly scaled,
# so use a transparent numerical tolerance alongside optimizer convergence and
# the substantive covariance-space admissibility checks below.
correlated_gradient_tolerance <- 1e-2
correlated_solution_admissible <- (
  correlated_optimizer$convergence == 0L &&
    correlated_maximum_gradient <= correlated_gradient_tolerance &&
    correlated_trait_minimum_eigenvalue > 1e-6 &&
    min(correlated_decoded$residual_variances) > 1e-8
)
correlated_result <- list(
  warnings = character(),
  converged = correlated_optimizer$convergence == 0L,
  post_check = (
    correlated_trait_minimum_eigenvalue > 1e-6 &&
      min(correlated_decoded$residual_variances) > 1e-8
  ),
  maximum_absolute_gradient = correlated_maximum_gradient,
  gradient_tolerance = correlated_gradient_tolerance,
  gradient_check_passed =
    correlated_maximum_gradient <= correlated_gradient_tolerance,
  admissible = correlated_solution_admissible,
  iterations = correlated_optimizer$iterations,
  free_parameters = correlated_parameter_count,
  minimum_trait_correlation_eigenvalue =
    correlated_trait_minimum_eigenvalue
)

loadings <- parameter_estimates[
  parameter_estimates$operator == "=~",
  c(
    "lhs", "rhs",
    "estimate", "se", "z", "p_value",
    "ci_lower", "ci_upper",
    "std_lv", "std_all", "significant_05"
  )
]
names(loadings)[1:2] <- c("factor", "indicator")

standardized_loadings <- suppressWarnings(
  lavaan::standardizedSolution(
    fit,
    type = "std.all",
    se = TRUE,
    zstat = TRUE,
    pvalue = TRUE,
    ci = TRUE,
    level = 0.95,
    remove.eq = FALSE
  )
)
standardized_loadings <- standardized_loadings[
  standardized_loadings$op == "=~",
  c("lhs", "rhs", "est.std", "se"),
  drop = FALSE
]
names(standardized_loadings) <- c(
  "factor",
  "indicator",
  "estimate",
  "se"
)

r_squared <- lavaan::lavInspect(fit, "rsquare")
r_squared <- r_squared[short_names]
r_squared_table <- data.frame(
  indicator = names(r_squared),
  r_squared = as.numeric(r_squared),
  status = ifelse(
    is.finite(as.numeric(r_squared)),
    "Available",
    "Not defined because the fitted residual variance is inadmissible"
  ),
  stringsAsFactors = FALSE
)

fitted_matrix <- lavaan::fitted(fit)$cov
fitted_matrix <- fitted_matrix[short_names, short_names, drop = FALSE]
residual_matrix <- correlation_matrix - fitted_matrix
hierarchical_direct_chisq <- as.numeric(crossprod(
  residual_matrix[pair_indices],
  wls_objective_precision %*% residual_matrix[pair_indices]
))
if (abs(hierarchical_direct_chisq - fit_measures[["chisq"]]) > 1e-6) {
  stop(
    "The direct WLS objective (",
    format(hierarchical_direct_chisq, digits = 12),
    ") does not reproduce lavaan's chi-square (",
    format(fit_measures[["chisq"]], digits = 12),
    ").",
    call. = FALSE
  )
}
hierarchical_srmr_direct <- sqrt(
  2 * sum(residual_matrix[upper.tri(residual_matrix)]^2) /
    (nrow(residual_matrix) * (nrow(residual_matrix) + 1))
)
if (abs(hierarchical_srmr_direct - fit_measures[["srmr"]]) > 1e-8) {
  stop("The direct SRMR calculation does not reproduce lavaan.", call. = FALSE)
}

pair_diagnostics <- data.frame(
  pair = short_pair_names,
  indicator_1 = short_names[pair_indices[, 1L]],
  indicator_2 = short_names[pair_indices[, 2L]],
  observed_correlation = correlation_matrix[pair_indices],
  fitted_correlation = fitted_matrix[pair_indices],
  residual_correlation = residual_matrix[pair_indices],
  jackknife_se = sqrt(diag(sampling_vcov)),
  stringsAsFactors = FALSE
)
pair_diagnostics$marginal_residual_z <-
  pair_diagnostics$residual_correlation /
    pair_diagnostics$jackknife_se
pair_diagnostics$absolute_marginal_residual_z <-
  abs(pair_diagnostics$marginal_residual_z)
pair_diagnostics <- pair_diagnostics[
  order(
    pair_diagnostics$absolute_marginal_residual_z,
    decreasing = TRUE
  ),
]
rownames(pair_diagnostics) <- NULL

latent_correlation <- lavaan::lavInspect(fit, "cor.lv")
latent_correlation <- latent_correlation[
  latent_factor_names,
  latent_factor_names,
  drop = FALSE
]

correlated_pair_diagnostics <- data.frame(
  pair = short_pair_names,
  indicator_1 = short_names[pair_indices[, 1L]],
  indicator_2 = short_names[pair_indices[, 2L]],
  observed_correlation = correlation_matrix[pair_indices],
  fitted_correlation = correlated_fitted_matrix[pair_indices],
  residual_correlation = correlated_residual_matrix[pair_indices],
  jackknife_se = sqrt(diag(sampling_vcov)),
  stringsAsFactors = FALSE
)
correlated_pair_diagnostics$marginal_residual_z <-
  correlated_pair_diagnostics$residual_correlation /
    correlated_pair_diagnostics$jackknife_se
correlated_pair_diagnostics$absolute_marginal_residual_z <-
  abs(correlated_pair_diagnostics$marginal_residual_z)
correlated_pair_diagnostics <- correlated_pair_diagnostics[
  order(
    correlated_pair_diagnostics$absolute_marginal_residual_z,
    decreasing = TRUE
  ),
]
rownames(correlated_pair_diagnostics) <- NULL

hierarchical_trait_correlation <-
  lavaan::lavInspect(fit, "cor.lv")[
    trait_factor_order,
    trait_factor_order,
    drop = FALSE
  ]
correlated_trait_correlation <-
  correlated_decoded$trait_correlation
dimnames(hierarchical_trait_correlation) <- list(
  unname(trait_factor_labels[trait_factor_order]),
  unname(trait_factor_labels[trait_factor_order])
)
dimnames(correlated_trait_correlation) <- list(
  unname(trait_factor_labels[trait_factor_order]),
  unname(trait_factor_labels[trait_factor_order])
)

comparison_delta_chisq <-
  unname(fit_measures[["chisq"]]) -
    unname(correlated_fit_measures[["chisq"]])
comparison_delta_df <- as.integer(
  round(
    unname(fit_measures[["df"]]) -
      unname(correlated_fit_measures[["df"]])
  )
)
if (comparison_delta_chisq < -1e-8 || comparison_delta_df <= 0L) {
  stop(
    "The hierarchical and correlated models are not ordered as expected.",
    call. = FALSE
  )
}
if (comparison_delta_df != 9L) {
  stop(
    "Expected nine restrictions in the hierarchical model; found ",
    comparison_delta_df,
    ".",
    call. = FALSE
  )
}
comparison_p_value <- pchisq(
  comparison_delta_chisq,
  df = comparison_delta_df,
  lower.tail = FALSE
)
comparison_significant <- comparison_p_value < 0.05

model_fit_table <- data.frame(
  model = c(
    "Correlated-factors model",
    "Hierarchical general-factor model"
  ),
  chi_square = c(
    correlated_fit_measures[["chisq"]],
    fit_measures[["chisq"]]
  ),
  degrees_of_freedom = c(
    correlated_fit_measures[["df"]],
    fit_measures[["df"]]
  ),
  chi_square_p_value = c(
    correlated_fit_measures[["pvalue"]],
    fit_measures[["pvalue"]]
  ),
  cfi = c(
    correlated_fit_measures[["cfi"]],
    fit_measures[["cfi"]]
  ),
  tli = c(
    correlated_fit_measures[["tli"]],
    fit_measures[["tli"]]
  ),
  rmsea = c(
    correlated_fit_measures[["rmsea"]],
    fit_measures[["rmsea"]]
  ),
  rmsea_ci_90_lower = c(
    correlated_fit_measures[["rmsea.ci.lower"]],
    fit_measures[["rmsea.ci.lower"]]
  ),
  rmsea_ci_90_upper = c(
    correlated_fit_measures[["rmsea.ci.upper"]],
    fit_measures[["rmsea.ci.upper"]]
  ),
  rmsea_p_close = c(
    correlated_fit_measures[["rmsea.pvalue"]],
    fit_measures[["rmsea.pvalue"]]
  ),
  srmr = c(
    correlated_fit_measures[["srmr"]],
    fit_measures[["srmr"]]
  ),
  effective_free_parameters = c(
    correlated_parameter_count,
    length(observed_pair_vector) -
      as.integer(round(fit_measures[["df"]]))
  ),
  admissible = c(
    correlated_solution_admissible,
    solution_admissible
  ),
  stringsAsFactors = FALSE
)

comparison_conclusion <- if (comparison_significant) {
  paste0(
    "The hierarchical model fits significantly worse than the ",
    "correlated-factors model."
  )
} else {
  paste0(
    "There is no evidence that the hierarchical model fits worse than the ",
    "correlated-factors model."
  )
}
model_comparison_table <- data.frame(
  comparison = "Hierarchical general factor vs correlated factors",
  less_restricted_model = "Correlated-factors model",
  restricted_model = "Hierarchical general-factor model",
  delta_chi_square = comparison_delta_chisq,
  delta_degrees_of_freedom = comparison_delta_df,
  p_value = comparison_p_value,
  significant_05 = comparison_significant,
  conclusion = comparison_conclusion,
  stringsAsFactors = FALSE
)

trait_pair_indices <- t(combn(seq_along(trait_factor_order), 2L))
correlated_covariance_rows <- lapply(
  seq_len(nrow(trait_pair_indices)),
  function(index) {
    left <- trait_factor_order[trait_pair_indices[index, 1L]]
    right <- trait_factor_order[trait_pair_indices[index, 2L]]
    rows <- correlated_parameter_estimates[
      correlated_parameter_estimates$operator == "~~" &
        (
          (
            correlated_parameter_estimates$lhs == left &
              correlated_parameter_estimates$rhs == right
          ) |
            (
              correlated_parameter_estimates$lhs == right &
                correlated_parameter_estimates$rhs == left
            )
        ),
      ,
      drop = FALSE
    ]
    if (nrow(rows) != 1L) {
      stop(
        "Could not identify exactly one correlated-factor estimate for ",
        left,
        " and ",
        right,
        ".",
        call. = FALSE
      )
    }
    rows
  }
)
correlated_covariance_table <- do.call(
  rbind,
  correlated_covariance_rows
)
structural_comparison_table <- data.frame(
  pair = paste(
    unname(
      trait_factor_labels[
        trait_factor_order[trait_pair_indices[, 1L]]
      ]
    ),
    unname(
      trait_factor_labels[
        trait_factor_order[trait_pair_indices[, 2L]]
      ]
    ),
    sep = " - "
  ),
  trait_1 = unname(
    trait_factor_labels[
      trait_factor_order[trait_pair_indices[, 1L]]
    ]
  ),
  trait_2 = unname(
    trait_factor_labels[
      trait_factor_order[trait_pair_indices[, 2L]]
    ]
  ),
  correlated_estimate = correlated_covariance_table$std_lv,
  correlated_se = correlated_covariance_table$se,
  correlated_z = correlated_covariance_table$z,
  correlated_p_value = correlated_covariance_table$p_value,
  correlated_ci_lower = correlated_covariance_table$ci_lower,
  correlated_ci_upper = correlated_covariance_table$ci_upper,
  hierarchical_implied = hierarchical_trait_correlation[
    trait_pair_indices
  ],
  correlated_minus_hierarchical =
    correlated_covariance_table$std_lv -
      hierarchical_trait_correlation[trait_pair_indices],
  stringsAsFactors = FALSE
)

fit_statistics_table <- data.frame(
  metric = c(
    "Chi-square",
    "Degrees of freedom",
    "Chi-square p-value",
    "CFI",
    "TLI",
    "RMSEA",
    "RMSEA 90% CI lower",
    "RMSEA 90% CI upper",
    "RMSEA p(close)",
    "SRMR"
  ),
  value = as.numeric(fit_measures),
  stringsAsFactors = FALSE
)

observed_residual_rows <-
  parameter_estimates$operator == "~~" &
    parameter_estimates$lhs == parameter_estimates$rhs &
    parameter_estimates$lhs %in% short_names
observed_residual_variances <-
  parameter_estimates$estimate[observed_residual_rows]
negative_residual_variance_count <- sum(
  observed_residual_variances < 0,
  na.rm = TRUE
)
minimum_observed_residual_variance <- min(
  observed_residual_variances,
  na.rm = TRUE
)

correlated_observed_residual_rows <-
  correlated_parameter_estimates$operator == "~~" &
    correlated_parameter_estimates$lhs ==
      correlated_parameter_estimates$rhs &
    correlated_parameter_estimates$lhs %in% short_names
correlated_observed_residual_variances <-
  correlated_parameter_estimates$estimate[
    correlated_observed_residual_rows
  ]
correlated_negative_residual_variance_count <- sum(
  correlated_observed_residual_variances < 0,
  na.rm = TRUE
)
correlated_minimum_observed_residual_variance <- min(
  correlated_observed_residual_variances,
  na.rm = TRUE
)

correlated_loadings <- correlated_parameter_estimates[
  correlated_parameter_estimates$operator == "=~",
  c(
    "lhs", "rhs",
    "estimate", "se", "z", "p_value",
    "ci_lower", "ci_upper",
    "std_lv", "std_all", "significant_05"
  )
]
names(correlated_loadings)[1:2] <- c("factor", "indicator")

comparison_loadings <- rbind(
  data.frame(
    model = "Hierarchical general-factor model",
    loadings,
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "Correlated-factors model",
    correlated_loadings,
    stringsAsFactors = FALSE
  )
)
comparison_parameter_estimates <- rbind(
  data.frame(
    model = "Hierarchical general-factor model",
    parameter_estimates,
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "Correlated-factors model",
    correlated_parameter_estimates,
    stringsAsFactors = FALSE
  )
)
comparison_pair_diagnostics <- rbind(
  data.frame(
    model = "Hierarchical general-factor model",
    pair_diagnostics,
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "Correlated-factors model",
    correlated_pair_diagnostics,
    stringsAsFactors = FALSE
  )
)

model_diagnostics_table <- data.frame(
  model = c(
    "Correlated-factors model",
    "Hierarchical general-factor model"
  ),
  admissible = c(
    correlated_solution_admissible,
    solution_admissible
  ),
  converged = c(
    correlated_result$converged,
    converged
  ),
  gradient_check_passed = c(
    correlated_result$gradient_check_passed,
    gradient_check_passed
  ),
  maximum_absolute_gradient = c(
    correlated_result$maximum_absolute_gradient,
    maximum_absolute_gradient
  ),
  gradient_tolerance = c(
    correlated_result$gradient_tolerance,
    gradient_tolerance
  ),
  post_estimation_check = c(
    correlated_result$post_check,
    post_check
  ),
  optimizer_iterations = c(
    correlated_result$iterations,
    lavaan::lavInspect(fit, "iterations")
  ),
  effective_free_parameters = c(
    correlated_result$free_parameters,
    length(observed_pair_vector) -
      as.integer(round(fit_measures[["df"]]))
  ),
  negative_observed_residual_variances = c(
    correlated_negative_residual_variance_count,
    negative_residual_variance_count
  ),
  minimum_observed_residual_variance = c(
    correlated_minimum_observed_residual_variance,
    minimum_observed_residual_variance
  ),
  minimum_trait_correlation_eigenvalue = c(
    correlated_result$minimum_trait_correlation_eigenvalue,
    min(eigen(
      hierarchical_trait_correlation,
      symmetric = TRUE,
      only.values = TRUE
    )$values)
  ),
  maximum_absolute_residual_correlation = c(
    max(
      abs(
        correlated_residual_matrix[
          upper.tri(correlated_residual_matrix)
        ]
      )
    ),
    max(abs(residual_matrix[upper.tri(residual_matrix)]))
  ),
  maximum_absolute_marginal_residual_z = c(
    max(correlated_pair_diagnostics$absolute_marginal_residual_z),
    max(pair_diagnostics$absolute_marginal_residual_z)
  ),
  captured_warnings = c(
    if (length(correlated_result$warnings) == 0L) {
      "None"
    } else {
      paste(correlated_result$warnings, collapse = " | ")
    },
    if (length(fit_warnings) == 0L) {
      "None"
    } else {
      paste(fit_warnings, collapse = " | ")
    }
  ),
  stringsAsFactors = FALSE
)

diagnostics_table <- data.frame(
  diagnostic = c(
    "Solution status",
    "Optimizer terminal solution retained",
    "Gradient check passed",
    "Maximum absolute gradient",
    "Gradient tolerance",
    "Post-estimation admissibility check",
    "Optimizer iterations",
    paste0(
      "lavaan free parameters (including observed residual variances)"
    ),
    "Negative observed residual variances",
    "Minimum observed residual variance",
    "Maximum absolute residual correlation",
    "Maximum absolute marginal residual z",
    "Correlation matrix minimum eigenvalue",
    "Correlation matrix condition number",
    "Sampling VCOV minimum eigenvalue",
    "Sampling VCOV condition number",
    "Sampling VCOV ridge added",
    "Captured lavaan warnings"
  ),
  value = c(
    analysis_status_short,
    as.character(converged),
    as.character(gradient_check_passed),
    format(maximum_absolute_gradient, digits = 8),
    format(gradient_tolerance, digits = 8),
    as.character(post_check),
    as.character(lavaan::lavInspect(fit, "iterations")),
    as.character(lavaan::lavInspect(fit, "npar")),
    as.character(negative_residual_variance_count),
    format(minimum_observed_residual_variance, digits = 8),
    format(max(abs(residual_matrix[upper.tri(residual_matrix)])), digits = 8),
    format(max(pair_diagnostics$absolute_marginal_residual_z), digits = 8),
    format(min(correlation_eigenvalues), digits = 8),
    format(max(correlation_eigenvalues) / min(correlation_eigenvalues), digits = 8),
    format(min(sampling_vcov_eigenvalues), digits = 8),
    format(
      max(sampling_vcov_eigenvalues) / min(sampling_vcov_eigenvalues),
      digits = 8
    ),
    "0",
    if (length(fit_warnings) == 0L) "None" else paste(fit_warnings, collapse = " | ")
  ),
  stringsAsFactors = FALSE
)

indicator_map <- data.frame(
  indicator = short_names,
  source_variable = long_names,
  construct = c(
    "Agreeableness", "Agreeableness",
    "Conscientiousness", "Conscientiousness",
    "Extraversion", "Extraversion",
    "Neuroticism", "Neuroticism",
    "Openness", "Openness",
    "IQ", "IQ"
  ),
  trait_factor = c(
    "agree", "agree",
    "consc", "consc",
    "extra", "extra",
    "neurot", "neurot",
    "open", "open",
    "iq", "iq"
  ),
  sample_factor = if (include_sample_two_factor) {
    c(
      "sample_one", "sample_two",
      "sample_one", "sample_two",
      "sample_one", "sample_two",
      "sample_one", "sample_two",
      "sample_one", "sample_two",
      "None", "None"
    )
  } else {
    c(
      "sample_one", "None",
      "sample_one", "None",
      "sample_one", "None",
      "sample_one", "None",
      "sample_one", "None",
      "None", "None"
    )
  },
  stringsAsFactors = FALSE
)

paper_results_mode <- identical(
  Sys.getenv("FULL_JACK_PAPER_RESULTS_MODE", unset = "0"),
  "1"
)
if (paper_results_mode) {
  paper_standardized_solution <- suppressWarnings(
    lavaan::standardizedSolution(
      fit,
      type = "std.all",
      se = TRUE,
      zstat = TRUE,
      pvalue = TRUE,
      ci = TRUE,
      level = 0.95,
      remove.eq = FALSE
    )
  )
  paper_standardized_solution <- paper_standardized_solution[
    paper_standardized_solution$op == "=~",
    c(
      "lhs", "rhs", "est.std", "se", "z", "pvalue",
      "ci.lower", "ci.upper"
    ),
    drop = FALSE
  ]

  paper_factor_labels <- c(
    gene_set_enrichment = general_factor_display_label,
    agree = "Agreeableness",
    consc = "Conscientiousness",
    extra = "Extraversion",
    neurot = "Neuroticism",
    open = "Openness",
    iq = "IQ",
    sample_one = "Sample-specific factor 1",
    sample_two = "Sample-specific factor 2"
  )
  paper_indicator_labels <- c(
    agree = "Agreeableness",
    consc = "Conscientiousness",
    extra = "Extraversion",
    neurot = "Neuroticism",
    open = "Openness",
    iq = "IQ",
    agr1 = "Agreeableness (half 1)",
    agr2 = "Agreeableness (half 2)",
    con1 = "Conscientiousness (half 1)",
    con2 = "Conscientiousness (half 2)",
    ext1 = "Extraversion (half 1)",
    ext2 = "Extraversion (half 2)",
    neu1 = "Neuroticism (half 1)",
    neu2 = "Neuroticism (half 2)",
    ope1 = "Openness (half 1)",
    ope2 = "Openness (half 2)",
    iq_female = "IQ (female)",
    iq_male = "IQ (male)"
  )
  equality_labels <- c(
    agree = "equal_agr",
    consc = "equal_con",
    extra = "equal_ext",
    neurot = "equal_neu",
    open = "equal_ope",
    iq = "equal_iq"
  )

  hierarchical_loadings_table <- data.frame(
    loading_type = ifelse(
      paper_standardized_solution$lhs == "gene_set_enrichment",
      "General factor to trait factor",
      ifelse(
        paper_standardized_solution$lhs %in% c(
          "sample_one",
          "sample_two"
        ),
        "Sample-specific factor to indicator",
        "Trait factor to split-half indicator"
      )
    ),
    factor = unname(
      paper_factor_labels[paper_standardized_solution$lhs]
    ),
    indicator = unname(
      paper_indicator_labels[paper_standardized_solution$rhs]
    ),
    standardized_loading = paper_standardized_solution$est.std,
    standard_error = paper_standardized_solution$se,
    ci_95_lower = paper_standardized_solution$ci.lower,
    ci_95_upper = paper_standardized_solution$ci.upper,
    z = paper_standardized_solution$z,
    p_value = paper_standardized_solution$pvalue,
    equality_constraint = ifelse(
      paper_standardized_solution$lhs %in% names(equality_labels),
      unname(equality_labels[paper_standardized_solution$lhs]),
      ""
    ),
    stringsAsFactors = FALSE
  )
  loading_type_order <- c(
    "General factor to trait factor",
    "Trait factor to split-half indicator",
    "Sample-specific factor to indicator"
  )
  hierarchical_loadings_table <- hierarchical_loadings_table[
    order(
      match(
        hierarchical_loadings_table$loading_type,
        loading_type_order
      ),
      hierarchical_loadings_table$factor,
      hierarchical_loadings_table$indicator
    ),
    ,
    drop = FALSE
  ]
  rownames(hierarchical_loadings_table) <- NULL

  implied_trait_pair_definitions <- data.frame(
    parameter = c(
      "corr_agree_consc",
      "corr_agree_extra",
      "corr_agree_neurot",
      "corr_agree_open",
      "corr_agree_iq",
      "corr_consc_extra",
      "corr_consc_neurot",
      "corr_consc_open",
      "corr_consc_iq",
      "corr_extra_neurot",
      "corr_extra_open",
      "corr_extra_iq",
      "corr_neurot_open",
      "corr_neurot_iq",
      "corr_open_iq"
    ),
    factor_1 = c(
      "agree", "agree", "agree", "agree", "agree",
      "consc", "consc", "consc", "consc",
      "extra", "extra", "extra",
      "neurot", "neurot",
      "open"
    ),
    factor_2 = c(
      "consc", "extra", "neurot", "open", "iq",
      "extra", "neurot", "open", "iq",
      "neurot", "open", "iq",
      "open", "iq",
      "iq"
    ),
    stringsAsFactors = FALSE
  )
  implied_latent_correlation_matrix <- lavaan::lavInspect(fit, "cor.lv")
  full_sample_implied_correlations <- vapply(
    seq_len(nrow(implied_trait_pair_definitions)),
    function(index) {
      implied_latent_correlation_matrix[
        implied_trait_pair_definitions$factor_1[[index]],
        implied_trait_pair_definitions$factor_2[[index]]
      ]
    },
    numeric(1)
  )
  names(full_sample_implied_correlations) <-
    implied_trait_pair_definitions$parameter
  if (any(!is.finite(full_sample_implied_correlations)) ||
      any(abs(full_sample_implied_correlations) > 1 + 1e-8)) {
    stop(
      "The full-sample model implies an invalid trait-factor correlation.",
      call. = FALSE
    )
  }
  full_sample_correlated_correlations <- vapply(
    seq_len(nrow(implied_trait_pair_definitions)),
    function(index) {
      correlated_decoded$trait_correlation[
        implied_trait_pair_definitions$factor_1[[index]],
        implied_trait_pair_definitions$factor_2[[index]]
      ]
    },
    numeric(1)
  )
  names(full_sample_correlated_correlations) <-
    implied_trait_pair_definitions$parameter
  if (any(!is.finite(full_sample_correlated_correlations)) ||
      any(abs(full_sample_correlated_correlations) > 1 + 1e-8)) {
    stop(
      paste0(
        "The full-sample correlated-factors model implies an invalid ",
        "trait-factor correlation."
      ),
      call. = FALSE
    )
  }

  hierarchical_defined_parameters <- lavaan::parameterEstimates(
    fit,
    standardized = TRUE,
    ci = FALSE
  )
  hierarchical_defined_parameters <- hierarchical_defined_parameters[
    hierarchical_defined_parameters$op == ":=",
    ,
    drop = FALSE
  ]
  definition_match <- match(
    implied_trait_pair_definitions$parameter,
    hierarchical_defined_parameters$lhs
  )
  if (anyNA(definition_match)) {
    stop(
      "The hierarchical model lacks one or more implied trait correlations.",
      call. = FALSE
    )
  }
  hierarchical_defined_parameters <- hierarchical_defined_parameters[
    definition_match,
    ,
    drop = FALSE
  ]
  standardized_defined_estimate <- hierarchical_defined_parameters$std.lv
  use_std_all <- !is.finite(standardized_defined_estimate)
  standardized_defined_estimate[use_std_all] <-
    hierarchical_defined_parameters$std.all[use_std_all]
  use_raw_estimate <- !is.finite(standardized_defined_estimate)
  standardized_defined_estimate[use_raw_estimate] <-
    hierarchical_defined_parameters$est[use_raw_estimate]
  if (
    any(!is.finite(standardized_defined_estimate)) ||
      max(abs(
        standardized_defined_estimate -
          full_sample_implied_correlations
      )) > 1e-8
  ) {
    stop(
      paste0(
        "The standardized defined parameters do not reproduce the ",
        "model-implied trait-factor correlations."
      ),
      call. = FALSE
    )
  }

  build_jackknife_correlation_matrix <- function(pair_values) {
    result <- diag(length(short_names))
    dimnames(result) <- list(short_names, short_names)
    result[pair_indices] <- pair_values
    result[pair_indices[, 2:1, drop = FALSE]] <- pair_values
    result
  }

  fit_jackknife_correlated_model <- function(
    replicate_correlation_matrix,
    block
  ) {
    replicate_pair_vector <- replicate_correlation_matrix[pair_indices]
    replicate_objective <- function(parameters) {
      implied_pair_vector <-
        decode_correlated_parameters(parameters)$implied_correlation[
          pair_indices
        ]
      residual <- replicate_pair_vector - implied_pair_vector
      as.numeric(crossprod(
        residual,
        wls_objective_precision %*% residual
      ))
    }

    captured_warnings <- character()
    optimizer <- tryCatch(
      withCallingHandlers(
        nlminb(
          start = correlated_parameters,
          objective = replicate_objective,
          lower = rep(-10, correlated_parameter_count),
          upper = rep(10, correlated_parameter_count),
          control = list(
            iter.max = 10000L,
            eval.max = 25000L,
            rel.tol = 1e-11,
            x.tol = 1e-9
          )
        ),
        warning = function(warning_condition) {
          captured_warnings <<- unique(c(
            captured_warnings,
            conditionMessage(warning_condition)
          ))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error_condition) error_condition
    )

    if (inherits(optimizer, "error")) {
      return(list(
        block = block,
        estimates = rep(
          NA_real_,
          nrow(implied_trait_pair_definitions)
        ),
        converged = FALSE,
        post_check = FALSE,
        gradient_check_passed = FALSE,
        maximum_absolute_gradient = NA_real_,
        optimizer_iterations = NA_integer_,
        admissible = FALSE,
        optimizer_message = conditionMessage(optimizer),
        warnings = if (length(captured_warnings) > 0L) {
          paste(captured_warnings, collapse = " | ")
        } else {
          "None"
        }
      ))
    }
    if (!isTRUE(optimizer$convergence == 0L)) {
      nlminb_optimizer <- optimizer
      fallback_optimizer <- tryCatch(
        withCallingHandlers(
          optim(
            par = nlminb_optimizer$par,
            fn = replicate_objective,
            method = "L-BFGS-B",
            lower = rep(-10, correlated_parameter_count),
            upper = rep(10, correlated_parameter_count),
            control = list(
              maxit = 10000L,
              factr = 1e7,
              pgtol = 1e-8
            )
          ),
          warning = function(warning_condition) {
            captured_warnings <<- unique(c(
              captured_warnings,
              conditionMessage(warning_condition)
            ))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(error_condition) error_condition
      )
      fallback_tolerance <- 1e-7 * max(
        1,
        abs(nlminb_optimizer$objective)
      )
      if (
        !inherits(fallback_optimizer, "error") &&
          isTRUE(fallback_optimizer$convergence == 0L) &&
          is.finite(fallback_optimizer$value) &&
          fallback_optimizer$value <=
            nlminb_optimizer$objective + fallback_tolerance
      ) {
        optimizer <- list(
          par = fallback_optimizer$par,
          objective = fallback_optimizer$value,
          convergence = 0L,
          iterations = nlminb_optimizer$iterations,
          message = paste0(
            "L-BFGS-B refinement after nlminb: ",
            nlminb_optimizer$message
          )
        )
      } else if (inherits(fallback_optimizer, "error")) {
        captured_warnings <- unique(c(
          captured_warnings,
          paste0(
            "L-BFGS-B refinement failed: ",
            conditionMessage(fallback_optimizer)
          )
        ))
      }
    }

    replicate_decoded <- decode_correlated_parameters(optimizer$par)
    replicate_estimates <- vapply(
      seq_len(nrow(implied_trait_pair_definitions)),
      function(index) {
        replicate_decoded$trait_correlation[
          implied_trait_pair_definitions$factor_1[[index]],
          implied_trait_pair_definitions$factor_2[[index]]
        ]
      },
      numeric(1)
    )
    maximum_absolute_gradient <- tryCatch(
      max(abs(numDeriv::grad(replicate_objective, optimizer$par))),
      error = function(error_condition) NA_real_
    )
    minimum_trait_eigenvalue <- min(eigen(
      replicate_decoded$trait_correlation,
      symmetric = TRUE,
      only.values = TRUE
    )$values)
    finite_estimates <- all(is.finite(replicate_estimates))
    bounded_estimates <- finite_estimates &&
      all(abs(replicate_estimates) <= 1 + 1e-8)
    replicate_converged <- isTRUE(optimizer$convergence == 0L)
    replicate_post_check <- (
      is.finite(minimum_trait_eigenvalue) &&
        minimum_trait_eigenvalue > 1e-6 &&
        all(is.finite(replicate_decoded$residual_variances)) &&
        min(replicate_decoded$residual_variances) > 1e-8 &&
        bounded_estimates
    )
    gradient_check_passed <- is.finite(maximum_absolute_gradient) &&
      maximum_absolute_gradient <= correlated_gradient_tolerance

    list(
      block = block,
      estimates = replicate_estimates,
      converged = replicate_converged,
      post_check = replicate_post_check,
      gradient_check_passed = gradient_check_passed,
      maximum_absolute_gradient = maximum_absolute_gradient,
      optimizer_iterations = if (
        length(optimizer$iterations) == 1L &&
          is.finite(optimizer$iterations)
      ) {
        as.integer(optimizer$iterations)
      } else {
        NA_integer_
      },
      admissible = (
        replicate_converged &&
          replicate_post_check &&
          gradient_check_passed
      ),
      optimizer_message = if (
        length(optimizer$message) == 1L &&
          nzchar(optimizer$message)
      ) {
        optimizer$message
      } else {
        "None"
      },
      warnings = if (length(captured_warnings) > 0L) {
        paste(captured_warnings, collapse = " | ")
      } else {
        "None"
      }
    )
  }

  fit_jackknife_hierarchical_model <- function(
    replicate_correlation_matrix,
    block
  ) {
    captured_warnings <- character()
    fit_result <- tryCatch(
      withCallingHandlers(
        lavaan::cfa(
          model = model_syntax,
          sample.cov = replicate_correlation_matrix,
          sample.nobs = jackknife_blocks,
          sample.cov.rescale = FALSE,
          std.lv = TRUE,
          estimator = "WLS",
          correlation = TRUE,
          WLS.V = wls_weight,
          NACOV = nacov,
          se = "none",
          test = "none",
          check.gradient = FALSE
        ),
        warning = function(warning_condition) {
          captured_warnings <<- unique(c(
            captured_warnings,
            conditionMessage(warning_condition)
          ))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error_condition) error_condition
    )

    if (inherits(fit_result, "error")) {
      return(list(
        block = block,
        estimates = rep(
          NA_real_,
          nrow(implied_trait_pair_definitions)
        ),
        converged = FALSE,
        post_check = FALSE,
        gradient_check_passed = FALSE,
        maximum_absolute_gradient = NA_real_,
        optimizer_iterations = NA_integer_,
        admissible = FALSE,
        warnings = paste(
          c(captured_warnings, conditionMessage(fit_result)),
          collapse = " | "
        )
      ))
    }

    replicate_converged <- isTRUE(
      lavaan::lavInspect(fit_result, "converged")
    )
    replicate_post_check <- withCallingHandlers(
      isTRUE(lavaan::lavInspect(fit_result, "post.check")),
      warning = function(warning_condition) {
        captured_warnings <<- unique(c(
          captured_warnings,
          conditionMessage(warning_condition)
        ))
        invokeRestart("muffleWarning")
      }
    )
    maximum_absolute_gradient <- if (
      length(fit_result@optim$dx) > 0L &&
        all(is.finite(fit_result@optim$dx))
    ) {
      max(abs(fit_result@optim$dx))
    } else {
      NA_real_
    }
    gradient_check_passed <- is.finite(maximum_absolute_gradient) &&
      maximum_absolute_gradient <= lavaan::lavOptions()$optim.dx.tol

    replicate_estimates <- rep(
      NA_real_,
      nrow(implied_trait_pair_definitions)
    )
    if (replicate_converged) {
      replicate_latent_correlations <- lavaan::lavInspect(
        fit_result,
        "cor.lv"
      )
      replicate_estimates <- vapply(
        seq_len(nrow(implied_trait_pair_definitions)),
        function(index) {
          replicate_latent_correlations[
            implied_trait_pair_definitions$factor_1[[index]],
            implied_trait_pair_definitions$factor_2[[index]]
          ]
        },
        numeric(1)
      )
    }
    finite_estimates <- all(is.finite(replicate_estimates))
    bounded_estimates <- finite_estimates &&
      all(abs(replicate_estimates) <= 1 + 1e-8)
    replicate_admissible <- (
      replicate_converged &&
        replicate_post_check &&
        gradient_check_passed &&
        bounded_estimates
    )

    list(
      block = block,
      estimates = replicate_estimates,
      converged = replicate_converged,
      post_check = replicate_post_check,
      gradient_check_passed = gradient_check_passed,
      maximum_absolute_gradient = maximum_absolute_gradient,
      optimizer_iterations = if (
        length(fit_result@optim$iterations) == 1L &&
          is.finite(fit_result@optim$iterations)
      ) {
        as.integer(fit_result@optim$iterations)
      } else {
        NA_integer_
      },
      admissible = replicate_admissible,
      warnings = if (length(captured_warnings) > 0L) {
        paste(captured_warnings, collapse = " | ")
      } else {
        "None"
      }
    )
  }

  jackknife_hierarchical_results <- vector(
    "list",
    jackknife_blocks
  )
  jackknife_correlated_results <- vector(
    "list",
    jackknife_blocks
  )
  for (block in seq_len(jackknife_blocks)) {
    replicate_correlation_matrix <-
      build_jackknife_correlation_matrix(
        jackknife_correlation_estimates[block, ]
      )
    if (min(eigen(
      replicate_correlation_matrix,
      symmetric = TRUE,
      only.values = TRUE
    )$values) <= 0) {
      stop(
        "Jackknife correlation matrix is not positive definite for block ",
        block,
        ".",
        call. = FALSE
      )
    }
    jackknife_correlated_results[[block]] <-
      fit_jackknife_correlated_model(
        replicate_correlation_matrix,
        block
      )
    jackknife_hierarchical_results[[block]] <-
      fit_jackknife_hierarchical_model(
        replicate_correlation_matrix,
        block
      )

    if (block %% 25L == 0L || block == jackknife_blocks) {
      message(
        "Completed ",
        analysis_output_label,
        " correlated and hierarchical SEM jackknife refit ",
        block,
        " of ",
        jackknife_blocks
      )
    }
  }

  jackknife_implied_correlation_matrix <- do.call(
    rbind,
    lapply(
      jackknife_hierarchical_results,
      `[[`,
      "estimates"
    )
  )
  colnames(jackknife_implied_correlation_matrix) <-
    implied_trait_pair_definitions$parameter
  jackknife_correlated_correlation_matrix <- do.call(
    rbind,
    lapply(
      jackknife_correlated_results,
      `[[`,
      "estimates"
    )
  )
  colnames(jackknife_correlated_correlation_matrix) <-
    implied_trait_pair_definitions$parameter
  jackknife_diagnostics <- data.frame(
    analysis = analysis_output_label,
    block = seq_len(jackknife_blocks),
    converged = vapply(
      jackknife_hierarchical_results,
      `[[`,
      logical(1),
      "converged"
    ),
    post_estimation_check = vapply(
      jackknife_hierarchical_results,
      `[[`,
      logical(1),
      "post_check"
    ),
    gradient_check_passed = vapply(
      jackknife_hierarchical_results,
      `[[`,
      logical(1),
      "gradient_check_passed"
    ),
    maximum_absolute_gradient = vapply(
      jackknife_hierarchical_results,
      `[[`,
      numeric(1),
      "maximum_absolute_gradient"
    ),
    optimizer_iterations = vapply(
      jackknife_hierarchical_results,
      `[[`,
      integer(1),
      "optimizer_iterations"
    ),
    admissible = vapply(
      jackknife_hierarchical_results,
      `[[`,
      logical(1),
      "admissible"
    ),
    captured_warnings = vapply(
      jackknife_hierarchical_results,
      `[[`,
      character(1),
      "warnings"
    ),
    stringsAsFactors = FALSE
  )
  jackknife_correlated_diagnostics <- data.frame(
    analysis = analysis_output_label,
    block = seq_len(jackknife_blocks),
    converged = vapply(
      jackknife_correlated_results,
      `[[`,
      logical(1),
      "converged"
    ),
    post_estimation_check = vapply(
      jackknife_correlated_results,
      `[[`,
      logical(1),
      "post_check"
    ),
    gradient_check_passed = vapply(
      jackknife_correlated_results,
      `[[`,
      logical(1),
      "gradient_check_passed"
    ),
    maximum_absolute_gradient = vapply(
      jackknife_correlated_results,
      `[[`,
      numeric(1),
      "maximum_absolute_gradient"
    ),
    optimizer_iterations = vapply(
      jackknife_correlated_results,
      `[[`,
      integer(1),
      "optimizer_iterations"
    ),
    admissible = vapply(
      jackknife_correlated_results,
      `[[`,
      logical(1),
      "admissible"
    ),
    optimizer_message = vapply(
      jackknife_correlated_results,
      `[[`,
      character(1),
      "optimizer_message"
    ),
    captured_warnings = vapply(
      jackknife_correlated_results,
      `[[`,
      character(1),
      "warnings"
    ),
    stringsAsFactors = FALSE
  )
  if (
    !all(jackknife_diagnostics$converged) ||
      any(!is.finite(jackknife_implied_correlation_matrix))
  ) {
    stop(
      paste0(
        "At least one hierarchical SEM jackknife refit failed to ",
        "converge or return finite implied correlations."
      ),
      call. = FALSE
    )
  }
  if (any(!is.finite(jackknife_correlated_correlation_matrix))) {
    stop(
      paste0(
        "At least one correlated-factors SEM jackknife refit failed to ",
        "return finite implied correlations."
      ),
      call. = FALSE
    )
  }

  jackknife_implied_correlation_means <- colMeans(
    jackknife_implied_correlation_matrix
  )
  centered_jackknife_implied_correlations <- sweep(
    jackknife_implied_correlation_matrix,
    MARGIN = 2L,
    STATS = jackknife_implied_correlation_means,
    FUN = "-"
  )
  jackknife_implied_standard_errors <- sqrt(
    (
      (jackknife_blocks - 1) / jackknife_blocks
    ) * colSums(centered_jackknife_implied_correlations^2)
  )
  jackknife_correlated_correlation_means <- colMeans(
    jackknife_correlated_correlation_matrix
  )
  centered_jackknife_correlated_correlations <- sweep(
    jackknife_correlated_correlation_matrix,
    MARGIN = 2L,
    STATS = jackknife_correlated_correlation_means,
    FUN = "-"
  )
  jackknife_correlated_standard_errors <- sqrt(
    (
      (jackknife_blocks - 1) / jackknife_blocks
    ) * colSums(centered_jackknife_correlated_correlations^2)
  )
  if (
    any(!is.finite(jackknife_implied_standard_errors)) ||
      any(jackknife_implied_standard_errors < 0) ||
      any(!is.finite(jackknife_correlated_standard_errors)) ||
      any(jackknife_correlated_standard_errors < 0)
  ) {
    stop(
      "The SEM jackknife produced invalid implied-correlation standard errors.",
      call. = FALSE
    )
  }

  hierarchical_implied_trait_correlations <- data.frame(
    analysis = analysis_output_label,
    trait_1 = unname(
      trait_factor_labels[implied_trait_pair_definitions$factor_1]
    ),
    trait_2 = unname(
      trait_factor_labels[implied_trait_pair_definitions$factor_2]
    ),
    correlation = full_sample_implied_correlations,
    standard_error = jackknife_implied_standard_errors,
    ci_95_lower = full_sample_implied_correlations -
      qnorm(0.975) * jackknife_implied_standard_errors,
    ci_95_upper = full_sample_implied_correlations +
      qnorm(0.975) * jackknife_implied_standard_errors,
    z = ifelse(
      jackknife_implied_standard_errors > 0,
      full_sample_implied_correlations /
        jackknife_implied_standard_errors,
      NA_real_
    ),
    p_value = ifelse(
      jackknife_implied_standard_errors > 0,
      2 * pnorm(
        abs(
          full_sample_implied_correlations /
            jackknife_implied_standard_errors
        ),
        lower.tail = FALSE
      ),
      NA_real_
    ),
    jackknife_blocks = jackknife_blocks,
    converged_replicates = sum(jackknife_diagnostics$converged),
    admissible_replicates = sum(jackknife_diagnostics$admissible),
    stringsAsFactors = FALSE
  )
  correlated_implied_trait_correlations <- data.frame(
    analysis = analysis_output_label,
    model = "Correlated-factors model",
    trait_1 = unname(
      trait_factor_labels[implied_trait_pair_definitions$factor_1]
    ),
    trait_2 = unname(
      trait_factor_labels[implied_trait_pair_definitions$factor_2]
    ),
    correlation = full_sample_correlated_correlations,
    standard_error = jackknife_correlated_standard_errors,
    ci_95_lower = full_sample_correlated_correlations -
      qnorm(0.975) * jackknife_correlated_standard_errors,
    ci_95_upper = full_sample_correlated_correlations +
      qnorm(0.975) * jackknife_correlated_standard_errors,
    z = ifelse(
      jackknife_correlated_standard_errors > 0,
      full_sample_correlated_correlations /
        jackknife_correlated_standard_errors,
      NA_real_
    ),
    p_value = ifelse(
      jackknife_correlated_standard_errors > 0,
      2 * pnorm(
        abs(
          full_sample_correlated_correlations /
            jackknife_correlated_standard_errors
        ),
        lower.tail = FALSE
      ),
      NA_real_
    ),
    jackknife_blocks = jackknife_blocks,
    converged_replicates = sum(
      jackknife_correlated_diagnostics$converged
    ),
    admissible_replicates = sum(
      jackknife_correlated_diagnostics$admissible
    ),
    stringsAsFactors = FALSE
  )
  model_implied_trait_correlations <- rbind(
    correlated_implied_trait_correlations,
    data.frame(
      analysis = hierarchical_implied_trait_correlations$analysis,
      model = "Hierarchical general-factor model",
      hierarchical_implied_trait_correlations[
        ,
        setdiff(
          names(hierarchical_implied_trait_correlations),
          "analysis"
        ),
        drop = FALSE
      ],
      stringsAsFactors = FALSE
    )
  )
  rownames(model_implied_trait_correlations) <- NULL

  jackknife_implied_correlation_output <- data.frame(
    analysis = analysis_output_label,
    block = seq_len(jackknife_blocks),
    jackknife_implied_correlation_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  jackknife_correlated_correlation_output <- data.frame(
    analysis = analysis_output_label,
    block = seq_len(jackknife_blocks),
    jackknife_correlated_correlation_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  model_implied_correlation_jackknife_output <- rbind(
    data.frame(
      analysis = jackknife_correlated_correlation_output$analysis,
      model = "Correlated-factors model",
      jackknife_correlated_correlation_output[
        ,
        setdiff(
          names(jackknife_correlated_correlation_output),
          "analysis"
        ),
        drop = FALSE
      ],
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    data.frame(
      analysis = jackknife_implied_correlation_output$analysis,
      model = "Hierarchical general-factor model",
      jackknife_implied_correlation_output[
        ,
        setdiff(
          names(jackknife_implied_correlation_output),
          "analysis"
        ),
        drop = FALSE
      ],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
  model_implied_correlation_jackknife_diagnostics <- rbind(
    data.frame(
      analysis = jackknife_correlated_diagnostics$analysis,
      model = "Correlated-factors model",
      jackknife_correlated_diagnostics[
        ,
        setdiff(
          names(jackknife_correlated_diagnostics),
          "analysis"
        ),
        drop = FALSE
      ],
      stringsAsFactors = FALSE
    ),
    data.frame(
      analysis = jackknife_diagnostics$analysis,
      model = "Hierarchical general-factor model",
      jackknife_diagnostics[
        ,
        setdiff(
          names(jackknife_diagnostics),
          "analysis"
        ),
        drop = FALSE
      ],
      optimizer_message = NA_character_,
      stringsAsFactors = FALSE
    )
  )
  diagnostic_column_order <- c(
    "analysis",
    "model",
    "block",
    "converged",
    "post_estimation_check",
    "gradient_check_passed",
    "maximum_absolute_gradient",
    "optimizer_iterations",
    "admissible",
    "optimizer_message",
    "captured_warnings"
  )
  model_implied_correlation_jackknife_diagnostics <-
    model_implied_correlation_jackknife_diagnostics[
      ,
      diagnostic_column_order,
      drop = FALSE
    ]

  write.csv(
    model_fit_table,
    model_fit_csv_path,
    row.names = FALSE,
    na = ""
  )
  write.csv(
    model_comparison_table,
    model_comparison_csv_path,
    row.names = FALSE,
    na = ""
  )
  write.csv(
    model_diagnostics_table,
    file.path(
      output_directory,
      paste0(artifact_base_name, "_model_diagnostics.csv")
    ),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    hierarchical_loadings_table,
    file.path(
      output_directory,
      paste0(artifact_base_name, "_hierarchical_loadings.csv")
    ),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    hierarchical_implied_trait_correlations,
    file.path(
      output_directory,
      paste0(
        artifact_base_name,
        "_hierarchical_implied_trait_correlations.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    model_implied_trait_correlations,
    file.path(
      output_directory,
      paste0(
        artifact_base_name,
        "_implied_trait_correlations.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    jackknife_implied_correlation_output,
    file.path(
      output_directory,
      paste0(
        artifact_base_name,
        "_hierarchical_implied_trait_correlation_jackknife.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    jackknife_diagnostics,
    file.path(
      output_directory,
      paste0(
        artifact_base_name,
        "_hierarchical_implied_trait_correlation_",
        "jackknife_diagnostics.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    model_implied_correlation_jackknife_output,
    file.path(
      output_directory,
      paste0(
        artifact_base_name,
        "_implied_trait_correlation_jackknife.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    model_implied_correlation_jackknife_diagnostics,
    file.path(
      output_directory,
      paste0(
        artifact_base_name,
        "_implied_trait_correlation_jackknife_diagnostics.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )

  message("Paper-results SEM tables written to: ", output_directory)
  quit(save = "no", status = 0L, runLast = FALSE)
}
