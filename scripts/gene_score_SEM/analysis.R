#!/usr/bin/env Rscript

# ==============================================================================
# 1. Setup and analysis configuration
# ==============================================================================

required_packages <- c(
  "data.table",
  "jsonlite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) == 1L) {
  script_path <- normalizePath(
    sub("^--file=", "", script_argument),
    mustWork = TRUE
  )
} else {
  script_path <- normalizePath("scripts/analysis.R", mustWork = TRUE)
}

project_directory <- dirname(dirname(script_path))
jackknife_directory <- file.path(
  project_directory,
  "data",
  "split_half_output",
  "output"
)
full_sample_archive <- file.path(
  project_directory,
  "data",
  "split_half_geneset_output.tar.gz"
)
output_directory <- file.path(project_directory, "output")
correlation_directory <- file.path(
  output_directory,
  "correlation_matrices"
)
regression_directory <- file.path(output_directory, "regression")
sampling_directory <- file.path(
  output_directory,
  "sampling_covariances"
)
gene_set_output_directory <- file.path(
  output_directory,
  "gene_sets"
)

excluded_collections <- "neural_brain"
excluded_gene_sets <- paste0(
  "GOBP_NEGATIVE_REGULATION_OF_INTRACELLULAR_",
  "LIPID_TRANSPORT"
)

trait_names <- c(
  "iq_female_dir",
  "iq_male_dir",
  "ReGPC_agr_half_one_no23_dir",
  "ReGPC_agr_half_two_no23_dir",
  "ReGPC_con_half_one_no23_dir",
  "ReGPC_con_half_two_no23_dir",
  "ReGPC_ext_half_one_no23_dir",
  "ReGPC_ext_half_two_no23_dir",
  "ReGPC_neu_half_one_no23_dir",
  "ReGPC_neu_half_two_no23_dir",
  "ReGPC_ope_half_one_no23_dir",
  "ReGPC_ope_half_two_no23_dir"
)

trait_halves <- data.table(
  trait = c("agree", "consc", "extra", "neurot", "open", "iq"),
  half_1 = c(
    "ReGPC_agr_half_one_no23_dir",
    "ReGPC_con_half_one_no23_dir",
    "ReGPC_ext_half_one_no23_dir",
    "ReGPC_neu_half_one_no23_dir",
    "ReGPC_ope_half_one_no23_dir",
    "iq_female_dir"
  ),
  half_2 = c(
    "ReGPC_agr_half_two_no23_dir",
    "ReGPC_con_half_two_no23_dir",
    "ReGPC_ext_half_two_no23_dir",
    "ReGPC_neu_half_two_no23_dir",
    "ReGPC_ope_half_two_no23_dir",
    "iq_male_dir"
  )
)

jackknife_blocks <- 200L
combined_result_name <- "all_gene_sets"
duplicate_tolerance <- 1e-12

if (!dir.exists(jackknife_directory)) {
  stop("Jackknife input directory not found: ", jackknife_directory)
}
if (!file.exists(full_sample_archive)) {
  stop("Full-sample archive not found: ", full_sample_archive)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(correlation_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(regression_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(sampling_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(
  gene_set_output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Helper functions
# ==============================================================================

identify_gene_set_column <- function(path) {
  column_names <- names(
    fread(
      path,
      skip = "VARIABLE",
      nrows = 0L,
      showProgress = FALSE
    )
  )

  if ("FULL_NAME" %in% column_names) {
    return("FULL_NAME")
  }
  if ("VARIABLE" %in% column_names) {
    return("VARIABLE")
  }

  stop("No gene-set identifier column found in: ", path)
}


read_gene_set_results <- function(path, identifier_column) {
  result <- fread(
    path,
    skip = "VARIABLE",
    select = c(identifier_column, "BETA_STD"),
    showProgress = FALSE
  )
  setnames(result, identifier_column, "gene_set")

  result <- result[!gene_set %chin% excluded_gene_sets]

  if (anyDuplicated(result$gene_set)) {
    stop("Duplicated gene-set identifiers within file: ", path)
  }
  if (anyNA(result$gene_set) || anyNA(result$BETA_STD)) {
    stop("Missing gene-set identifiers or BETA_STD values in: ", path)
  }
  if (!all(is.finite(result$BETA_STD))) {
    stop("Non-finite BETA_STD values in: ", path)
  }

  result
}


build_trait_matrix <- function(files, identifier_column, context) {
  if (!identical(names(files), trait_names)) {
    stop("Trait files are not in the configured order for: ", context)
  }

  trait_results <- lapply(
    files,
    read_gene_set_results,
    identifier_column = identifier_column
  )

  common_gene_sets <- Reduce(
    intersect,
    lapply(trait_results, `[[`, "gene_set")
  )

  if (length(common_gene_sets) < 2L) {
    stop("Fewer than two common gene sets for: ", context)
  }

  beta_matrix <- vapply(
    trait_results,
    function(result) {
      result$BETA_STD[match(common_gene_sets, result$gene_set)]
    },
    numeric(length(common_gene_sets))
  )
  colnames(beta_matrix) <- trait_names

  if (anyNA(beta_matrix) || !all(is.finite(beta_matrix))) {
    stop("Failed to align finite BETA_STD values for: ", context)
  }

  list(
    gene_set = common_gene_sets,
    beta_std = beta_matrix
  )
}


combine_gene_set_matrices <- function(group_results, context) {
  pooled_gene_sets <- unlist(
    lapply(group_results, `[[`, "gene_set"),
    use.names = FALSE
  )
  pooled_betas <- do.call(
    rbind,
    lapply(group_results, `[[`, "beta_std")
  )

  duplicated_names <- unique(
    pooled_gene_sets[
      duplicated(pooled_gene_sets) |
        duplicated(pooled_gene_sets, fromLast = TRUE)
    ]
  )

  if (length(duplicated_names) > 0L) {
    for (gene_set in duplicated_names) {
      duplicate_rows <- pooled_betas[pooled_gene_sets == gene_set, , drop = FALSE]
      value_ranges <- apply(duplicate_rows, 2L, function(x) max(x) - min(x))

      if (any(value_ranges > duplicate_tolerance)) {
        stop(
          "Conflicting duplicate gene set after exclusions in ",
          context,
          ": ",
          gene_set
        )
      }
    }

    keep <- !duplicated(pooled_gene_sets)
    pooled_gene_sets <- pooled_gene_sets[keep]
    pooled_betas <- pooled_betas[keep, , drop = FALSE]
  }

  list(
    gene_set = pooled_gene_sets,
    beta_std = pooled_betas
  )
}


estimate_correlation_matrix <- function(group_result, context) {
  correlation_matrix <- cor(
    group_result$beta_std,
    method = "pearson",
    use = "everything"
  )

  if (!all(is.finite(correlation_matrix))) {
    stop("Non-finite correlation estimate for: ", context)
  }

  correlation_matrix
}


write_named_matrix <- function(matrix_object, path) {
  write.csv(
    matrix_object,
    file = path,
    row.names = TRUE,
    quote = TRUE
  )
}


validate_symmetric_matrix <- function(matrix_object, context) {
  if (!all(is.finite(matrix_object))) {
    stop("Non-finite matrix values for: ", context)
  }
  if (!isTRUE(all.equal(matrix_object, t(matrix_object), tolerance = 1e-10))) {
    stop("Matrix is not symmetric for: ", context)
  }
  invisible(TRUE)
}


# ==============================================================================
# 3. Discover collections and index jackknife files
# ==============================================================================

all_collections <- sort(
  basename(
    list.dirs(
      jackknife_directory,
      recursive = FALSE,
      full.names = TRUE
    )
  )
)
collections <- setdiff(all_collections, excluded_collections)

if (length(collections) == 0L) {
  stop("No gene-set collections remain after exclusions.")
}

message(
  "Including ",
  length(collections),
  " collections; excluding: ",
  paste(excluded_collections, collapse = ", ")
)

jackknife_index <- setNames(vector("list", length(collections)), collections)

for (collection in collections) {
  collection_index <- setNames(vector("list", length(trait_names)), trait_names)

  for (trait in trait_names) {
    trait_directory <- file.path(
      jackknife_directory,
      collection,
      paste0("jack_", trait)
    )
    files <- list.files(
      trait_directory,
      pattern = "_jk[0-9]{3}[.]gsa[.]out$",
      full.names = TRUE
    )

    block_ids <- sub(
      "^.*_jk([0-9]{3})[.]gsa[.]out$",
      "\\1",
      files
    )
    expected_ids <- sprintf("%03d", seq_len(jackknife_blocks))

    if (
      length(files) != jackknife_blocks ||
        anyDuplicated(block_ids) ||
        !setequal(block_ids, expected_ids)
    ) {
      stop(
        "Expected exactly ",
        jackknife_blocks,
        " jackknife results for ",
        collection,
        " / ",
        trait
      )
    }

    collection_index[[trait]] <- setNames(files, block_ids)[expected_ids]
  }

  jackknife_index[[collection]] <- collection_index
}


# ==============================================================================
# 4. Extract and index full-sample results
# ==============================================================================

full_sample_temp_directory <- tempfile("full_jack_full_sample_")
dir.create(full_sample_temp_directory)

full_sample_members <- untar(full_sample_archive, list = TRUE)
full_sample_result_members <- full_sample_members[
  grepl("[.]gsa[.]out$", full_sample_members)
]

untar_status <- untar(
  full_sample_archive,
  files = full_sample_result_members,
  exdir = full_sample_temp_directory
)
if (!identical(untar_status, 0L)) {
  stop("Failed to extract the full-sample archive.")
}

full_sample_directory <- file.path(
  full_sample_temp_directory,
  "geneset_output"
)
if (!dir.exists(full_sample_directory)) {
  stop("Expected geneset_output directory was not found in the archive.")
}

identifier_columns <- setNames(character(length(collections)), collections)
full_sample_index <- setNames(vector("list", length(collections)), collections)

for (collection in collections) {
  files <- file.path(
    full_sample_directory,
    collection,
    paste0(trait_names, ".gsa.out")
  )
  names(files) <- trait_names

  if (!all(file.exists(files))) {
    stop("Missing full-sample results for collection: ", collection)
  }

  identifier_columns[[collection]] <- identify_gene_set_column(files[[1L]])
  full_sample_index[[collection]] <- files
}


# ==============================================================================
# 5. Estimate full-sample point correlation matrices
# ==============================================================================

result_names <- c(collections, combined_result_name)
point_correlations <- setNames(
  vector("list", length(result_names)),
  result_names
)
n_gene_sets_by_result <- setNames(
  integer(length(result_names)),
  result_names
)
full_sample_group_results <- setNames(
  vector("list", length(collections)),
  collections
)

for (collection in collections) {
  group_result <- build_trait_matrix(
    full_sample_index[[collection]],
    identifier_columns[[collection]],
    paste0("full sample / ", collection)
  )

  full_sample_group_results[[collection]] <- group_result
  n_gene_sets_by_result[[collection]] <- length(group_result$gene_set)
  point_correlations[[collection]] <- estimate_correlation_matrix(
    group_result,
    paste0("full sample / ", collection)
  )
}

combined_full_sample <- combine_gene_set_matrices(
  full_sample_group_results,
  "full sample / all gene sets"
)
n_gene_sets_by_result[[combined_result_name]] <-
  length(combined_full_sample$gene_set)
point_correlations[[combined_result_name]] <- estimate_correlation_matrix(
  combined_full_sample,
  "full sample / all gene sets"
)

rm(full_sample_group_results, combined_full_sample)
unlink(full_sample_temp_directory, recursive = TRUE, force = TRUE)
invisible(gc())


# ==============================================================================
# 6. Estimate correlations in each jackknife replicate
# ==============================================================================

correlation_pair_indices <- t(combn(seq_along(trait_names), 2L))
correlation_pair_names <- paste(
  trait_names[correlation_pair_indices[, 1L]],
  trait_names[correlation_pair_indices[, 2L]],
  sep = "__"
)

jackknife_correlations <- setNames(
  lapply(
    result_names,
    function(x) {
      matrix(
        NA_real_,
        nrow = jackknife_blocks,
        ncol = length(correlation_pair_names),
        dimnames = list(
          sprintf("jk%03d", seq_len(jackknife_blocks)),
          correlation_pair_names
        )
      )
    }
  ),
  result_names
)

for (block in seq_len(jackknife_blocks)) {
  block_id <- sprintf("%03d", block)
  block_group_results <- setNames(
    vector("list", length(collections)),
    collections
  )

  for (collection in collections) {
    files <- vapply(
      trait_names,
      function(trait) {
        jackknife_index[[collection]][[trait]][[block_id]]
      },
      character(1)
    )
    names(files) <- trait_names

    group_result <- build_trait_matrix(
      files,
      identifier_columns[[collection]],
      paste0("jackknife ", block_id, " / ", collection)
    )
    block_group_results[[collection]] <- group_result

    correlation_matrix <- estimate_correlation_matrix(
      group_result,
      paste0("jackknife ", block_id, " / ", collection)
    )
    jackknife_correlations[[collection]][block, ] <-
      correlation_matrix[correlation_pair_indices]
  }

  combined_block <- combine_gene_set_matrices(
    block_group_results,
    paste0("jackknife ", block_id, " / all gene sets")
  )
  combined_correlation <- estimate_correlation_matrix(
    combined_block,
    paste0("jackknife ", block_id, " / all gene sets")
  )
  jackknife_correlations[[combined_result_name]][block, ] <-
    combined_correlation[correlation_pair_indices]

  if (block %% 10L == 0L || block == jackknife_blocks) {
    message(
      "Completed jackknife block ",
      block,
      " of ",
      jackknife_blocks
    )
  }
}


# ==============================================================================
# 7. Calculate jackknife sampling variance-covariance matrices
# ==============================================================================

sampling_vcov <- lapply(
  jackknife_correlations,
  function(replicate_estimates) {
    replicate_means <- colMeans(replicate_estimates)
    centered_estimates <- sweep(
      replicate_estimates,
      MARGIN = 2L,
      STATS = replicate_means,
      FUN = "-"
    )

    covariance_matrix <- (
      (jackknife_blocks - 1) / jackknife_blocks
    ) * crossprod(centered_estimates)

    dimnames(covariance_matrix) <- list(
      correlation_pair_names,
      correlation_pair_names
    )
    covariance_matrix
  }
)


# ==============================================================================
# 8. Fit the collection-level regression in every jackknife replicate
# ==============================================================================

regression_collections <- collections
trait_pair_indices <- t(combn(seq_len(nrow(trait_halves)), 2L))
iq_trait_index <- match("iq", trait_halves$trait)
iq_pair_indicator <- apply(
  trait_pair_indices,
  MARGIN = 1L,
  FUN = function(indices) iq_trait_index %in% indices
)

pair_column_name <- function(trait_1, trait_2) {
  index_1 <- match(trait_1, trait_names)
  index_2 <- match(trait_2, trait_names)

  if (is.na(index_1) || is.na(index_2) || index_1 == index_2) {
    stop("Invalid trait pair: ", trait_1, " / ", trait_2)
  }

  ordered_indices <- sort(c(index_1, index_2))
  paste(trait_names[ordered_indices], collapse = "__")
}


summarize_point_correlation <- function(correlation_matrix) {
  reliability <- vapply(
    seq_len(nrow(trait_halves)),
    function(i) {
      correlation_matrix[
        trait_halves$half_1[[i]],
        trait_halves$half_2[[i]]
      ]
    },
    numeric(1)
  )

  cross_trait <- vapply(
    seq_len(nrow(trait_pair_indices)),
    function(i) {
      trait_1 <- trait_pair_indices[i, 1L]
      trait_2 <- trait_pair_indices[i, 2L]
      mean(c(
        correlation_matrix[
          trait_halves$half_1[[trait_1]],
          trait_halves$half_2[[trait_2]]
        ],
        correlation_matrix[
          trait_halves$half_2[[trait_1]],
          trait_halves$half_1[[trait_2]]
        ]
      ))
    },
    numeric(1)
  )

  c(
    avg_reliability = mean(reliability),
    avg_trait_correlation = mean(cross_trait),
    avg_iq_correlation = mean(cross_trait[iq_pair_indicator]),
    avg_big5_correlation = mean(cross_trait[!iq_pair_indicator]),
    min_reliability = min(reliability),
    max_reliability = max(reliability)
  )
}


reliability_columns <- vapply(
  seq_len(nrow(trait_halves)),
  function(i) {
    pair_column_name(
      trait_halves$half_1[[i]],
      trait_halves$half_2[[i]]
    )
  },
  character(1)
)

cross_trait_columns <- lapply(
  seq_len(nrow(trait_pair_indices)),
  function(i) {
    trait_1 <- trait_pair_indices[i, 1L]
    trait_2 <- trait_pair_indices[i, 2L]
    c(
      pair_column_name(
        trait_halves$half_1[[trait_1]],
        trait_halves$half_2[[trait_2]]
      ),
      pair_column_name(
        trait_halves$half_2[[trait_1]],
        trait_halves$half_1[[trait_2]]
      )
    )
  }
)

collection_regression_full_sample <- rbindlist(
  lapply(
    regression_collections,
    function(collection) {
      summary_values <- summarize_point_correlation(
        point_correlations[[collection]]
      )

      data.table(
        collection = collection,
        avg_reliability = unname(summary_values[["avg_reliability"]]),
        avg_trait_correlation = unname(
          summary_values[["avg_trait_correlation"]]
        ),
        avg_iq_correlation = unname(
          summary_values[["avg_iq_correlation"]]
        ),
        avg_big5_correlation = unname(
          summary_values[["avg_big5_correlation"]]
        ),
        min_reliability = unname(summary_values[["min_reliability"]]),
        max_reliability = unname(summary_values[["max_reliability"]])
      )
    }
  )
)

collection_regression_jackknife <- rbindlist(
  lapply(
    regression_collections,
    function(collection) {
      replicate_estimates <- jackknife_correlations[[collection]]

      reliability <- replicate_estimates[
        ,
        reliability_columns,
        drop = FALSE
      ]
      cross_trait <- vapply(
        cross_trait_columns,
        function(columns) {
          rowMeans(replicate_estimates[, columns, drop = FALSE])
        },
        numeric(nrow(replicate_estimates))
      )

      data.table(
        block = seq_len(nrow(replicate_estimates)),
        collection = collection,
        avg_reliability = rowMeans(reliability),
        avg_trait_correlation = rowMeans(cross_trait),
        avg_iq_correlation = rowMeans(
          cross_trait[, iq_pair_indicator, drop = FALSE]
        ),
        avg_big5_correlation = rowMeans(
          cross_trait[, !iq_pair_indicator, drop = FALSE]
        )
      )
    }
  )
)

if (
  any(!is.finite(
    collection_regression_full_sample$avg_reliability
  )) ||
    any(!is.finite(
      collection_regression_full_sample$avg_trait_correlation
    )) ||
    any(!is.finite(
      collection_regression_jackknife$avg_reliability
    )) ||
    any(!is.finite(
      collection_regression_jackknife$avg_trait_correlation
    )) ||
    any(!is.finite(
      collection_regression_full_sample$avg_iq_correlation
    )) ||
    any(!is.finite(
      collection_regression_full_sample$avg_big5_correlation
    )) ||
    any(!is.finite(
      collection_regression_jackknife$avg_iq_correlation
    )) ||
    any(!is.finite(
      collection_regression_jackknife$avg_big5_correlation
    ))
) {
  stop("Non-finite collection-level regression inputs.")
}

collections_per_block <- collection_regression_jackknife[
  ,
  .N,
  by = block
]
if (
  nrow(collections_per_block) != jackknife_blocks ||
    any(collections_per_block$N != length(regression_collections))
) {
  stop("Jackknife regression blocks do not contain every collection.")
}

full_regression_fit <- lm(
  avg_trait_correlation ~ avg_reliability,
  data = collection_regression_full_sample
)
full_regression_coefficients <- c(
  intercept = unname(coef(full_regression_fit)[[1L]]),
  slope = unname(coef(full_regression_fit)[[2L]])
)

jackknife_regression <- collection_regression_jackknife[
  ,
  {
    fit <- lm(
      avg_trait_correlation ~ avg_reliability,
      data = .SD
    )
    list(
      intercept = unname(coef(fit)[[1L]]),
      slope = unname(coef(fit)[[2L]]),
      r_squared = summary(fit)$r.squared,
      residual_sd = summary(fit)$sigma,
      n_collections = .N
    )
  },
  by = block
]

jackknife_coefficient_estimates <- as.matrix(
  jackknife_regression[, .(intercept, slope)]
)
jackknife_coefficient_means <- colMeans(jackknife_coefficient_estimates)
centered_coefficient_estimates <- sweep(
  jackknife_coefficient_estimates,
  MARGIN = 2L,
  STATS = jackknife_coefficient_means,
  FUN = "-"
)
regression_coefficient_vcov <- (
  (jackknife_blocks - 1) / jackknife_blocks
) * crossprod(centered_coefficient_estimates)
dimnames(regression_coefficient_vcov) <- list(
  c("intercept", "slope"),
  c("intercept", "slope")
)

normal_critical_value <- qnorm(0.975)
regression_coefficient_se <- sqrt(diag(regression_coefficient_vcov))
regression_coefficient_summary <- data.table(
  term = c("intercept", "slope"),
  full_sample_estimate = unname(
    full_regression_coefficients[c("intercept", "slope")]
  ),
  mean_leave_one_out_estimate = unname(
    jackknife_coefficient_means[c("intercept", "slope")]
  ),
  jackknife_se = unname(
    regression_coefficient_se[c("intercept", "slope")]
  )
)
regression_coefficient_summary[
  ,
  `:=`(
    normal_95_low = full_sample_estimate -
      normal_critical_value * jackknife_se,
    normal_95_high = full_sample_estimate +
      normal_critical_value * jackknife_se
  )
]

jackknife_standard_error <- function(x) {
  centered <- x - mean(x)
  sqrt(((length(x) - 1) / length(x)) * sum(centered^2))
}

collection_point_uncertainty <- collection_regression_jackknife[
  ,
  .(
    se_avg_reliability = jackknife_standard_error(avg_reliability),
    se_avg_trait_correlation = jackknife_standard_error(
      avg_trait_correlation
    ),
    se_avg_iq_correlation = jackknife_standard_error(
      avg_iq_correlation
    ),
    se_avg_big5_correlation = jackknife_standard_error(
      avg_big5_correlation
    ),
    sampling_covariance = (
      (.N - 1) / .N
    ) * sum(
      (avg_reliability - mean(avg_reliability)) *
        (avg_trait_correlation - mean(avg_trait_correlation))
    ),
    sampling_covariance_reliability_iq = (
      (.N - 1) / .N
    ) * sum(
      (avg_reliability - mean(avg_reliability)) *
        (avg_iq_correlation - mean(avg_iq_correlation))
    ),
    sampling_covariance_reliability_big5 = (
      (.N - 1) / .N
    ) * sum(
      (avg_reliability - mean(avg_reliability)) *
        (avg_big5_correlation - mean(avg_big5_correlation))
    )
  ),
  by = collection
]

collection_plot_data <- merge(
  collection_regression_full_sample,
  collection_point_uncertainty,
  by = "collection",
  all.x = TRUE,
  sort = FALSE
)
collection_plot_data[
  ,
  collection_order := match(collection, regression_collections)
]
setorder(collection_plot_data, collection_order)
collection_plot_data[, collection_order := NULL]
collection_plot_data[
  ,
  `:=`(
    reliability_95_low = avg_reliability -
      normal_critical_value * se_avg_reliability,
    reliability_95_high = avg_reliability +
      normal_critical_value * se_avg_reliability,
    trait_correlation_95_low = avg_trait_correlation -
      normal_critical_value * se_avg_trait_correlation,
    trait_correlation_95_high = avg_trait_correlation +
      normal_critical_value * se_avg_trait_correlation,
    iq_correlation_95_low = avg_iq_correlation -
      normal_critical_value * se_avg_iq_correlation,
    iq_correlation_95_high = avg_iq_correlation +
      normal_critical_value * se_avg_iq_correlation,
    big5_correlation_95_low = avg_big5_correlation -
      normal_critical_value * se_avg_big5_correlation,
    big5_correlation_95_high = avg_big5_correlation +
      normal_critical_value * se_avg_big5_correlation
  )
]

regression_band <- data.table(
  avg_reliability = seq(
    min(collection_plot_data$avg_reliability),
    max(collection_plot_data$avg_reliability),
    length.out = 300L
  )
)
band_design_matrix <- cbind(
  intercept = 1,
  slope = regression_band$avg_reliability
)
regression_band[
  ,
  `:=`(
    fitted_value = full_regression_coefficients[["intercept"]] +
      full_regression_coefficients[["slope"]] * avg_reliability,
    fitted_se = sqrt(
      pmax(
        rowSums(
          (band_design_matrix %*% regression_coefficient_vcov) *
            band_design_matrix
        ),
        0
      )
    )
  )
]

two_series_ids <- c("iq_big5", "big5_big5")
two_series_labels <- c(
  iq_big5 = "IQ - Big Five",
  big5_big5 = "Big Five - Big Five"
)

two_series_full_sample <- rbindlist(list(
  collection_regression_full_sample[
    ,
    .(
      collection,
      avg_reliability,
      correlation = avg_iq_correlation,
      series = "iq_big5"
    )
  ],
  collection_regression_full_sample[
    ,
    .(
      collection,
      avg_reliability,
      correlation = avg_big5_correlation,
      series = "big5_big5"
    )
  ]
))

two_series_jackknife <- rbindlist(list(
  collection_regression_jackknife[
    ,
    .(
      block,
      collection,
      avg_reliability,
      correlation = avg_iq_correlation,
      series = "iq_big5"
    )
  ],
  collection_regression_jackknife[
    ,
    .(
      block,
      collection,
      avg_reliability,
      correlation = avg_big5_correlation,
      series = "big5_big5"
    )
  ]
))
two_series_full_sample[
  ,
  series_label := unname(two_series_labels[series])
]
two_series_jackknife[
  ,
  series_label := unname(two_series_labels[series])
]

two_series_full_fits <- setNames(
  lapply(
    two_series_ids,
    function(series_id) {
      lm(
        correlation ~ avg_reliability,
        data = two_series_full_sample[series == series_id]
      )
    }
  ),
  two_series_ids
)

two_series_jackknife_regression <- two_series_jackknife[
  ,
  {
    fit <- lm(correlation ~ avg_reliability, data = .SD)
    list(
      intercept = unname(coef(fit)[[1L]]),
      slope = unname(coef(fit)[[2L]]),
      r_squared = summary(fit)$r.squared,
      residual_sd = summary(fit)$sigma,
      n_collections = .N
    )
  },
  by = .(series, block)
]

two_series_coefficient_vcov <- setNames(
  lapply(
    two_series_ids,
    function(series_id) {
      coefficient_estimates <- as.matrix(
        two_series_jackknife_regression[
          series == series_id,
          .(intercept, slope)
        ]
      )
      coefficient_means <- colMeans(coefficient_estimates)
      centered_estimates <- sweep(
        coefficient_estimates,
        MARGIN = 2L,
        STATS = coefficient_means,
        FUN = "-"
      )
      covariance_matrix <- (
        (jackknife_blocks - 1) / jackknife_blocks
      ) * crossprod(centered_estimates)
      dimnames(covariance_matrix) <- list(
        c("intercept", "slope"),
        c("intercept", "slope")
      )
      covariance_matrix
    }
  ),
  two_series_ids
)

two_series_coefficient_summary <- rbindlist(
  lapply(
    two_series_ids,
    function(series_id) {
      full_fit <- two_series_full_fits[[series_id]]
      full_coefficients <- c(
        intercept = unname(coef(full_fit)[[1L]]),
        slope = unname(coef(full_fit)[[2L]])
      )
      covariance_matrix <- two_series_coefficient_vcov[[series_id]]
      coefficient_se <- sqrt(diag(covariance_matrix))

      data.table(
        series = series_id,
        series_label = unname(two_series_labels[[series_id]]),
        term = c("intercept", "slope"),
        full_sample_estimate = unname(
          full_coefficients[c("intercept", "slope")]
        ),
        jackknife_se = unname(
          coefficient_se[c("intercept", "slope")]
        ),
        normal_95_low = unname(
          full_coefficients[c("intercept", "slope")] -
            normal_critical_value *
              coefficient_se[c("intercept", "slope")]
        ),
        normal_95_high = unname(
          full_coefficients[c("intercept", "slope")] +
            normal_critical_value *
              coefficient_se[c("intercept", "slope")]
        ),
        full_sample_r_squared = summary(full_fit)$r.squared
      )
    }
  )
)

two_series_regression_band <- rbindlist(
  lapply(
    two_series_ids,
    function(series_id) {
      full_fit <- two_series_full_fits[[series_id]]
      full_coefficients <- c(
        intercept = unname(coef(full_fit)[[1L]]),
        slope = unname(coef(full_fit)[[2L]])
      )
      covariance_matrix <- two_series_coefficient_vcov[[series_id]]
      reliability_values <- seq(
        min(two_series_full_sample$avg_reliability),
        max(two_series_full_sample$avg_reliability),
        length.out = 300L
      )
      design_matrix <- cbind(
        intercept = 1,
        slope = reliability_values
      )
      fitted_se <- sqrt(
        pmax(
          rowSums(
            (design_matrix %*% covariance_matrix) * design_matrix
          ),
          0
        )
      )
      fitted_value <- full_coefficients[["intercept"]] +
        full_coefficients[["slope"]] * reliability_values

      data.table(
        series = series_id,
        series_label = unname(two_series_labels[[series_id]]),
        avg_reliability = reliability_values,
        fitted_value = fitted_value,
        fitted_se = fitted_se,
        normal_95_low = fitted_value -
          normal_critical_value * fitted_se,
        normal_95_high = fitted_value +
          normal_critical_value * fitted_se
      )
    }
  )
)
regression_band[
  ,
  `:=`(
    normal_95_low = fitted_value -
      normal_critical_value * fitted_se,
    normal_95_high = fitted_value +
      normal_critical_value * fitted_se
  )
]


# ==============================================================================
# 9. Validate and write CSV outputs
# ==============================================================================

for (result_name in result_names) {
  correlation_matrix <- point_correlations[[result_name]]
  covariance_matrix <- sampling_vcov[[result_name]]

  validate_symmetric_matrix(
    correlation_matrix,
    paste0(result_name, " point correlation")
  )
  validate_symmetric_matrix(
    covariance_matrix,
    paste0(result_name, " sampling covariance")
  )

  if (!isTRUE(all.equal(
    unname(diag(correlation_matrix)),
    rep(1, length(trait_names))
  ))) {
    stop("Correlation diagonal is not one for: ", result_name)
  }

  write_named_matrix(
    correlation_matrix,
    file.path(
      correlation_directory,
      paste0(result_name, "_correlation_matrix.csv")
    )
  )
  write_named_matrix(
    covariance_matrix,
    file.path(
      sampling_directory,
      paste0(result_name, "_correlation_sampling_vcov.csv")
    )
  )
}

write.csv(
  data.frame(
    block = seq_len(jackknife_blocks),
    jackknife_correlations[[combined_result_name]],
    check.names = FALSE
  ),
  file.path(
    gene_set_output_directory,
    "all_gene_sets_jackknife_correlations.csv"
  ),
  row.names = FALSE
)

write.csv(
  collection_regression_full_sample,
  file.path(regression_directory, "collection_regression_full_sample.csv"),
  row.names = FALSE
)
write.csv(
  collection_regression_jackknife,
  file.path(regression_directory, "collection_regression_jackknife.csv"),
  row.names = FALSE
)
write.csv(
  jackknife_regression,
  file.path(regression_directory, "jackknife_regression_estimates.csv"),
  row.names = FALSE
)
write.csv(
  regression_coefficient_summary,
  file.path(regression_directory, "regression_coefficient_summary.csv"),
  row.names = FALSE
)
write_named_matrix(
  regression_coefficient_vcov,
  file.path(regression_directory, "regression_coefficient_vcov.csv")
)
write.csv(
  collection_plot_data,
  file.path(regression_directory, "collection_point_uncertainty.csv"),
  row.names = FALSE
)
write.csv(
  regression_band,
  file.path(regression_directory, "regression_confidence_band.csv"),
  row.names = FALSE
)
write.csv(
  two_series_full_sample,
  file.path(regression_directory, "two_series_full_sample.csv"),
  row.names = FALSE
)
write.csv(
  two_series_jackknife,
  file.path(regression_directory, "two_series_jackknife.csv"),
  row.names = FALSE
)
write.csv(
  two_series_jackknife_regression,
  file.path(
    regression_directory,
    "two_series_jackknife_regression_estimates.csv"
  ),
  row.names = FALSE
)
write.csv(
  two_series_coefficient_summary,
  file.path(
    regression_directory,
    "two_series_regression_coefficient_summary.csv"
  ),
  row.names = FALSE
)
write.csv(
  two_series_regression_band,
  file.path(
    regression_directory,
    "two_series_regression_confidence_band.csv"
  ),
  row.names = FALSE
)
for (series_id in two_series_ids) {
  write_named_matrix(
    two_series_coefficient_vcov[[series_id]],
    file.path(
      regression_directory,
      paste0("two_series_", series_id, "_coefficient_vcov.csv")
    )
  )
}


# ==============================================================================
# 10. Create the retained paired-series figure
# ==============================================================================

slope_result <- regression_coefficient_summary[term == "slope"]

figure_script <- file.path(
  project_directory,
  "scripts",
  "render_paired_figure.R"
)
figure_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  shQuote(figure_script)
)
if (!identical(figure_status, 0L)) {
  stop("The retained paired-series figure failed to render.", call. = FALSE)
}


# ==============================================================================
# 11. Fit the primary SEMs and create the gene-set workbook
# ==============================================================================

sem_script <- file.path(project_directory, "scripts", "gene_set_enrichment_sem.R")
sem_output_directory <- file.path(
  output_directory,
  "sem",
  "gene_set_enrichment_sem"
)
dir.create(
  sem_output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

sem_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    shQuote(sem_script),
    shQuote(sem_output_directory),
    "both_samples"
  ),
  env = "FULL_JACK_PAPER_RESULTS_MODE=1"
)
if (!identical(sem_status, 0L)) {
  stop(
    "The primary correlated-factor and hierarchical SEM workflow failed.",
    call. = FALSE
  )
}

sem_file <- function(suffix) {
  file.path(
    sem_output_directory,
    paste0("gene_set_enrichment_sem_", suffix, ".csv")
  )
}
sem_required_files <- c(
  sem_file("model_fit"),
  sem_file("model_comparison"),
  sem_file("model_diagnostics"),
  sem_file("hierarchical_loadings"),
  sem_file("implied_trait_correlations"),
  sem_file("implied_trait_correlation_jackknife"),
  sem_file("implied_trait_correlation_jackknife_diagnostics"),
  sem_file("hierarchical_implied_trait_correlations"),
  sem_file("hierarchical_implied_trait_correlation_jackknife"),
  sem_file(
    paste0(
      "hierarchical_implied_trait_correlation_",
      "jackknife_diagnostics"
    )
  )
)
sem_missing_files <- sem_required_files[!file.exists(sem_required_files)]
if (length(sem_missing_files) > 0L) {
  stop(
    "The SEM workflow did not produce: ",
    paste(sem_missing_files, collapse = ", "),
    call. = FALSE
  )
}

paper_model_fit <- read.csv(
  sem_file("model_fit"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
paper_model_comparison <- read.csv(
  sem_file("model_comparison"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
paper_model_diagnostics <- read.csv(
  sem_file("model_diagnostics"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
paper_hierarchical_loadings <- read.csv(
  sem_file("hierarchical_loadings"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
paper_implied_trait_correlations <- read.csv(
  sem_file("implied_trait_correlations"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expected_implied_correlation_columns <- c(
  "analysis",
  "model",
  "trait_1",
  "trait_2",
  "correlation",
  "standard_error",
  "ci_95_lower",
  "ci_95_upper",
  "z",
  "p_value",
  "jackknife_blocks",
  "converged_replicates",
  "admissible_replicates"
)
if (
  !identical(
    names(paper_implied_trait_correlations),
    expected_implied_correlation_columns
  ) ||
    nrow(paper_implied_trait_correlations) != 30L ||
    !setequal(
      paper_implied_trait_correlations$model,
      c(
        "Correlated-factors model",
        "Hierarchical general-factor model"
      )
    ) ||
    any(!is.finite(paper_implied_trait_correlations$correlation)) ||
    any(abs(paper_implied_trait_correlations$correlation) > 1 + 1e-8) ||
    any(!is.finite(
      paper_implied_trait_correlations$standard_error
    )) ||
    any(paper_implied_trait_correlations$standard_error < 0) ||
    any(
      paper_implied_trait_correlations$jackknife_blocks !=
        jackknife_blocks
    ) ||
    anyDuplicated(
      paper_implied_trait_correlations[
        ,
        c("model", "trait_1", "trait_2")
      ]
    )
) {
  stop(
    "The pooled gene-set implied-correlation table is invalid.",
    call. = FALSE
  )
}

diagnostic_columns <- c(
  "model",
  "converged",
  "gradient_check_passed",
  "post_estimation_check",
  "negative_observed_residual_variances",
  "minimum_observed_residual_variance"
)
if (!all(diagnostic_columns %in% names(paper_model_diagnostics))) {
  stop("The SEM diagnostics table lacks required columns.", call. = FALSE)
}
model_order <- paper_model_fit$model
paper_model_fit <- merge(
  paper_model_fit,
  paper_model_diagnostics[, diagnostic_columns, drop = FALSE],
  by = "model",
  all.x = TRUE,
  sort = FALSE
)
paper_model_fit <- paper_model_fit[
  match(model_order, paper_model_fit$model),
  ,
  drop = FALSE
]
rownames(paper_model_fit) <- NULL

trait_display_labels_half_1 <- c(
  agree = "Agreeableness (half 1)",
  consc = "Conscientiousness (half 1)",
  extra = "Extraversion (half 1)",
  neurot = "Neuroticism (half 1)",
  open = "Openness (half 1)",
  iq = "IQ (female)"
)
trait_display_labels_half_2 <- c(
  agree = "Agreeableness (half 2)",
  consc = "Conscientiousness (half 2)",
  extra = "Extraversion (half 2)",
  neurot = "Neuroticism (half 2)",
  open = "Openness (half 2)",
  iq = "IQ (male)"
)

paper_correlation_rows <- lapply(
  result_names,
  function(result_name) {
    correlation_matrix <- point_correlations[[result_name]]
    covariance_matrix <- sampling_vcov[[result_name]]
    rows <- vector(
      "list",
      nrow(trait_halves) * nrow(trait_halves)
    )
    row_index <- 0L

    for (trait_1_index in seq_len(nrow(trait_halves))) {
      for (trait_2_index in seq_len(nrow(trait_halves))) {
        row_index <- row_index + 1L
        half_1_name <- trait_halves$half_1[[trait_1_index]]
        half_2_name <- trait_halves$half_2[[trait_2_index]]
        covariance_name <- pair_column_name(half_1_name, half_2_name)

        rows[[row_index]] <- data.frame(
          gene_set_collection = result_name,
          trait_1 = unname(
            trait_display_labels_half_1[
              trait_halves$trait[[trait_1_index]]
            ]
          ),
          trait_2 = unname(
            trait_display_labels_half_2[
              trait_halves$trait[[trait_2_index]]
            ]
          ),
          correlation = correlation_matrix[
            half_1_name,
            half_2_name
          ],
          n_gene_sets = unname(n_gene_sets_by_result[[result_name]]),
          standard_error = sqrt(
            covariance_matrix[covariance_name, covariance_name]
          ),
          stringsAsFactors = FALSE
        )
      }
    }

    do.call(rbind, rows)
  }
)
paper_split_half_correlations <- do.call(
  rbind,
  paper_correlation_rows
)
rownames(paper_split_half_correlations) <- NULL

expected_correlation_rows <-
  length(result_names) * nrow(trait_halves)^2
if (
  nrow(paper_split_half_correlations) != expected_correlation_rows ||
    any(!is.finite(paper_split_half_correlations$correlation)) ||
    any(!is.finite(paper_split_half_correlations$standard_error)) ||
    any(paper_split_half_correlations$n_gene_sets < 2L)
) {
  stop(
    "The paper split-half correlation table failed validation.",
    call. = FALSE
  )
}

overall_regression_results <- data.frame(
  analysis = "All cross-trait correlations",
  term = regression_coefficient_summary$term,
  estimate = regression_coefficient_summary$full_sample_estimate,
  standard_error = regression_coefficient_summary$jackknife_se,
  ci_95_lower = regression_coefficient_summary$normal_95_low,
  ci_95_upper = regression_coefficient_summary$normal_95_high,
  r_squared = summary(full_regression_fit)$r.squared,
  n_gene_set_collections = length(regression_collections),
  stringsAsFactors = FALSE
)
separate_regression_results <- data.frame(
  analysis = two_series_coefficient_summary$series_label,
  term = two_series_coefficient_summary$term,
  estimate = two_series_coefficient_summary$full_sample_estimate,
  standard_error = two_series_coefficient_summary$jackknife_se,
  ci_95_lower = two_series_coefficient_summary$normal_95_low,
  ci_95_upper = two_series_coefficient_summary$normal_95_high,
  r_squared = two_series_coefficient_summary$full_sample_r_squared,
  n_gene_set_collections = length(regression_collections),
  stringsAsFactors = FALSE
)
paper_regression_results <- rbind(
  overall_regression_results,
  separate_regression_results
)
rownames(paper_regression_results) <- NULL

paper_metadata <- data.frame(
  item = c(
    "Workbook purpose",
    "Canonical analysis entry point",
    "Included correlation block",
    "Correlation rows",
    "Gene-set result sets",
    "Jackknife blocks",
    "Excluded gene-set collection",
    "Excluded gene set",
    "Primary SEM input",
    "SEM estimator",
    "SEM models",
    "Loading standardization",
    "Jackknife scaling note",
    "Implied trait correlations",
    "Implied correlation method",
    "Implied correlation rows",
    "R version",
    "lavaan version",
    "Generated"
  ),
  value = c(
    "Key pooled gene-set paper and supplement results",
    "scripts/analysis.R",
    paste0(
      "All 6 x 6 half-1-by-half-2 correlations, including the six ",
      "same-trait split-half reliabilities"
    ),
    as.character(nrow(paper_split_half_correlations)),
    as.character(length(result_names)),
    as.character(jackknife_blocks),
    excluded_collections,
    excluded_gene_sets,
    "Pooled all_gene_sets correlation matrix and jackknife sampling VCOV",
    "WLS with the full delete-one jackknife sampling covariance",
    paste0(
      "Correlated six-factor model and nested hierarchical general-factor ",
      "model; both include the two orthogonal sample-specific factors"
    ),
    "Fully standardized (std.all) loadings with delta-method standard errors",
    paste0(
      "N=200 converts the direct jackknife VCOV to lavaan NACOV; it is not ",
      "an independent-observation count"
    ),
    "Both correlated-factor and hierarchical general-factor models",
    paste0(
      "Full-sample latent trait correlations with SEs from 200 ",
      "delete-one-block SEM refits using the fixed full-sample WLS weight ",
      "matrix"
    ),
    as.character(nrow(paper_implied_trait_correlations)),
    R.version.string,
    as.character(utils::packageVersion("lavaan")),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  stringsAsFactors = FALSE
)

records_from_data_frame <- function(x) {
  rownames(x) <- NULL
  lapply(
    seq_len(nrow(x)),
    function(index) as.list(x[index, , drop = FALSE])
  )
}
paper_results_payload <- list(
  metadata = records_from_data_frame(paper_metadata),
  split_half_correlations = records_from_data_frame(
    paper_split_half_correlations
  ),
  model_fit = records_from_data_frame(paper_model_fit),
  model_comparison = records_from_data_frame(
    paper_model_comparison
  ),
  hierarchical_loadings = records_from_data_frame(
    paper_hierarchical_loadings
  ),
  implied_trait_correlations = records_from_data_frame(
    paper_implied_trait_correlations
  ),
  regression_results = records_from_data_frame(
    paper_regression_results
  )
)

paper_workbook_path <- Sys.getenv(
  "PAPER_RESULTS_WORKBOOK",
  unset = file.path(output_directory, "gene_set.xlsx")
)
paper_workbook_path <- path.expand(paper_workbook_path)
if (!grepl("^/", paper_workbook_path)) {
  paper_workbook_path <- file.path(
    project_directory,
    paper_workbook_path
  )
}
paper_preview_directory <- Sys.getenv(
  "PAPER_RESULTS_PREVIEW_DIRECTORY",
  unset = file.path(project_directory, "tmp", "gene_set_previews")
)
paper_preview_directory <- path.expand(paper_preview_directory)
if (!grepl("^/", paper_preview_directory)) {
  paper_preview_directory <- file.path(
    project_directory,
    paper_preview_directory
  )
}

artifact_node_modules <- Sys.getenv(
  "ARTIFACT_TOOL_NODE_MODULES",
  unset = ""
)
artifact_node <- Sys.getenv(
  "ARTIFACT_TOOL_NODE",
  unset = "node"
)
if (!nzchar(artifact_node_modules) ||
    !dir.exists(artifact_node_modules)) {
  stop(
    paste0(
      "ARTIFACT_TOOL_NODE_MODULES must point to the bundled Node ",
      "dependency directory so the gene-set workbook can be authored."
    ),
    call. = FALSE
  )
}

paper_temporary_directory <- tempfile("paper_results_")
dir.create(paper_temporary_directory, recursive = TRUE)
on.exit(
  unlink(
    paper_temporary_directory,
    recursive = TRUE,
    force = TRUE
  ),
  add = TRUE
)
paper_results_json <- file.path(
  paper_temporary_directory,
  "paper_results.json"
)
jsonlite::write_json(
  paper_results_payload,
  paper_results_json,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA,
  na = "null"
)

paper_builder_source <- file.path(
  project_directory,
  "scripts",
  "build_paper_results_workbook.mjs"
)
paper_builder_runtime <- file.path(
  paper_temporary_directory,
  basename(paper_builder_source)
)
if (!file.copy(paper_builder_source, paper_builder_runtime)) {
  stop("Failed to prepare the paper-workbook builder.", call. = FALSE)
}
if (!file.symlink(
  artifact_node_modules,
  file.path(paper_temporary_directory, "node_modules")
)) {
  stop("Failed to link the bundled workbook dependencies.", call. = FALSE)
}

unlink(
  paper_preview_directory,
  recursive = TRUE,
  force = TRUE
)
dir.create(
  paper_preview_directory,
  recursive = TRUE,
  showWarnings = FALSE
)
paper_workbook_status <- system2(
  artifact_node,
  c(
    shQuote(paper_builder_runtime),
    shQuote(paper_results_json),
    shQuote(paper_workbook_path),
    shQuote(paper_preview_directory)
  )
)
if (!identical(paper_workbook_status, 0L) ||
    !file.exists(paper_workbook_path)) {
  stop("Gene-set workbook generation failed.", call. = FALSE)
}


# ==============================================================================
# 12. Estimate gene-score correlations, fit SEMs, and build their workbook
# ==============================================================================

gene_score_analysis_script <- file.path(
  project_directory,
  "scripts",
  "gene_score_analysis.R"
)
gene_score_workbook_path <- file.path(
  output_directory,
  "gene_score.xlsx"
)
gene_score_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  shQuote(gene_score_analysis_script),
  env = c(
    "FULL_JACK_GENE_SCORE_MODE=1",
    paste0("ARTIFACT_TOOL_NODE_MODULES=", artifact_node_modules),
    paste0("ARTIFACT_TOOL_NODE=", artifact_node)
  )
)
if (!identical(gene_score_status, 0L) ||
    !file.exists(gene_score_workbook_path)) {
  stop(
    "The gene-score correlation, SEM, or workbook workflow failed.",
    call. = FALSE
  )
}


# ==============================================================================
# 13. Completion summary
# ==============================================================================

message(
  "Wrote ",
  length(result_names),
  " correlation matrices and ",
  length(result_names),
  " sampling covariance matrices, plus collection regression results and ",
  "one retained figure, to ",
  output_directory,
  ". Regression used ",
  length(regression_collections),
  " collections; slope = ",
  format(full_regression_coefficients[["slope"]], digits = 8),
  " (jackknife SE ",
  format(slope_result$jackknife_se, digits = 8),
  "). Gene-set workbook: ",
  paper_workbook_path,
  ". Gene-score workbook: ",
  gene_score_workbook_path,
  "."
)
