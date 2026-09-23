#!/usr/bin/env Rscript

# Internal helper for scripts/analysis.R. Computes gene-score correlations,
# their 200-block delete-one jackknife covariance, the two primary SEMs, and
# the gene-score summary workbook.

required_packages <- c("data.table", "jsonlite")
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

if (!identical(
  Sys.getenv("FULL_JACK_GENE_SCORE_MODE", unset = "0"),
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

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) == 1L) {
  script_path <- normalizePath(
    sub("^--file=", "", script_argument),
    mustWork = TRUE
  )
} else {
  script_path <- normalizePath(
    "scripts/gene_score_analysis.R",
    mustWork = TRUE
  )
}

project_directory <- dirname(dirname(script_path))
gene_score_input_directory <- file.path(
  project_directory,
  "data",
  "gene_enrichments"
)
block_source_directories <- c(
  iq_male_dir = file.path(
    project_directory,
    "data",
    "split_half_output",
    "output",
    "GO_bio_process",
    "jack_iq_male_dir"
  ),
  ReGPC_ext_half_two_no23_dir = file.path(
    project_directory,
    "data",
    "split_half_output",
    "output",
    "akingbuwa_gs",
    "jack_ReGPC_ext_half_two_no23_dir"
  )
)
output_directory <- file.path(project_directory, "output")
correlation_directory <- file.path(output_directory, "correlation_matrices")
sampling_directory <- file.path(output_directory, "sampling_covariances")
gene_score_output_directory <- file.path(output_directory, "gene_scores")
sem_output_directory <- file.path(
  output_directory,
  "sem",
  "gene_score_sem"
)
gene_set_sem_output_directory <- file.path(
  output_directory,
  "sem",
  "gene_set_enrichment_sem"
)
latent_correlation_directory <- file.path(
  output_directory,
  "latent_factor_correlations"
)
workbook_path <- file.path(output_directory, "gene_score.xlsx")
preview_directory <- file.path(
  project_directory,
  "tmp",
  "gene_score_previews"
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
sem_trait_order <- c(
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
trait_labels <- c(
  ReGPC_agr_half_one_no23_dir = "Agreeableness (half 1)",
  ReGPC_agr_half_two_no23_dir = "Agreeableness (half 2)",
  ReGPC_con_half_one_no23_dir = "Conscientiousness (half 1)",
  ReGPC_con_half_two_no23_dir = "Conscientiousness (half 2)",
  ReGPC_ext_half_one_no23_dir = "Extraversion (half 1)",
  ReGPC_ext_half_two_no23_dir = "Extraversion (half 2)",
  ReGPC_neu_half_one_no23_dir = "Neuroticism (half 1)",
  ReGPC_neu_half_two_no23_dir = "Neuroticism (half 2)",
  ReGPC_ope_half_one_no23_dir = "Openness (half 1)",
  ReGPC_ope_half_two_no23_dir = "Openness (half 2)",
  iq_female_dir = "IQ (female)",
  iq_male_dir = "IQ (male)"
)
trait_constructs <- c(
  ReGPC_agr_half_one_no23_dir = "Agreeableness",
  ReGPC_agr_half_two_no23_dir = "Agreeableness",
  ReGPC_con_half_one_no23_dir = "Conscientiousness",
  ReGPC_con_half_two_no23_dir = "Conscientiousness",
  ReGPC_ext_half_one_no23_dir = "Extraversion",
  ReGPC_ext_half_two_no23_dir = "Extraversion",
  ReGPC_neu_half_one_no23_dir = "Neuroticism",
  ReGPC_neu_half_two_no23_dir = "Neuroticism",
  ReGPC_ope_half_one_no23_dir = "Openness",
  ReGPC_ope_half_two_no23_dir = "Openness",
  iq_female_dir = "IQ",
  iq_male_dir = "IQ"
)
jackknife_blocks <- 200L

if (!dir.exists(gene_score_input_directory)) {
  stop(
    "Gene-score input directory not found: ",
    gene_score_input_directory,
    call. = FALSE
  )
}
missing_block_source_directories <- block_source_directories[
  !dir.exists(block_source_directories)
]
if (length(missing_block_source_directories) > 0L) {
  stop(
    "Jackknife block-source directories not found: ",
    paste(missing_block_source_directories, collapse = ", "),
    call. = FALSE
  )
}
dir.create(correlation_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(sampling_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(
  gene_score_output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(sem_output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(
  latent_correlation_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 1. Helpers
# ==============================================================================

read_raw_gene_scores <- function(path, include_positions = FALSE) {
  fields <- if (include_positions) {
    "$1, $2, $3, $4, $9"
  } else {
    "$1, $9"
  }
  command <- paste0(
    "awk '!/^#/ {print ",
    fields,
    "}' ",
    shQuote(path)
  )
  result <- fread(
    cmd = command,
    header = FALSE,
    showProgress = FALSE
  )
  setnames(
    result,
    if (include_positions) {
      c("gene", "chromosome", "start", "stop", "z")
    } else {
      c("gene", "z")
    }
  )
  result[, gene := as.character(gene)]

  if (anyDuplicated(result$gene)) {
    stop("Duplicated gene IDs in: ", path, call. = FALSE)
  }
  if (
    nrow(result) < 2L ||
      anyNA(result$gene) ||
      anyNA(result$z) ||
      !all(is.finite(result$z))
  ) {
    stop("Invalid gene-score values in: ", path, call. = FALSE)
  }

  result
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
    stop(context, " contains non-finite values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(
    matrix_object,
    t(matrix_object),
    tolerance = 1e-10
  ))) {
    stop(context, " is not symmetric.", call. = FALSE)
  }
  invisible(TRUE)
}

records_from_data_frame <- function(x) {
  rownames(x) <- NULL
  lapply(
    seq_len(nrow(x)),
    function(index) as.list(x[index, , drop = FALSE])
  )
}

matrix_records <- function(matrix_object) {
  matrix_object <- matrix_object[
    sem_trait_order,
    sem_trait_order,
    drop = FALSE
  ]
  result <- data.frame(
    Indicator = unname(trait_labels[rownames(matrix_object)]),
    matrix_object,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(result)[-1L] <- unname(trait_labels[colnames(matrix_object)])
  result
}


# ==============================================================================
# 2. Read and align the raw MAGMA Z statistics
# ==============================================================================

gene_score_paths <- file.path(
  gene_score_input_directory,
  paste0(trait_names, ".genes.raw")
)
names(gene_score_paths) <- trait_names
if (!all(file.exists(gene_score_paths))) {
  stop(
    "Missing raw MAGMA gene-score files: ",
    paste(basename(gene_score_paths[!file.exists(gene_score_paths)]), collapse = ", "),
    call. = FALSE
  )
}

all_trait_scores <- lapply(
  gene_score_paths,
  read_raw_gene_scores,
  include_positions = TRUE
)
reference_scores <- all_trait_scores[["iq_male_dir"]]

position_rows <- rbindlist(
  lapply(
    all_trait_scores,
    function(x) x[, .(gene, chromosome, start, stop)]
  )
)
position_conflicts <- position_rows[
  ,
  .(
    chromosome_values = uniqueN(chromosome),
    start_values = uniqueN(start),
    stop_values = uniqueN(stop)
  ),
  by = gene
][
  chromosome_values != 1L |
    start_values != 1L |
    stop_values != 1L
]
if (nrow(position_conflicts) > 0L) {
  stop(
    "At least one gene has conflicting positions across MAGMA files.",
    call. = FALSE
  )
}
gene_positions <- unique(position_rows, by = "gene")

score_table <- gene_positions[, .(gene)]
for (trait in trait_names) {
  trait_scores <- all_trait_scores[[trait]][, .(gene, z)]
  score_table[
    trait_scores,
    on = "gene",
    (trait) := i.z
  ]
}

score_matrix <- as.matrix(score_table[, ..trait_names])
storage.mode(score_matrix) <- "double"
rownames(score_matrix) <- score_table$gene
if (any(colSums(is.finite(score_matrix)) < 2L)) {
  stop("One or more traits have fewer than two gene scores.", call. = FALSE)
}


# ==============================================================================
# 3. Recover the retained 200 gene-exclusion blocks
# ==============================================================================

expected_block_ids <- sprintf("%03d", seq_len(jackknife_blocks))
block_assignments <- rbindlist(
  lapply(
    names(block_source_directories),
    function(source_trait) {
      block_source_files <- list.files(
        block_source_directories[[source_trait]],
        pattern = "_jk[0-9]{3}[.]gsa[.]genes[.]out$",
        full.names = TRUE
      )
      block_ids <- sub(
        "^.*_jk([0-9]{3})[.]gsa[.]genes[.]out$",
        "\\1",
        block_source_files
      )
      if (
        length(block_source_files) != jackknife_blocks ||
          anyDuplicated(block_ids) ||
          !setequal(block_ids, expected_block_ids)
      ) {
        stop(
          "Expected exactly 200 retained gene-output files for block ",
          "recovery from ",
          source_trait,
          ".",
          call. = FALSE
        )
      }
      block_source_files <- setNames(
        block_source_files,
        block_ids
      )[expected_block_ids]
      source_gene_ids <- all_trait_scores[[source_trait]]$gene

      rbindlist(
        lapply(
          seq_len(jackknife_blocks),
          function(block) {
            included_genes <- fread(
              block_source_files[[block]],
              skip = "GENE",
              select = "GENE",
              showProgress = FALSE
            )[[1L]]
            excluded_genes <- setdiff(
              source_gene_ids,
              as.character(included_genes)
            )
            data.table(
              gene = excluded_genes,
              block = block,
              source_trait = source_trait
            )
          }
        )
      )
    }
  )
)
conflicting_block_assignments <- block_assignments[
  ,
  .(block_values = uniqueN(block)),
  by = gene
][block_values != 1L]
if (nrow(conflicting_block_assignments) > 0L) {
  stop(
    "Recovered jackknife sources disagree on at least one gene block.",
    call. = FALSE
  )
}
block_assignments <- unique(
  block_assignments[, .(gene, block)],
  by = "gene"
)
if (
  nrow(block_assignments) == 0L ||
    !all(block_assignments$block %in% seq_len(jackknife_blocks))
) {
  stop("Failed to recover jackknife gene blocks.", call. = FALSE)
}

gene_blocks <- block_assignments$block[
  match(score_table$gene, block_assignments$gene)
]
block_sizes <- tabulate(
  block_assignments$block,
  nbins = jackknife_blocks
)
if (any(block_sizes == 0L)) {
  stop("At least one recovered jackknife block is empty.", call. = FALSE)
}


# ==============================================================================
# 4. Point correlations and delete-one-block sampling covariance
# ==============================================================================

point_correlation <- cor(
  score_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)
validate_symmetric_matrix(point_correlation, "Gene-score correlation matrix")
if (max(abs(diag(point_correlation) - 1)) > 1e-10) {
  stop("Gene-score correlation diagonal is not one.", call. = FALSE)
}

pair_indices <- t(combn(seq_along(trait_names), 2L))
pair_names <- paste(
  trait_names[pair_indices[, 1L]],
  trait_names[pair_indices[, 2L]],
  sep = "__"
)
pair_gene_counts <- vapply(
  seq_len(nrow(pair_indices)),
  function(index) {
    sum(
      complete.cases(
        score_matrix[
          ,
          pair_indices[index, ],
          drop = FALSE
        ]
      )
    )
  },
  integer(1)
)

jackknife_correlations <- matrix(
  NA_real_,
  nrow = jackknife_blocks,
  ncol = length(pair_names),
  dimnames = list(
    sprintf("jk%03d", seq_len(jackknife_blocks)),
    pair_names
  )
)
for (block in seq_len(jackknife_blocks)) {
  keep <- is.na(gene_blocks) | gene_blocks != block
  replicate_correlation <- cor(
    score_matrix[keep, , drop = FALSE],
    use = "pairwise.complete.obs",
    method = "pearson"
  )
  jackknife_correlations[block, ] <- replicate_correlation[pair_indices]

  if (block %% 25L == 0L || block == jackknife_blocks) {
    message(
      "Completed gene-score jackknife block ",
      block,
      " of ",
      jackknife_blocks
    )
  }
}
if (!all(is.finite(jackknife_correlations))) {
  stop("Gene-score jackknife produced non-finite correlations.", call. = FALSE)
}

replicate_means <- colMeans(jackknife_correlations)
centered_estimates <- sweep(
  jackknife_correlations,
  MARGIN = 2L,
  STATS = replicate_means,
  FUN = "-"
)
sampling_vcov <- (
  (jackknife_blocks - 1) / jackknife_blocks
) * crossprod(centered_estimates)
dimnames(sampling_vcov) <- list(pair_names, pair_names)
validate_symmetric_matrix(
  sampling_vcov,
  "Gene-score correlation sampling covariance"
)

sampling_eigenvalues <- eigen(
  sampling_vcov,
  symmetric = TRUE,
  only.values = TRUE
)$values
if (min(sampling_eigenvalues) <= 0) {
  stop(
    "Gene-score sampling covariance is not positive definite.",
    call. = FALSE
  )
}

standard_error_matrix <- matrix(
  0,
  nrow = length(trait_names),
  ncol = length(trait_names),
  dimnames = list(trait_names, trait_names)
)
pair_standard_errors <- sqrt(diag(sampling_vcov))
for (index in seq_len(nrow(pair_indices))) {
  left <- pair_indices[index, 1L]
  right <- pair_indices[index, 2L]
  standard_error_matrix[left, right] <- pair_standard_errors[[index]]
  standard_error_matrix[right, left] <- pair_standard_errors[[index]]
}

correlation_path <- file.path(
  correlation_directory,
  "gene_score_correlation_matrix.csv"
)
sampling_vcov_path <- file.path(
  sampling_directory,
  "gene_score_correlation_sampling_vcov.csv"
)
standard_error_path <- file.path(
  gene_score_output_directory,
  "gene_score_correlation_standard_errors.csv"
)
jackknife_path <- file.path(
  gene_score_output_directory,
  "gene_score_jackknife_correlations.csv"
)
block_assignment_path <- file.path(
  gene_score_output_directory,
  "gene_score_block_assignments.csv"
)

write_named_matrix(point_correlation, correlation_path)
write_named_matrix(sampling_vcov, sampling_vcov_path)
write_named_matrix(standard_error_matrix, standard_error_path)
write.csv(
  data.frame(
    block = seq_len(jackknife_blocks),
    jackknife_correlations,
    check.names = FALSE
  ),
  jackknife_path,
  row.names = FALSE
)
block_output <- merge(
  gene_positions,
  block_assignments,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
setorder(block_output, chromosome, start, stop, gene)
write.csv(block_output, block_assignment_path, row.names = FALSE, na = "")


# ==============================================================================
# 5. Fit the correlated-factor and hierarchical SEMs
# ==============================================================================

sem_script <- file.path(
  project_directory,
  "scripts",
  "gene_set_enrichment_sem.R"
)
sem_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    shQuote(sem_script),
    shQuote(sem_output_directory),
    "both_samples"
  ),
  env = c(
    "FULL_JACK_PAPER_RESULTS_MODE=1",
    "FULL_JACK_SEM_ANALYSIS_KIND=gene_scores",
    paste0("FULL_JACK_SEM_CORRELATION_PATH=", correlation_path),
    paste0("FULL_JACK_SEM_SAMPLING_VCOV_PATH=", sampling_vcov_path),
    paste0(
      "FULL_JACK_SEM_JACKKNIFE_CORRELATION_PATH=",
      jackknife_path
    )
  )
)
if (!identical(sem_status, 0L)) {
  stop("The gene-score SEM workflow failed.", call. = FALSE)
}

sem_file <- function(suffix) {
  file.path(
    sem_output_directory,
    paste0("gene_score_sem_", suffix, ".csv")
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
    "The gene-score SEM did not produce: ",
    paste(sem_missing_files, collapse = ", "),
    call. = FALSE
  )
}

model_fit <- read.csv(
  sem_file("model_fit"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
model_comparison <- read.csv(
  sem_file("model_comparison"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
model_diagnostics <- read.csv(
  sem_file("model_diagnostics"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
hierarchical_loadings <- read.csv(
  sem_file("hierarchical_loadings"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gene_score_implied_trait_correlations <- read.csv(
  sem_file("hierarchical_implied_trait_correlations"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gene_score_model_implied_trait_correlations <- read.csv(
  sem_file("implied_trait_correlations"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_set_sem_file <- function(suffix) {
  file.path(
    gene_set_sem_output_directory,
    paste0("gene_set_enrichment_sem_", suffix, ".csv")
  )
}
gene_set_required_files <- c(
  gene_set_sem_file("implied_trait_correlations"),
  gene_set_sem_file("implied_trait_correlation_jackknife"),
  gene_set_sem_file(
    "implied_trait_correlation_jackknife_diagnostics"
  ),
  gene_set_sem_file("hierarchical_implied_trait_correlations"),
  gene_set_sem_file(
    "hierarchical_implied_trait_correlation_jackknife"
  ),
  gene_set_sem_file(
    paste0(
      "hierarchical_implied_trait_correlation_",
      "jackknife_diagnostics"
    )
  )
)
gene_set_missing_files <- gene_set_required_files[
  !file.exists(gene_set_required_files)
]
if (length(gene_set_missing_files) > 0L) {
  stop(
    paste0(
      "The pooled gene-set hierarchical jackknife outputs are missing. ",
      "Run the canonical scripts/analysis.R workflow so the gene-set ",
      "SEM precedes the gene-score analysis. Missing: ",
      paste(basename(gene_set_missing_files), collapse = ", ")
    ),
    call. = FALSE
  )
}
gene_set_implied_trait_correlations <- read.csv(
  gene_set_sem_file("hierarchical_implied_trait_correlations"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
hierarchical_implied_trait_correlations <- rbind(
  gene_set_implied_trait_correlations,
  gene_score_implied_trait_correlations
)
expected_implied_columns <- c(
  "analysis",
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
    names(hierarchical_implied_trait_correlations),
    expected_implied_columns
  ) ||
    nrow(hierarchical_implied_trait_correlations) != 30L ||
    any(!is.finite(
      hierarchical_implied_trait_correlations$correlation
    )) ||
    any(abs(
      hierarchical_implied_trait_correlations$correlation
    ) > 1 + 1e-8) ||
    any(!is.finite(
      hierarchical_implied_trait_correlations$standard_error
    )) ||
    any(
      hierarchical_implied_trait_correlations$standard_error < 0
    ) ||
    any(
      hierarchical_implied_trait_correlations$jackknife_blocks !=
        jackknife_blocks
    ) ||
    any(
      hierarchical_implied_trait_correlations$converged_replicates !=
        jackknife_blocks
    ) ||
    any(
      hierarchical_implied_trait_correlations$admissible_replicates < 0 |
        hierarchical_implied_trait_correlations$admissible_replicates >
          jackknife_blocks
    ) ||
    anyDuplicated(
      hierarchical_implied_trait_correlations[
        ,
        c("analysis", "trait_1", "trait_2")
      ]
    )
) {
  stop(
    "The combined hierarchical implied-correlation table is invalid.",
    call. = FALSE
  )
}
latent_correlation_path <- file.path(
  latent_correlation_directory,
  "hierarchical_implied_trait_correlations.csv"
)
write.csv(
  hierarchical_implied_trait_correlations,
  latent_correlation_path,
  row.names = FALSE,
  na = ""
)

gene_set_implied_correlation_jackknife <- read.csv(
  gene_set_sem_file(
    "hierarchical_implied_trait_correlation_jackknife"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gene_score_implied_correlation_jackknife <- read.csv(
  sem_file("hierarchical_implied_trait_correlation_jackknife"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
combined_implied_correlation_jackknife <- rbind(
  gene_set_implied_correlation_jackknife,
  gene_score_implied_correlation_jackknife
)
if (
  nrow(combined_implied_correlation_jackknife) !=
    2L * jackknife_blocks ||
    anyDuplicated(
      combined_implied_correlation_jackknife[
        ,
        c("analysis", "block")
      ]
    ) ||
    any(!is.finite(
      as.matrix(
        combined_implied_correlation_jackknife[
          ,
          setdiff(
            names(combined_implied_correlation_jackknife),
            c("analysis", "block")
          ),
          drop = FALSE
        ]
      )
    ))
) {
  stop(
    "The combined hierarchical implied-correlation replicates are invalid.",
    call. = FALSE
  )
}
write.csv(
  combined_implied_correlation_jackknife,
  file.path(
    latent_correlation_directory,
    paste0(
      "hierarchical_implied_trait_correlation_",
      "jackknife_estimates.csv"
    )
  ),
  row.names = FALSE,
  na = ""
)

gene_set_implied_correlation_diagnostics <- read.csv(
  gene_set_sem_file(
    paste0(
      "hierarchical_implied_trait_correlation_",
      "jackknife_diagnostics"
    )
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gene_score_implied_correlation_diagnostics <- read.csv(
  sem_file(
    paste0(
      "hierarchical_implied_trait_correlation_",
      "jackknife_diagnostics"
    )
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
combined_implied_correlation_diagnostics <- rbind(
  gene_set_implied_correlation_diagnostics,
  gene_score_implied_correlation_diagnostics
)
if (
  nrow(combined_implied_correlation_diagnostics) !=
    2L * jackknife_blocks ||
    anyDuplicated(
      combined_implied_correlation_diagnostics[
        ,
        c("analysis", "block")
      ]
    )
) {
  stop(
    "The combined hierarchical jackknife diagnostics are invalid.",
    call. = FALSE
  )
}
write.csv(
  combined_implied_correlation_diagnostics,
  file.path(
    latent_correlation_directory,
    paste0(
      "hierarchical_implied_trait_correlation_",
      "jackknife_diagnostics.csv"
    )
  ),
  row.names = FALSE,
  na = ""
)

gene_set_model_implied_trait_correlations <- read.csv(
  gene_set_sem_file("implied_trait_correlations"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
model_implied_trait_correlations <- rbind(
  gene_set_model_implied_trait_correlations,
  gene_score_model_implied_trait_correlations
)
expected_model_implied_columns <- c(
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
    names(model_implied_trait_correlations),
    expected_model_implied_columns
  ) ||
    nrow(model_implied_trait_correlations) != 60L ||
    !setequal(
      model_implied_trait_correlations$model,
      c(
        "Correlated-factors model",
        "Hierarchical general-factor model"
      )
    ) ||
    any(!is.finite(model_implied_trait_correlations$correlation)) ||
    any(abs(model_implied_trait_correlations$correlation) > 1 + 1e-8) ||
    any(!is.finite(
      model_implied_trait_correlations$standard_error
    )) ||
    any(model_implied_trait_correlations$standard_error < 0) ||
    any(
      model_implied_trait_correlations$jackknife_blocks !=
        jackknife_blocks
    ) ||
    anyDuplicated(
      model_implied_trait_correlations[
        ,
        c("analysis", "model", "trait_1", "trait_2")
      ]
    )
) {
  stop(
    "The combined model-implied trait-correlation table is invalid.",
    call. = FALSE
  )
}
write.csv(
  model_implied_trait_correlations,
  file.path(
    latent_correlation_directory,
    "implied_trait_correlations.csv"
  ),
  row.names = FALSE,
  na = ""
)

gene_set_model_implied_correlation_jackknife <- read.csv(
  gene_set_sem_file("implied_trait_correlation_jackknife"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gene_score_model_implied_correlation_jackknife <- read.csv(
  sem_file("implied_trait_correlation_jackknife"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
combined_model_implied_correlation_jackknife <- rbind(
  gene_set_model_implied_correlation_jackknife,
  gene_score_model_implied_correlation_jackknife
)
model_jackknife_numeric_columns <- setdiff(
  names(combined_model_implied_correlation_jackknife),
  c("analysis", "model")
)
if (
  nrow(combined_model_implied_correlation_jackknife) !=
    4L * jackknife_blocks ||
    anyDuplicated(
      combined_model_implied_correlation_jackknife[
        ,
        c("analysis", "model", "block")
      ]
    ) ||
    any(!is.finite(
      as.matrix(
        combined_model_implied_correlation_jackknife[
          ,
          model_jackknife_numeric_columns,
          drop = FALSE
        ]
      )
    ))
) {
  stop(
    "The combined model-implied correlation replicates are invalid.",
    call. = FALSE
  )
}
write.csv(
  combined_model_implied_correlation_jackknife,
  file.path(
    latent_correlation_directory,
    "implied_trait_correlation_jackknife_estimates.csv"
  ),
  row.names = FALSE,
  na = ""
)

gene_set_model_implied_correlation_diagnostics <- read.csv(
  gene_set_sem_file(
    "implied_trait_correlation_jackknife_diagnostics"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gene_score_model_implied_correlation_diagnostics <- read.csv(
  sem_file(
    "implied_trait_correlation_jackknife_diagnostics"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
combined_model_implied_correlation_diagnostics <- rbind(
  gene_set_model_implied_correlation_diagnostics,
  gene_score_model_implied_correlation_diagnostics
)
if (
  nrow(combined_model_implied_correlation_diagnostics) !=
    4L * jackknife_blocks ||
    anyDuplicated(
      combined_model_implied_correlation_diagnostics[
        ,
        c("analysis", "model", "block")
      ]
    )
) {
  stop(
    "The combined model-implied jackknife diagnostics are invalid.",
    call. = FALSE
  )
}
write.csv(
  combined_model_implied_correlation_diagnostics,
  file.path(
    latent_correlation_directory,
    "implied_trait_correlation_jackknife_diagnostics.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ==============================================================================
# 6. Build output/gene_score.xlsx
# ==============================================================================

pair_details <- data.frame(
  trait_1 = unname(trait_labels[trait_names[pair_indices[, 1L]]]),
  trait_2 = unname(trait_labels[trait_names[pair_indices[, 2L]]]),
  same_construct = unname(
    trait_constructs[trait_names[pair_indices[, 1L]]]
  ) == unname(
    trait_constructs[trait_names[pair_indices[, 2L]]]
  ),
  n_genes = pair_gene_counts,
  correlation = point_correlation[pair_indices],
  standard_error = pair_standard_errors,
  ci_95_lower = point_correlation[pair_indices] -
    qnorm(0.975) * pair_standard_errors,
  ci_95_upper = point_correlation[pair_indices] +
    qnorm(0.975) * pair_standard_errors,
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  item = c(
    "Workbook purpose",
    "Canonical analysis entry point",
    "Raw MAGMA inputs",
    "Gene-score field",
    "Correlation estimator",
    "Pairwise gene counts",
    "Jackknife strategy",
    "Genes assigned to blocks",
    "Genes retained in every replicate",
    "Recovered block sizes",
    "SEM estimator",
    "SEM models",
    "Loading standardization",
    "Jackknife scaling note",
    "Latent correlation analysis",
    "Latent correlation method",
    "Latent correlation output",
    "Latent jackknife diagnostics",
    "R version",
    "lavaan version",
    "Generated"
  ),
  value = c(
    "Gene-score correlations and structural equation models",
    "scripts/analysis.R",
    "data/gene_enrichments/*.genes.raw",
    "MAGMA ZSTAT (column 9 of each raw gene-results file)",
    "Pearson correlation after alignment by gene ID; pairwise complete genes",
    paste0(
      min(pair_gene_counts),
      " to ",
      max(pair_gene_counts),
      " genes across the 66 unique pairs"
    ),
    paste0(
      "The same 200 delete-one gene-exclusion blocks used by the gene-set ",
      "analysis; membership recovered from retained per-block MAGMA outputs"
    ),
    as.character(nrow(block_assignments)),
    as.character(nrow(score_table) - nrow(block_assignments)),
    paste0(min(block_sizes), " to ", max(block_sizes), " genes"),
    "WLS with the full delete-one jackknife sampling covariance",
    paste0(
      "Correlated six-factor model and nested hierarchical general-factor ",
      "model; both include two orthogonal sample-specific factors"
    ),
    paste0(
      "Fully standardized (std.all) loadings with delta-method standard ",
      "errors in the loading table only"
    ),
    paste0(
      "N=200 converts the direct jackknife VCOV to lavaan NACOV; it is a ",
      "scaling convention, not an independent-observation count"
    ),
    "Gene scores; both correlated-factor and hierarchical models",
    paste0(
      "Full-sample latent trait correlations; SEs from 200 delete-one-block ",
      "SEM refits using the fixed full-sample WLS weight matrix"
    ),
    "output/latent_factor_correlations/implied_trait_correlations.csv",
    paste0(
      "output/latent_factor_correlations/",
      "implied_trait_correlation_jackknife_diagnostics.csv"
    ),
    R.version.string,
    as.character(utils::packageVersion("lavaan")),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  stringsAsFactors = FALSE
)

workbook_payload <- list(
  metadata = records_from_data_frame(metadata),
  correlation_matrix = records_from_data_frame(
    matrix_records(point_correlation)
  ),
  standard_error_matrix = records_from_data_frame(
    matrix_records(standard_error_matrix)
  ),
  pair_details = records_from_data_frame(pair_details),
  model_fit = records_from_data_frame(model_fit),
  model_comparison = records_from_data_frame(model_comparison),
  model_diagnostics = records_from_data_frame(model_diagnostics),
  latent_factor_correlations = records_from_data_frame(
    gene_score_model_implied_trait_correlations
  ),
  hierarchical_loadings = records_from_data_frame(
    hierarchical_loadings
  )
)

artifact_node_modules <- Sys.getenv(
  "ARTIFACT_TOOL_NODE_MODULES",
  unset = ""
)
artifact_node <- Sys.getenv("ARTIFACT_TOOL_NODE", unset = "node")
if (!nzchar(artifact_node_modules) ||
    !dir.exists(artifact_node_modules)) {
  stop(
    paste0(
      "ARTIFACT_TOOL_NODE_MODULES must point to the bundled Node ",
      "dependency directory so the gene-score workbook can be authored."
    ),
    call. = FALSE
  )
}

temporary_directory <- tempfile("gene_score_workbook_")
dir.create(temporary_directory, recursive = TRUE)
on.exit(
  unlink(temporary_directory, recursive = TRUE, force = TRUE),
  add = TRUE
)
payload_path <- file.path(temporary_directory, "gene_score_results.json")
jsonlite::write_json(
  workbook_payload,
  payload_path,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA,
  na = "null"
)

builder_source <- file.path(
  project_directory,
  "scripts",
  "build_gene_score_workbook.mjs"
)
builder_runtime <- file.path(
  temporary_directory,
  basename(builder_source)
)
if (!file.copy(builder_source, builder_runtime)) {
  stop("Failed to prepare the gene-score workbook builder.", call. = FALSE)
}
if (!file.symlink(
  artifact_node_modules,
  file.path(temporary_directory, "node_modules")
)) {
  stop(
    "Failed to link the bundled workbook dependencies.",
    call. = FALSE
  )
}

unlink(preview_directory, recursive = TRUE, force = TRUE)
dir.create(preview_directory, recursive = TRUE, showWarnings = FALSE)
workbook_status <- system2(
  artifact_node,
  c(
    shQuote(builder_runtime),
    shQuote(payload_path),
    shQuote(workbook_path),
    shQuote(preview_directory)
  )
)
if (!identical(workbook_status, 0L) || !file.exists(workbook_path)) {
  stop("Gene-score workbook generation failed.", call. = FALSE)
}

message(
  "Gene-score analysis wrote the correlation matrix, sampling covariance, ",
  "200 jackknife replicates, two SEMs, two-model implied trait ",
  "correlations, and workbook to ",
  output_directory,
  ". Workbook: ",
  workbook_path
)
