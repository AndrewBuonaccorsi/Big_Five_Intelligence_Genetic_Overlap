#!/usr/bin/env Rscript

# ==============================================================================
# Setup and saved correlation inputs
# ==============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, lavaan, writexl)

cores <- 4L
if (.Platform$OS.type == "windows") cores <- 1L
traits <- c("agree", "consc", "extra", "neurot", "open", "iq")
indicators <- c(
  "agr1", "agr2", "con1", "con2", "ext1", "ext2",
  "neu1", "neu2", "ope1", "ope2", "iq_female", "iq_male"
)
indicator_trait <- rep(1:6, each = 2)
pairs <- t(combn(1:12, 2))
trait_pairs <- t(combn(1:6, 2))
residual_tolerance <- 1e-7
boundary_tolerance <- 1e-6
stationarity_tolerance <- 1e-5

# ==============================================================================
# Model specifications: free trait loadings and two orthogonal sample factors
# ==============================================================================

measurement <- "
  agree =~ l_agr1*agr1 + l_agr2*agr2
  consc =~ l_con1*con1 + l_con2*con2
  extra =~ l_ext1*ext1 + l_ext2*ext2
  neurot =~ l_neu1*neu1 + l_neu2*neu2
  open =~ l_ope1*ope1 + l_ope2*ope2
  iq =~ l_iq_female*iq_female + l_iq_male*iq_male
  sample_1 =~ m_agr1*agr1 + m_con1*con1 + m_ext1*ext1 + m_neu1*neu1 + m_ope1*ope1
  sample_2 =~ m_agr2*agr2 + m_con2*con2 + m_ext2*ext2 + m_neu2*neu2 + m_ope2*ope2
"
structures <- list(
  correlated = "
    agree ~~ p_agree_consc*consc + p_agree_extra*extra + p_agree_neurot*neurot + p_agree_open*open + p_agree_iq*iq
    consc ~~ p_consc_extra*extra + p_consc_neurot*neurot + p_consc_open*open + p_consc_iq*iq
    extra ~~ p_extra_neurot*neurot + p_extra_open*open + p_extra_iq*iq
    neurot ~~ p_neurot_open*open + p_neurot_iq*iq
    open ~~ p_open_iq*iq
  ",
  hierarchical = "general =~ g_agree*agree + g_consc*consc + g_extra*extra + g_neurot*neurot + g_open*open + g_iq*iq"
)
models <- imap(structures, \(structural, structure) {
  # correlation=TRUE derives indicator residual variances. std.lv fixes each
  # first-order latent residual to one, so its hierarchical total is 1 + g^2.
  constraints <- map_chr(1:12, \(i) {
    explained <- paste0("l_", indicators[i], "^2")
    if (structure == "hierarchical") {
      explained <- paste0(explained, "*(1 + g_", traits[indicator_trait[i]], "^2)")
    }
    if (i <= 10) explained <- paste0(explained, " + m_", indicators[i], "^2")
    # Lavaan implements < as a weak inequality: zero residuals are permitted.
    paste0(explained, " < 1")
  })
  paste(measurement, structural, paste(constraints, collapse = "\n"), sep = "\n")
})

# ==============================================================================
# Starting values and consistent factor signs
# ==============================================================================

make_starts <- function(r, structure, template) {
  reliability <- r[cbind(seq(1, 11, 2), seq(2, 12, 2))]
  a <- trait_pairs[, 1]
  b <- trait_pairs[, 2]
  cross_12 <- r[cbind(2 * a - 1, 2 * b)]
  cross_21 <- r[cbind(2 * a, 2 * b - 1)]
  rho <- suppressWarnings(sign(cross_12 + cross_21) * sqrt(
    cross_12 * cross_21 / (reliability[a] * reliability[b])
  ))
  rho[!is.finite(rho)] <- 0
  phi <- diag(6)
  phi[trait_pairs] <- phi[trait_pairs[, 2:1]] <- pmin(.8, pmax(-.8, rho))
  eig <- eigen(phi, symmetric = TRUE)
  phi <- cov2cor(eig$vectors %*% diag(pmax(eig$values, 1e-5)) %*% t(eig$vectors))

  # Bounds and positive-definite repair apply only to starting values. They do
  # not change the analytic estimates or constrain the fitted latent correlations.
  starts <- map(1:8, \(i) {
    ratio <- switch(as.character(i),
      `1` = rep(1, 6), `2` = c(rep(1.35, 5), .9),
      `3` = c(1.2, 1.3, 1.25, 1.32, 1.33, .9),
      `4` = c(rep(1.6, 5), .85), exp(rnorm(6, log(1.25), .2))
    )
    first <- sqrt(pmax(.005, abs(reliability)) / ratio)
    second <- ratio * first
    scale_down <- pmax(1, pmax(first, second) / .85)
    loading <- as.vector(rbind(first / scale_down, second / scale_down))
    method <- switch(as.character(i),
      `1` = rep(0, 10), `2` = rep(c(.15, .1), 5),
      `3` = rep(c(.25, -.1), 5), rnorm(10, sd = .18)
    )
    method <- .999 * sqrt(1 - loading[1:10]^2) * tanh(method)
    general <- if (i == 1) rep(.55, 6) else rep(.7, 6)
    if (i >= 3) general <- .999 * tanh(atanh(general / .999) + rnorm(6, sd = .12))
    latent_sd <- if (structure == "hierarchical") sqrt(1 - general^2) else rep(1, 6)
    values <- c(setNames(loading * latent_sd[indicator_trait], paste0("l_", indicators)),
                setNames(method, paste0("m_", indicators[1:10])))
    if (structure == "hierarchical") {
      values <- c(values, setNames(general / latent_sd, paste0("g_", traits)))
    } else {
      correlations <- if (i == 1) rep(0, 15) else phi[trait_pairs]
      values <- c(values, setNames(correlations,
                                   paste0("p_", traits[a], "_", traits[b])))
    }
    labels <- template$label[template$free > 0]
    stopifnot(all(labels %in% names(values)))
    unname(values[labels])
  })
  c(list("default"), starts)
}

align_signs <- function(estimates, reference = NULL) {
  keys <- list(paste0("sample_loading__", indicators[seq(1, 9, 2)]),
               paste0("sample_loading__", indicators[seq(2, 10, 2)]))
  if ("general_loading__iq" %in% names(estimates)) {
    keys <- c(keys, list(paste0("general_loading__", traits)))
  }
  for (group in keys) {
    direction <- if (is.null(reference)) sum(estimates[group]) else sum(estimates[group] * reference[group])
    if (direction < 0) estimates[group] <- -estimates[group]
  }
  estimates
}

# ==============================================================================
# Full WLS estimation and constrained-fit diagnostics
# ==============================================================================

fit_once <- function(r, input, structure, start = "default", do_fit = TRUE) {
  warnings <- character()
  fit <- withCallingHandlers(cfa(
    models[[structure]], sample.cov = r, sample.nobs = 200,
    sample.cov.rescale = FALSE, correlation = TRUE, std.lv = TRUE,
    orthogonal = TRUE, estimator = "WLS", WLS.V = input$precision / 199,
    NACOV = 200 * input$vcov, se = "none", test = "standard",
    start = start, do.fit = do_fit,
    control = list(iter.max = 20000L, eval.max = 50000L, control.outer = list(tol = 1e-10)),
    check.gradient = TRUE, check.post = TRUE
  ), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  if (!do_fit) return(parTable(fit))

  parameters <- coef(fit)
  implied <- fitted(fit)$cov[indicators, indicators]
  theta <- diag(lavInspect(fit, "theta"))[indicators]
  phi <- lavInspect(fit, "cor.lv")[traits, traits]
  delta <- lavInspect(fit, "delta")
  if (is.list(delta)) delta <- delta[[1]]
  # This check protects the off-diagonal order of both the observations and V.
  stopifnot(max(abs(lavInspect(fit, "wls.obs") - r[pairs])) < 1e-10)
  residual <- r[pairs] - implied[pairs]
  q <- drop(crossprod(residual, input$precision %*% residual))
  eigenvalues <- eigen(crossprod(delta, input$precision %*% delta), symmetric = TRUE)$values
  rank <- sum(eigenvalues > max(eigenvalues) * 1e-9)

  # Constrained lavaan can return a zero-filled gradient. Calculate the WLS
  # gradient independently and check the active-constraint KKT conditions.
  gradient <- drop(-crossprod(delta, input$precision %*% residual) / 199)
  active <- which(theta <= boundary_tolerance)
  constraint_gradient <- matrix(0, length(active), length(parameters),
                                dimnames = list(indicators[active], names(parameters)))
  for (row in seq_along(active)) {
    i <- active[row]
    loading <- parameters[paste0("l_", indicators[i])]
    variance <- 1
    if (structure == "hierarchical") {
      general_name <- paste0("g_", traits[indicator_trait[i]])
      general <- parameters[general_name]
      variance <- 1 + general^2
      constraint_gradient[row, general_name] <- -2 * loading^2 * general
    }
    constraint_gradient[row, paste0("l_", indicators[i])] <- -2 * loading * variance
    if (i <= 10) {
      method_name <- paste0("m_", indicators[i])
      constraint_gradient[row, method_name] <- -2 * parameters[method_name]
    }
  }
  multipliers <- if (length(active)) qr.solve(t(constraint_gradient), gradient) else numeric()
  projected <- if (length(active)) gradient - drop(t(constraint_gradient) %*% multipliers) else gradient
  complementarity <- if (length(active)) max(abs(multipliers * theta[active])) else 0
  feasible <- min(theta) >= -residual_tolerance
  kkt <- feasible && max(abs(projected)) <= stationarity_tolerance &&
    all(multipliers >= -stationarity_tolerance) && complementarity <= stationarity_tolerance
  min_phi <- min(eigen(phi, symmetric = TRUE)$values)
  covariance_ok <- feasible && min_phi > 1e-6 && min(eigen(implied, symmetric = TRUE)$values) > 1e-8
  converged <- isTRUE(lavInspect(fit, "converged"))

  standardized <- standardizedSolution(fit, se = FALSE, zstat = FALSE, pvalue = FALSE, ci = FALSE)
  loading <- map_dbl(1:12, \(i) {
    standardized$est.std[standardized$lhs == traits[indicator_trait[i]] &
                          standardized$op == "=~" & standardized$rhs == indicators[i]]
  })
  trait_sign <- sign(loading[seq(1, 11, 2)])
  trait_sign[trait_sign == 0] <- 1
  loading <- loading * trait_sign[indicator_trait]
  phi <- phi * outer(trait_sign, trait_sign)
  method <- map_dbl(1:10, \(i) {
    standardized$est.std[standardized$lhs == paste0("sample_", if (i %% 2) 1 else 2) &
                          standardized$op == "=~" & standardized$rhs == indicators[i]]
  })
  general <- if (structure == "hierarchical") {
    map_dbl(traits, \(trait) {
      standardized$est.std[standardized$lhs == "general" & standardized$op == "=~" &
                            standardized$rhs == trait]
    }) * trait_sign
  } else numeric()
  rho <- phi[trait_pairs]
  iq_pair <- trait_pairs[, 2] == 6
  estimates <- c(
    setNames(loading, paste0("loading__", indicators)),
    setNames(method, paste0("sample_loading__", indicators[1:10])),
    setNames(general, if (length(general)) paste0("general_loading__", traits) else character()),
    setNames(rho, paste0("correlation__", traits[trait_pairs[, 1]], "__", traits[trait_pairs[, 2]])),
    mean_big_five = mean(rho[!iq_pair]), mean_iq_big_five = mean(rho[iq_pair]),
    difference_iq_minus_big_five = mean(rho[iq_pair]) - mean(rho[!iq_pair])
  )
  valid <- converged && covariance_ok && kkt && rank == length(parameters) && all(is.finite(estimates))
  diagnostics <- tibble(
    q = q, df = 66 - length(parameters), valid = valid, converged = converged,
    feasible = feasible, kkt = kkt, boundary = length(active) > 0,
    boundary_residuals = paste(indicators[active], collapse = ","),
    covariance_ok = covariance_ok, information_rank = rank, parameters = length(parameters),
    max_gradient = max(abs(projected)), min_residual = min(theta), min_trait_eigenvalue = min_phi,
    min_multiplier = if (length(active)) min(multipliers) else NA_real_,
    complementarity = complementarity,
    lavaan_post_check = suppressWarnings(isTRUE(lavInspect(fit, "post.check"))),
    message = paste(unique(warnings), collapse = " | ")
  )
  list(fit = fit, parameters = parameters, estimates = align_signs(estimates),
       implied = implied, diagnostics = diagnostics)
}

fit_model <- function(r, input, structure, starts, parallel_starts = FALSE) {
  fit_start <- function(start) tryCatch(fit_once(r, input, structure, start), error = \(e) e)
  candidates <- if (parallel_starts) {
    parallel::mclapply(starts, fit_start, mc.cores = cores)
  } else map(starts, fit_start)
  attempts <- imap_dfr(candidates, \(x, i) {
    if (inherits(x, "error")) {
      tibble(start = i, q = Inf, converged = FALSE, feasible = FALSE, kkt = FALSE,
             valid = FALSE, message = conditionMessage(x))
    } else mutate(x$diagnostics, start = i, .before = 1)
  })
  if (!any(is.finite(attempts$q))) stop(paste(attempts$message, collapse = " | "))
  eligible <- attempts$converged & attempts$feasible & attempts$kkt
  objective <- if (any(eligible)) replace(attempts$q, !eligible, Inf) else attempts$q
  # Select on convergence and optimality, never on a preferable covariance or
  # information-rank diagnosis at a worse-fitting local solution.
  best <- candidates[[which.min(objective)]]
  best$attempts <- attempts
  best
}

# ==============================================================================
# Full-data fits, with resumable fixed-weight jackknife refits
# ==============================================================================

# Tests read definitions above this line without running the paper analysis.
paths <- c(gene_scores = "output/gene_scores/inputs.rds", gene_sets = "output/gene_sets/inputs.rds")
inputs <- map(paths, readRDS)
stopifnot(identical(inputs$gene_scores$block_hashes, inputs$gene_sets$block_hashes))
inputs <- map(inputs, \(x) {
  stopifnot(identical(x$indicators, indicators), identical(x$pairs, pairs),
            nrow(x$jackknife) == 200, all(is.finite(x$jackknife)))
  centered <- sweep(x$jackknife, 2, colMeans(x$jackknife), "-")
  stopifnot(max(abs(x$vcov - (199 / 200) * crossprod(centered))) < 1e-10,
            min(eigen(x$correlation, symmetric = TRUE)$values) > 0,
            min(eigen(x$vcov, symmetric = TRUE)$values) > 0)
  x$precision <- solve(x$vcov)
  x
})
signature <- list(inputs = unname(tools::md5sum(paths)),
                  code = unname(tools::md5sum("scripts/04_sem.R")),
                  lavaan = as.character(packageVersion("lavaan")))
dir.create("output/sem", recursive = TRUE, showWarnings = FALSE)
full_path <- "output/sem/full_fits.rds"
if (file.exists(full_path)) {
  saved <- readRDS(full_path)
  if (!identical(signature, saved$signature)) {
    stop("SEM inputs or fitting code changed. Move the existing output/sem directory before rerunning.")
  }
  fits <- saved$fits
} else {
  fits <- list()
  for (analysis in names(inputs)) {
    input <- inputs[[analysis]]
    for (structure in names(models)) {
      key <- paste(analysis, structure, sep = "__")
      set.seed(3700 + match(analysis, names(inputs)) * 100 + match(structure, names(models)) * 10 + 2)
      template <- fit_once(input$correlation, input, structure, do_fit = FALSE)
      starts <- make_starts(input$correlation, structure, template)
      fits[[key]] <- fit_model(input$correlation, input, structure, starts, parallel_starts = TRUE)
      # A sample factor dominated by one indicator can form a separate optimum.
      # Check each possible dominant indicator using this fit as the starting
      # point; no archived fit or second estimator is needed to find that basin.
      concentrated <- map(1:10, \(i) {
        p <- fits[[key]]$parameters
        half <- if (i %% 2) seq(1, 9, 2) else seq(2, 10, 2)
        p[paste0("m_", indicators[half])] <- .02
        p[paste0("m_", indicators[i])] <- .85
        variance <- if (structure == "hierarchical") {
          1 + p[paste0("g_", traits[indicator_trait[i]])]^2
        } else 1
        name <- paste0("l_", indicators[i])
        p[name] <- sign(p[name]) * min(abs(p[name]), .9 * sqrt((1 - .85^2) / variance))
        unname(p)
      })
      fits[[key]] <- fit_model(input$correlation, input, structure,
                               c(list(unname(fits[[key]]$parameters)), concentrated),
                               parallel_starts = TRUE)
      writeLines(models[[structure]], paste0("output/sem/", key, "__model.lav"))
      write_csv(fits[[key]]$attempts, paste0("output/sem/", key, "__starts.csv"))
      message(key, ": Q = ", signif(fits[[key]]$diagnostics$q, 7),
              "; valid = ", fits[[key]]$diagnostics$valid)
    }
  }
  saveRDS(list(signature = signature, fits = fits), full_path)
}

fit_summary <- list()
estimates <- list()
for (analysis in names(inputs)) {
  input <- inputs[[analysis]]
  for (structure in names(models)) {
    key <- paste(analysis, structure, sep = "__")
    full <- fits[[key]]
    point <- full$estimates
    null_q <- drop(crossprod(input$correlation[pairs], input$precision %*% input$correlation[pairs]))
    residual <- input$correlation[pairs] - full$implied[pairs]
    fit_summary[[key]] <- full$diagnostics |>
      mutate(analysis = analysis, structure = structure, .before = 1,
             p_value = if_else(valid & !boundary, pchisq(q, df, lower.tail = FALSE), NA_real_),
             cfi = 1 - pmax(q - df, 0) / pmax(null_q - 66, q - df, 0),
             tli = (null_q / 66 - q / df) / (null_q / 66 - 1),
             srmr = sqrt(sum(residual^2) / 78))
    write.csv(full$implied, paste0("output/sem/", key, "__implied.csv"))
    write.csv(input$correlation - full$implied, paste0("output/sem/", key, "__residuals.csv"))

    checkpoint_dir <- file.path("output/sem/replicates", key)
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    replicates <- parallel::mclapply(1:200, function(block) {
      checkpoint <- file.path(checkpoint_dir, sprintf("block_%03d.rds", block))
      if (file.exists(checkpoint)) return(readRDS(checkpoint))
      result <- list(estimates = point * NA_real_, diagnostics = tibble(
        valid = FALSE, boundary = FALSE, message = "Skipped: invalid full-data fit"
      ))
      if (full$diagnostics$valid) {
        r <- input$correlation
        r[pairs] <- r[pairs[, 2:1]] <- input$jackknife[block, ]
        result <- tryCatch({
          best <- fit_model(r, input, structure, list(unname(full$parameters)))
          if (!best$diagnostics$valid || best$diagnostics$boundary) {
            set.seed(4200 + block)
            starts <- c(list(unname(best$parameters)), make_starts(r, structure, parTable(full$fit)))
            best <- fit_model(r, input, structure, starts)
          }
          best$estimates <- align_signs(best$estimates, point)
          best
        }, error = \(e) list(estimates = point * NA_real_, diagnostics = tibble(
          valid = FALSE, boundary = FALSE, message = conditionMessage(e)
        )))
      }
      saveRDS(result, checkpoint)
      if (block %% 20 == 0) message(key, ": ", block, "/200 deletions")
      result
    }, mc.cores = cores, mc.preschedule = FALSE)
    stopifnot(!any(map_lgl(replicates, inherits, "try-error")))
    diagnostics <- map(replicates, "diagnostics") |>
      list_rbind() |>
      mutate(block = 1:200, .before = 1)
    deleted <- do.call(rbind, map(replicates, "estimates"))
    complete <- full$diagnostics$valid && all(diagnostics$valid) && all(is.finite(deleted))
    se <- rep(NA_real_, length(point))
    if (complete) {
      centered <- sweep(deleted, 2, colMeans(deleted), "-")
      covariance <- (199 / 200) * crossprod(centered)
      se <- sqrt(diag(covariance))
      write.csv(covariance, paste0("output/sem/", key, "__estimand_covariance.csv"))
    }
    estimates[[key]] <- tibble(
      analysis = analysis, structure = structure, quantity = names(point), estimate = unname(point),
      standard_error = se, valid_replicates = sum(diagnostics$valid),
      boundary_replicates = sum(diagnostics$boundary),
      se_status = if (!complete) "unavailable_invalid_fit_or_replicate" else if (
        full$diagnostics$boundary || any(diagnostics$boundary)
      ) "complete_200_blocks_with_boundaries" else "complete_200_blocks"
    )
    write_csv(diagnostics, paste0("output/sem/", key, "__jackknife_diagnostics.csv"))
    write_csv(as_tibble(deleted) |> mutate(block = 1:200, .before = 1),
              paste0("output/sem/", key, "__jackknife_estimates.csv"))
    message(key, ": ", sum(diagnostics$valid), "/200 valid; SEs ", if (complete) "reported" else "withheld")
  }
}

# ==============================================================================
# Paper tables and model comparisons
# ==============================================================================

fit_summary <- list_rbind(fit_summary)
estimates <- list_rbind(estimates)
comparison <- map_dfr(names(inputs), \(analysis) {
  x <- filter(fit_summary, .data$analysis == .env$analysis)
  delta_q <- x$q[x$structure == "hierarchical"] - x$q[x$structure == "correlated"]
  delta_df <- x$df[x$structure == "hierarchical"] - x$df[x$structure == "correlated"]
  regular <- all(x$valid & !x$boundary) && delta_q >= 0
  tibble(analysis = analysis, delta_q = delta_q, delta_df = delta_df,
         p_value = if (regular) pchisq(delta_q, delta_df, lower.tail = FALSE) else NA_real_)
})
groups <- estimates |>
  filter(quantity %in% c("mean_big_five", "mean_iq_big_five", "difference_iq_minus_big_five")) |>
  mutate(lower = estimate - qnorm(.975) * standard_error,
         upper = estimate + qnorm(.975) * standard_error,
         p_value = if_else(quantity == "difference_iq_minus_big_five" & standard_error > 0,
                           2 * pnorm(-abs(estimate / standard_error)), NA_real_))
write_csv(fit_summary, "output/sem/fit_summary.csv")
write_csv(estimates, "output/sem/estimates.csv")
write_csv(comparison, "output/sem/model_comparison.csv")
write_csv(groups, "output/sem/group_comparison.csv")
write_csv(tibble(path = unname(paths), md5 = unname(tools::md5sum(paths))), "output/sem/sources.csv")
write_xlsx(list(fit = fit_summary, estimates = estimates, group_comparison = groups,
                model_comparison = comparison), "output/sem/sem.xlsx")
writeLines(capture.output(sessionInfo()), "output/sem/session_info.txt")
print(fit_summary)
print(groups)
