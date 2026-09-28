#!/usr/bin/env Rscript

# Load model definitions without running the full-data or 200-block analyses.
expressions <- parse("scripts/04_sem.R")
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("paths"))) break
  eval(expr)
}
for (expr in parse("scripts/03_analytic_correlations.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]]) %in% c("analytic_correlations", "group_correlations")) eval(expr)
}
pair_names <- paste(traits[trait_pairs[, 1]], traits[trait_pairs[, 2]], sep = "__")
iq_pair <- trait_pairs[, 2] == 6
pacman::p_load(numDeriv)

set.seed(842)
for (structure in names(models)) {
  loading <- seq(.3, .6, length.out = 12)
  method <- seq(.08, .25, length.out = 10)
  general <- seq(.4, .7, length.out = 6)
  phi <- outer(general, general)
  if (structure == "correlated") {
    perturbation <- matrix(rnorm(36, sd = .2), 6)
    phi <- cov2cor(diag(6) + tcrossprod(perturbation) + outer(general, general))
  }
  diag(phi) <- 1
  observed <- outer(loading, loading) * phi[indicator_trait, indicator_trait]
  for (indices in list(seq(1, 9, 2), seq(2, 10, 2))) {
    observed[indices, indices] <- observed[indices, indices] + outer(method[indices], method[indices])
  }
  diag(observed) <- 1
  dimnames(observed) <- list(indicators, indicators)
  noise <- matrix(rnorm(66 * 66), 66)
  variance <- (diag(66) + crossprod(noise) / 66) * 1e-4
  input <- list(vcov = variance, precision = solve(variance))

  # Unequal paired loadings and arbitrary same-half sample effects must not
  # change the population analytic cross-half correlation.
  analytic <- analytic_correlations(observed)
  stopifnot(max(abs(analytic$estimate - phi[trait_pairs])) < 1e-12)
  invalid <- observed
  invalid[1, 4] <- invalid[4, 1] <- -invalid[1, 4]
  stopifnot(is.na(analytic_correlations(invalid)$estimate[1]))
  unbounded <- observed
  unbounded[1, 4] <- unbounded[4, 1] <- .9
  unbounded[2, 3] <- unbounded[3, 2] <- .9
  stopifnot(analytic_correlations(unbounded)$estimate[1] > 1,
            analytic_correlations(unbounded)$status[1] == "outside_correlation_bounds")

  template <- fit_once(observed, input, structure, do_fit = FALSE)
  latent_sd <- if (structure == "hierarchical") sqrt(1 - general^2) else rep(1, 6)
  values <- c(setNames(loading * latent_sd[indicator_trait], paste0("l_", indicators)),
              setNames(method, paste0("m_", indicators[1:10])))
  values <- if (structure == "hierarchical") {
    c(values, setNames(general / latent_sd, paste0("g_", traits)))
  } else {
    c(values, setNames(phi[trait_pairs], paste0("p_", traits[trait_pairs[, 1]], "_", traits[trait_pairs[, 2]])))
  }
  start <- unname(values[template$label[template$free > 0]])
  fit <- fit_once(observed, input, structure, start)
  stopifnot(fit$diagnostics$valid, !fit$diagnostics$boundary,
            fit$diagnostics$parameters == if (structure == "correlated") 37 else 28,
            max(abs(fit$implied - observed)) < 1e-7,
            sum(template$op == "<") == 12,
            all(template$free[template$op == "=~"] > 0))
  stopifnot(!any(template$lhs %in% c("sample_1", "sample_2") &
                  template$op == "=~" & template$rhs %in% indicators[11:12]))

  # Force an observed negative-uniqueness solution and verify enforcement of the
  # zero-variance boundary, including an independently calculated KKT gradient.
  half <- seq(2, 10, 2)
  new_method <- method[half]
  new_method[5] <- 1.1
  boundary_data <- observed
  boundary_data[half, half] <- boundary_data[half, half] - outer(method[half], method[half]) +
    outer(new_method, new_method)
  diag(boundary_data) <- 1
  constrained <- fit_once(boundary_data, input, structure, start)
  stopifnot(constrained$diagnostics$converged, constrained$diagnostics$kkt,
            constrained$diagnostics$boundary, constrained$diagnostics$min_residual >= -residual_tolerance)
  residual <- boundary_data[pairs] - constrained$implied[pairs]
  stopifnot(abs(unname(fitMeasures(constrained$fit, "chisq")) -
                  drop(crossprod(residual, input$precision %*% residual))) < 1e-6)
  objective <- function(p) {
    latent <- diag(6)
    if (structure == "hierarchical") {
      g <- p[paste0("g_", traits)]
      latent <- latent + outer(g, g)
    } else {
      latent[trait_pairs] <- latent[trait_pairs[, 2:1]] <- p[
        paste0("p_", traits[trait_pairs[, 1]], "_", traits[trait_pairs[, 2]])
      ]
    }
    l <- p[paste0("l_", indicators)]
    implied <- outer(l, l) * latent[indicator_trait, indicator_trait]
    for (indices in list(seq(1, 9, 2), seq(2, 10, 2))) {
      m <- p[paste0("m_", indicators[indices])]
      implied[indices, indices] <- implied[indices, indices] + outer(m, m)
    }
    error <- boundary_data[pairs] - implied[pairs]
    drop(crossprod(error, input$precision %*% error)) / 398
  }
  delta <- lavInspect(constrained$fit, "delta")
  gradient <- drop(-crossprod(delta, input$precision %*% residual) / 199)
  stopifnot(max(abs(numDeriv::grad(objective, constrained$parameters) - gradient)) < 1e-6)

  reference <- fit$estimates
  reversed <- reference
  keys <- paste0("sample_loading__", indicators[half])
  reversed[keys] <- -reversed[keys]
  stopifnot(max(abs(align_signs(reversed, reference) - reference)) < 1e-12)
}
message("PASS: analytic recovery, unequal loadings, both SEMs, WLS scaling, residual constraints and sign alignment.")
