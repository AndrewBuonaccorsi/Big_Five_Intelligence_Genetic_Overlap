
traits <- list(
  Openness          = list(es = c(0.48,  0.19,  0.37,  0.17), n = c(646, 424, 658, 1140)),
  Neuroticism       = list(es = c(-0.18, -0.14, -0.18, -0.20), n = c(646, 424, 658, 1140)),
  Extraversion      = list(es = c(0.10, -0.36, -0.08),         n = c(424, 658, 1140)),
  Conscientiousness = list(es = c(0.00, -0.15,  0.16),         n = c(424, 658, 1140)),
  Agreeableness     = list(es = c(0.42, -0.12, -0.14, -0.08),  n = c(646, 424, 658, 1140))
)

results <- do.call(rbind, lapply(names(traits), function(trait) {
  es <- traits[[trait]]$es
  n <- traits[[trait]]$n
  stopifnot(length(es) == length(n))
  data.frame(
    Trait = trait,
    k = length(es),
    Total_N = sum(n),
    Weighted_rg = sum(es * n) / sum(n),
    stringsAsFactors = FALSE
  )
}))

results <- results[order(-abs(results$Weighted_rg)), ]
results$Weighted_rg <- round(results$Weighted_rg, 3)
rownames(results) <- NULL

print(results, row.names = FALSE)