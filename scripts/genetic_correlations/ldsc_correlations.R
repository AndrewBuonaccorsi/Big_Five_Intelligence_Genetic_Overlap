# install.packages("devtools")
# devtools::install_github("GenomicSEM/GenomicSEM")

suppressMessages(library(GenomicSEM))

munged_dir <- "munged"          # folder holding the *.sumstats.gz files
ld_dir     <- "/users/0/buona008/andrew/reference_files/eur_w_ld_chr"    # LD scores directory (EUR); change if not EUR
wld_dir    <- ld_dir            # regression weights (same dir for standard EUR)
out_prefix <- "ldsc_rg"         # prefix for all output files

traits <- sort(list.files(
  path       = munged_dir,
  pattern    = "\\.sumstats$",
  full.names = TRUE
))

if (length(traits) < 2L) {
  stop("Need >= 2 '.sumstats.gz' files in '", munged_dir,
       "' to estimate genetic correlations (found ", length(traits), ").")
}
if (!dir.exists(ld_dir)) {
  stop("LD score directory not found: '", ld_dir, "'. See notes at top of script.")
}

trait.names <- sub("\\.sumstats$", "", basename(traits))

message("Found ", length(traits), " munged files:")
message(paste0("  - ", trait.names, collapse = "\n"))

sample.prev     <- rep(NA, length(traits))
population.prev <- rep(NA, length(traits))

LDSCoutput <- ldsc(
  traits          = traits,
  sample.prev     = sample.prev,
  population.prev = population.prev,
  ld              = ld_dir,
  wld             = wld_dir,
  trait.names     = trait.names,
  stand           = TRUE,          # also return S_Stand / V_Stand (the rg matrix)
  ldsc.log        = out_prefix     # writes <out_prefix>_ldsc.log
)

save(LDSCoutput, file = paste0(out_prefix, ".RData"))

k <- nrow(LDSCoutput$S)

SE_cov <- matrix(0, k, k)
SE_cov[lower.tri(SE_cov, diag = TRUE)] <- sqrt(diag(LDSCoutput$V))

h2    <- diag(LDSCoutput$S)
h2_se <- diag(SE_cov)

h2_tab <- data.frame(
  trait = trait.names,
  h2    = round(h2, 4),
  se    = round(h2_se, 4),
  Z     = round(h2 / h2_se, 3),
  stringsAsFactors = FALSE
)
write.csv(h2_tab, paste0(out_prefix, "_h2.csv"), row.names = FALSE)

rg_mat <- LDSCoutput$S_Stand                 # == cov2cor(LDSCoutput$S)
dimnames(rg_mat) <- list(trait.names, trait.names)

SE_stand <- matrix(0, k, k)
SE_stand[lower.tri(SE_stand, diag = TRUE)] <- sqrt(diag(LDSCoutput$V_Stand))
SE_stand[upper.tri(SE_stand)] <- t(SE_stand)[upper.tri(SE_stand)]
dimnames(SE_stand) <- list(trait.names, trait.names)

write.csv(round(rg_mat, 4),   paste0(out_prefix, "_rg_matrix.csv"))
write.csv(round(SE_stand, 4), paste0(out_prefix, "_rg_se_matrix.csv"))

pairs <- t(utils::combn(k, 2))
rg_long <- data.frame(
  trait1 = trait.names[pairs[, 1]],
  trait2 = trait.names[pairs[, 2]],
  rg     = rg_mat[pairs],
  se     = SE_stand[pairs],
  stringsAsFactors = FALSE
)
rg_long$Z <- rg_long$rg / rg_long$se
rg_long$p <- 2 * pnorm(-abs(rg_long$Z))          # tests rg = 0 (not rg = 1)
rg_long$rg <- round(rg_long$rg, 4)
rg_long$se <- round(rg_long$se, 4)
rg_long$Z  <- round(rg_long$Z, 3)
rg_long$p  <- signif(rg_long$p, 3)
rg_long <- rg_long[order(-abs(rg_long$Z)), ]
write.csv(rg_long, paste0(out_prefix, "_rg_pairs.csv"), row.names = FALSE)

I_mat <- LDSCoutput$I
dimnames(I_mat) <- list(trait.names, trait.names)
write.csv(round(I_mat, 4), paste0(out_prefix, "_intercepts.csv"))

print(h2_tab, row.names = FALSE)

print(round(rg_mat, 3))

print(rg_long, row.names = FALSE)

cat("\nOutputs written with prefix '", out_prefix, "':\n", sep = "")
cat("  ", out_prefix, ".RData           full LDSCoutput object\n", sep = "")
cat("  ", out_prefix, "_rg_matrix.csv   rg matrix\n", sep = "")
cat("  ", out_prefix, "_rg_se_matrix.csv rg standard errors\n", sep = "")
cat("  ", out_prefix, "_rg_pairs.csv    tidy pairwise rg / SE / Z / p\n", sep = "")
cat("  ", out_prefix, "_h2.csv          SNP heritabilities\n", sep = "")
cat("  ", out_prefix, "_intercepts.csv  LDSC intercept matrix (overlap check)\n", sep = "")
cat("  ", out_prefix, "_ldsc.log        run log\n", sep = "")
