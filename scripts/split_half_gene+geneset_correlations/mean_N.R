library(data.table)

munged_dir <- "munged"

files <- list.files(munged_dir, pattern = "\\.sumstats(\\.gz)?$", full.names = TRUE)

results <- rbindlist(lapply(files, function(f) {
  d <- fread(f, select = "N")
  data.table(
    trait   = sub("\\.sumstats(\\.gz)?$", "", basename(f)),
    n_snps  = nrow(d),
    mean_N  = mean(d$N, na.rm = TRUE),
    median_N = median(d$N, na.rm = TRUE),
    min_N   = min(d$N, na.rm = TRUE),
    max_N   = max(d$N, na.rm = TRUE)
  )
}))

results[, mean_N := round(mean_N, 1)]
print(results)

fwrite(results, "mean_N_by_trait.csv")