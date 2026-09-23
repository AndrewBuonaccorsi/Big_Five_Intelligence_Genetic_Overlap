library(dplyr)



group     <- Sys.getenv("GROUP")

group_dir <- Sys.getenv("GROUP_DIR")

out_file  <- Sys.getenv("OUT_FILE")


trait_files <- list.files(group_dir, pattern = "\\.gsa\\.out$", full.names = TRUE)

short_names <- sub("\\..*$", "", basename(trait_files))



if (length(trait_files) == 0) {

  stop(sprintf("No .gsa.out files found in %s", group_dir))

}

traits_list <- lapply(seq_along(trait_files), function(k) {
  f  <- trait_files[k]
  sn <- short_names[k]
  df <- read.table(
    f,
    header = TRUE,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"FULL_NAME" %in% colnames(df)) {
    df$FULL_NAME <- df$VARIABLE
  }

  df %>%

    select(FULL_NAME, BETA_STD, SE) %>%

    rename(

      !!paste0("BETA_STD_", sn) := BETA_STD,

      !!paste0("SE_", sn)       := SE

    )

})

names(traits_list) <- short_names
merged <- Reduce(

  function(x, y) full_join(x, y, by = "FULL_NAME"),

  traits_list

)



merged$FULL_NAME <- trimws(tolower(merged$FULL_NAME))



cat("Group:", group, "-- retained", nrow(merged), "gene sets.\n")

beta_cols <- grep("^BETA_STD_", colnames(merged), value = TRUE)

se_cols   <- grep("^SE_", colnames(merged), value = TRUE)



beta_mat <- merged[, c("FULL_NAME", beta_cols)]

se_mat   <- merged[, c("FULL_NAME", se_cols)]



colnames(beta_mat) <- c("FULL_NAME", short_names)

colnames(se_mat)   <- c("FULL_NAME", short_names)



valid_idx <- complete.cases(beta_mat, se_mat)

beta_mat <- beta_mat[valid_idx, ]

se_mat   <- se_mat[valid_idx, ]

n_traits <- length(short_names)

cor_mat <- matrix(NA_real_, n_traits, n_traits,

                   dimnames = list(short_names, short_names))



for (i in 1:n_traits) {

  for (j in i:n_traits) {

    bi <- beta_mat[[short_names[i]]]

    bj <- beta_mat[[short_names[j]]]



    ok    <- is.finite(bi) & is.finite(bj)

    bi_ok <- bi[ok]

    bj_ok <- bj[ok]



    cv <- if (length(bi_ok) < 5) NA else cor(bi_ok, bj_ok, method = "pearson")



    cor_mat[i, j] <- cv

    cor_mat[j, i] <- cv

  }

}

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

write.csv(cor_mat, out_file)

