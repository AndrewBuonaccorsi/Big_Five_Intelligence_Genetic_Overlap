library(dplyr)



block <- Sys.getenv("BLOCK")

group <- Sys.getenv("GROUP")

setwd("/users/0/buona008/andrew/personality_iq_project/split_half/temp_dir")

trait_files <- list.files(pattern = "\\.gsa\\.out$")

short_names <- sub("\\..*$", "", trait_files)

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

cat("Retained", nrow(merged), "gene sets for group:", group, "\n")

beta_cols <- grep("^BETA_STD_", colnames(merged), value = TRUE)

se_cols   <- grep("^SE_", colnames(merged), value = TRUE)



beta_mat <- merged[, c("FULL_NAME", beta_cols)]

se_mat   <- merged[, c("FULL_NAME", se_cols)]



colnames(beta_mat) <- c("FULL_NAME", short_names)

colnames(se_mat)   <- c("FULL_NAME", short_names)



valid_idx <- complete.cases(beta_mat, se_mat)

beta_mat <- beta_mat[valid_idx, ]

se_mat   <- se_mat[valid_idx, ]

pearson_corr <- function(bi, bj) {

  ok <- is.finite(bi) & is.finite(bj)

  bi <- bi[ok]

  bj <- bj[ok]



  n <- length(bi)

  if (n < 5) return(NA)



  cor(bi, bj, method = "pearson")

}

n_traits <- length(short_names)

cor_vals <- c()

col_names <- c()



for (i in 1:n_traits) {

  for (j in i:n_traits) {

    bi <- beta_mat[[short_names[i]]]

    bj <- beta_mat[[short_names[j]]]



    cv <- pearson_corr(bi, bj)



    col_names <- c(col_names, paste(short_names[i], short_names[j], sep = "_vs_"))

    cor_vals  <- c(cor_vals, cv)

  }

}



cor_df <- as.data.frame(t(cor_vals))

colnames(cor_df) <- col_names

cor_df <- cbind(GROUP = group, BLOCK = block, cor_df)



print(col_names)

print(cor_vals)



setwd("/users/0/buona008/andrew/personality_iq_project/split_half/unweighted_geneset_corr")

out_file <- paste0(group, "_unweighted.txt")



write.table(

  cor_df,

  file = out_file,

  append = TRUE,

  row.names = FALSE,

  col.names = !file.exists(out_file),

  sep = "\t",

  quote = FALSE

)
