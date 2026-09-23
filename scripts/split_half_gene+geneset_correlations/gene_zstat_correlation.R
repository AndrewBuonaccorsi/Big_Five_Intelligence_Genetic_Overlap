args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript gene_zstat_correlation.R <file_list.txt> <cor_output_prefix> <dice_output_prefix>")
}

file_list_path     <- args[1]
cor_output_prefix  <- args[2]
dice_output_prefix <- args[3]

files <- readLines(file_list_path)


read_gene_file <- function(f, exclude_genes = NULL) {
  dat <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
  
  if (!all(c("GENE", "ZSTAT") %in% colnames(dat))) {
    stop(paste("Missing required columns in file:", f))
  }
  
  if (!is.null(exclude_genes)) {
    dat <- dat[!(dat$GENE %in% exclude_genes), ]
  }
  
  out <- dat[, c("GENE", "ZSTAT")]
  colnames(out)[2] <- tools::file_path_sans_ext(basename(f))
  
  return(out)
}

get_sig_genes <- function(f, exclude_genes = NULL) {
  dat <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
  
  if (!all(c("GENE", "ZSTAT") %in% colnames(dat))) {
    stop(paste("Missing required columns in file:", f))
  }
  
  if (!is.null(exclude_genes)) {
    dat <- dat[!(dat$GENE %in% exclude_genes), ]
  }
  
  sig_genes <- dat$GENE[dat$ZSTAT > 1]
  return(unique(sig_genes))
}

exclude_files <- list.files(
  path = "excluded_genes",
  pattern = "^exclude_genes_[0-9]+\\.txt$",
  full.names = TRUE
)

if (length(exclude_files) == 0) {
  stop("No exclusion files found in excluded_genes/")
}

exclude_files <- sort(exclude_files)
n_iter <- length(exclude_files)

cat("Found", n_iter, "jackknife exclusion files.\n")

trait_names <- tools::file_path_sans_ext(basename(files))
n_traits <- length(files)

cor_list  <- vector("list", n_iter)
dice_list <- vector("list", n_iter)

for (k in seq_along(exclude_files)) {
  
  cat("Running jackknife iteration", k, "of", n_iter, "\n")
  
  exclude_genes <- readLines(exclude_files[k])
  
  data_list <- lapply(files, read_gene_file, exclude_genes = exclude_genes)
  
  merged_data <- Reduce(function(x, y)
    merge(x, y, by = "GENE", all = FALSE),
    data_list)
  
  z_matrix <- merged_data[, -1]
  
  cor_matrix <- cor(z_matrix,
                    use = "pairwise.complete.obs",
                    method = "pearson")
  
  cor_list[[k]] <- cor_matrix
  
  gene_lists <- lapply(files, get_sig_genes, exclude_genes = exclude_genes)
  names(gene_lists) <- trait_names
  
  dice_matrix <- matrix(0, nrow = n_traits, ncol = n_traits)
  rownames(dice_matrix) <- trait_names
  colnames(dice_matrix) <- trait_names
  
  for (i in 1:n_traits) {
    for (j in i:n_traits) {
      
      A <- gene_lists[[i]]
      B <- gene_lists[[j]]
      
      intersection_size <- length(intersect(A, B))
      size_A <- length(A)
      size_B <- length(B)
      
      if ((size_A + size_B) == 0) {
        dice_value <- NA
      } else {
        dice_value <- (2 * intersection_size) / (size_A + size_B)
      }
      
      dice_matrix[i, j] <- dice_value
      dice_matrix[j, i] <- dice_value
    }
  }
  
  dice_list[[k]] <- dice_matrix
}

cor_array  <- simplify2array(cor_list)
dice_array <- simplify2array(dice_list)

cor_mean  <- apply(cor_array,  c(1,2), mean, na.rm = TRUE)
dice_mean <- apply(dice_array, c(1,2), mean, na.rm = TRUE)

cor_var <- apply(cor_array, c(1,2), function(x) {
  m <- mean(x, na.rm = TRUE)
  ((n_iter - 1) / n_iter) * sum((x - m)^2, na.rm = TRUE)
})

dice_var <- apply(dice_array, c(1,2), function(x) {
  m <- mean(x, na.rm = TRUE)
  ((n_iter - 1) / n_iter) * sum((x - m)^2, na.rm = TRUE)
})

cor_se  <- sqrt(cor_var)
dice_se <- sqrt(dice_var)

get_upper <- function(mat) {
  mat[upper.tri(mat)]
}

n_pairs <- n_traits * (n_traits - 1) / 2

cor_theta  <- matrix(NA, nrow = n_iter, ncol = n_pairs)
dice_theta <- matrix(NA, nrow = n_iter, ncol = n_pairs)

for (k in 1:n_iter) {
  cor_theta[k, ]  <- get_upper(cor_array[,,k])
  dice_theta[k, ] <- get_upper(dice_array[,,k])
}

cor_cov  <- (n_iter - 1) * cov(cor_theta)
dice_cov <- (n_iter - 1) * cov(dice_theta)

pair_names <- combn(trait_names, 2,
                    FUN = function(x) paste(x[1], x[2], sep = "_"))

colnames(cor_cov)  <- pair_names
rownames(cor_cov)  <- pair_names
colnames(dice_cov) <- pair_names
rownames(dice_cov) <- pair_names

colnames(cor_theta)  <- pair_names
colnames(dice_theta) <- pair_names

rownames(cor_theta)  <- paste0("JK_", seq_len(n_iter))
rownames(dice_theta) <- paste0("JK_", seq_len(n_iter))

write.table(cor_theta,
            file = paste0(cor_output_prefix, "_jackknife_estimates.txt"),
            quote = FALSE,
            sep = "\t",
            col.names = NA)

write.table(cor_mean,
            file = paste0(cor_output_prefix, "_mean.txt"),
            quote = FALSE,
            sep = "\t",
            col.names = NA)

write.table(cor_se,
            file = paste0(cor_output_prefix, "_jackknife_se.txt"),
            quote = FALSE,
            sep = "\t",
            col.names = NA)

write.table(cor_cov,
            file = paste0(cor_output_prefix, "_jackknife_cov.txt"),
            quote = FALSE,
            sep = "\t")

cat("Jackknife analysis complete (mean, SE, full covariance matrices).\n")
