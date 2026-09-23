args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript gene_zstat_correlation.R <file_list.txt> <cor_output.txt> <dice_output.txt>")
}

file_list_path <- args[1]
cor_output     <- args[2]
dice_output    <- args[3]

files <- readLines(file_list_path)


read_gene_file <- function(f) {
  dat <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
  
  if (!all(c("GENE", "ZSTAT") %in% colnames(dat))) {
    stop(paste("Missing required columns in file:", f))
  }
  
  out <- dat[, c("GENE", "ZSTAT")]
  colnames(out)[2] <- tools::file_path_sans_ext(basename(f))
  
  return(out)
}

get_sig_genes <- function(f) {
  dat <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
  
  if (!all(c("GENE", "ZSTAT") %in% colnames(dat))) {
    stop(paste("Missing required columns in file:", f))
  }
  
  sig_genes <- dat$GENE[dat$ZSTAT > 1]
  return(unique(sig_genes))
}


trait_names <- tools::file_path_sans_ext(basename(files))
n_traits <- length(files)


data_list <- lapply(files, read_gene_file)

merged_data <- Reduce(function(x, y)
  merge(x, y, by = "GENE", all = FALSE),
  data_list)

z_matrix <- merged_data[, -1]

cor_matrix <- cor(
  z_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)

gene_lists <- lapply(files, get_sig_genes)
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

write.table(
  cor_matrix,
  file = cor_output,
  quote = FALSE,
  sep = "\t",
  col.names = NA
)

#write.table(
#  dice_matrix,
#  file = dice_output,
#  quote = FALSE,
#  sep = "\t",
#  col.names = NA
#)

cat("Analysis complete.\n")
cat("Correlation matrix written to:", cor_output, "\n")
cat("Dice matrix written to:", dice_output, "\n")
