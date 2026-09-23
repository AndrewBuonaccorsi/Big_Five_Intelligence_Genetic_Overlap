lines <- readLines("snps.genes.annot")

lines <- lines[!grepl("^#", lines)]

split_lines <- strsplit(lines, "\\s+")

gene <- sapply(split_lines, `[`, 1)

coord_field <- sapply(split_lines, `[`, 2)

coords <- strsplit(coord_field, ":")

out <- data.frame(
  GENE  = gene,
  CHR   = sapply(coords, `[`, 1),
  START = sapply(coords, `[`, 2),
  END   = sapply(coords, `[`, 3),
  stringsAsFactors = FALSE
)

write.table(
  out,
  file = "gene_coords.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)