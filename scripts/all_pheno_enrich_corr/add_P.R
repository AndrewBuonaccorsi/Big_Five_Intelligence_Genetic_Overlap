#load packages
library(vroom)

files <- list.files("phenotypes", pattern = "\\.sumstats$", full.names = TRUE)

for (filename in files) {
sumstats <- vroom(filename)

sumstats$P <- 2 * pnorm(-abs(sumstats$Z))

vroom_write(sumstats, filename)
}
