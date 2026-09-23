#Load packages
options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages("remotes", lib="~/R/libs")
#install.packages("MASS", lib = "~/R/libs")
#install.packages("lavaan", lib = "~/R/libs")
#install.packages("tzdb", lib = "~/R/libs")
#install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "forcats"), dependencies = TRUE)
#install.packages("vroom", lib = "~/R/libs")
#library("remotes", lib = "~/R/libs")
#remotes::install_github("GenomicSEM/GenomicSEM", lib = "~/R/libs")
#library("tidyverse", lib = "~/R/libs")
library("readr", lib = "~/R/libs")
library("vroom", lib = "~/R/libs")
library("tzdb", lib = "~/R/libs")
library("lavaan", lib = "~/R/libs")
library("GenomicSEM", lib = "~/R/libs")

load("personality-iq-correlations_mid.Rdata")

#creating variable names in each category
bf = sub("\\.sumstats\\.gz$", "", big_five)
phe = sub("\\.sumstats\\.gz$", "", phenotypes)
intel = sub("\\.sumstats\\.gz$", "", intelligence) 

#Instantiating the data frame
output_table <- data.frame(
  PREDICTOR = character(),
  PHENOTYPE = character(),
  BETA = numeric(),
  SE = numeric(),
  MULTIPLE_BETA = numeric(),
  MULTIPLE_SE = numeric(),
  stringsAsFactors = FALSE
)


#Define function
make_change_model <- function(x, y, z) {
    sprintf("
    %s ~ cprime*%s + b*%s
    # label moments among predictors
    %s ~~ vX*%s
    %s ~~ vZ*%s
    %s ~~ cXZ*%s
    # implied (unadjusted) effect of %s on %s
    beta_marg := cprime + b*(cXZ/vX)
    # change when adding %s (unadjusted - adjusted)
    delta := beta_marg - cprime
  ",
          y, x, z,
          x, x,
          z, z,
          z, x,
          x, y,
          z)
}

#Calculating betas and SEs for all combinations of Big Five domain and phenotype
for (x in bf) {
  for (y in phe) {
fit_change <- usermodel(LDSCoutput, estimation="DWLS", make_change_model(x, y, intel))
new_row = c(x,y,fit_change$results[1,6],fit_change$results[1,7],fit_change$results[1,9],fit_change$results[7,6],fit_change$results[7,7],fit_change$results[7,9],fit_change$results[8,6],fit_change$results[8,7],fit_change$results[8,9])

print(new_row)
print(fit_change$results)

output_table = rbind(output_table, new_row)
}}
  
write.table(output_table, file = "personality-gencorr-table.txt", sep = "\t", row.names = FALSE, quote = FALSE)
