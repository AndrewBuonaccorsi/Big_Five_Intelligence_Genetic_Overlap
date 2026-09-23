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
library("data.table", lib = "~/R/libs")
library(mvtnorm)
source("prepare_functions.R")
source("gcor_PB_adjust.R")

setwd("/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/R_scripts")
big_five = c("open","consc","extra","agree","neurot")

if(T) {
for (trait in big_five) {
  
  setwd("/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/R_scripts")
  df = data.frame(trait1 = character(), trait2 = character(), mean_shift = numeric(), adj_rg = numeric(), se = numeric())
  
  for (n in 1:401) {
    meanshift = ((n-201)/100)
    res_gcor <- gcor_PB_adjust(path = '/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/',
                               mean_shift1 = meanshift, mean_shift2 = 0.44,
                               trait_name1 = trait, trait_name2 = 'ea', shift_sd1 = 0)
    
    newrow = data.frame(trait1 = 'ea',
                        trait2 = trait,
                        mean_shift = meanshift,
                        adj_rg = res_gcor$phig_adj,
                        se = res_gcor$phig_se_ajd)
    df = rbind(df, newrow)
  }
  setwd("/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/output")
  write.table(df, file = paste0(trait,"_ea_adj_0.44.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
}
}


levels = c(-1.5,-1,-0.5,-0.44,-0.15,0,0.15,0.44,0.5,1,1.5)

for (shift in levels) {
for (trait in big_five) {
  
  setwd("/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/R_scripts")
  df = data.frame(trait1 = character(), trait2 = character(), mean_shift = numeric(), adj_rg = numeric(), se = numeric())
  
  for (n in 1:401) {
    meanshift = ((n-201)/100)
    res_gcor <- gcor_PB_adjust(path = '/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/',
                               mean_shift1 = meanshift, mean_shift2 = shift,
                               trait_name1 = trait, trait_name2 = 'iq',shift_sd1 = 0, shift_sd2 = 0)
    
    newrow = data.frame(trait1 = 'iq',

                        trait2 = trait,

                        mean_shift = meanshift,

                        adj_rg = res_gcor$phig_adj,

                        se = res_gcor$phig_se_ajd)
    df = rbind(df, newrow)
    
  }
  
  setwd("/users/0/buona008/andrew/personality_iq_project/participation_bias_correction/output")
  write.table(df, file = paste0(trait,"_iq_adj_", as.character(shift), ".txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  }
}
