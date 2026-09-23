
setwd("/projects/standard/leej5/edwa0506/andrew_project/gsa")
library(pacman)

p_load(tidyverse,
       purrr)

#library(readr)
#library(dplyr)
#library(stringr)
#library(purrr)

df <- read_delim("figures/RELIG_vs_ATTEND.csv", delim = "\t", show_col_types = FALSE, progress = FALSE, trim_ws = TRUE)

# --------------------------------------------------
# Add non-null SNP counts (pi * M) with 95% CIs
# --------------------------------------------------

library(dplyr)

# Your total number of SNPs
M <- 11956373
df2 <- df %>%
  mutate(
    # Trait 1 (RELIG): pi1 * M
    nn1_mean  = `pi1 (mean)` * M,
    nn1_ci_lo = pmax(0, (`pi1 (mean)` - 1.96 * `pi1 (std)`) * M),
    nn1_ci_hi = (`pi1 (mean)` + 1.96 * `pi1 (std)`) * M,
    
    # Trait 2 (ATTEND): pi2 * M
    nn2_mean  = `pi2 (mean)` * M,
    nn2_ci_lo = pmax(0, (`pi2 (mean)` - 1.96 * `pi2 (std)`) * M),
    nn2_ci_hi = (`pi2 (mean)` + 1.96 * `pi2 (std)`) * M,
    
    # Shared component: pi12 * M
    nn12_mean  = `pi12 (mean)` * M,
    nn12_ci_lo = pmax(0, (`pi12 (mean)` - 1.96 * `pi12 (std)`) * M),
    nn12_ci_hi = (`pi12 (mean)` + 1.96 * `pi12 (std)`) * M
  )

dir.create("figures", showWarnings = FALSE)
write.csv(df2, "figures/cleaned_relig_attend.csv", row.names = FALSE)


# personality
pairs <- c(
  "OPEN_vs_IQ","CONSC_vs_IQ","EXTRA_vs_IQ","AGREE_vs_IQ","NEURO_vs_IQ",
  "OPEN_vs_HEIGHT","CONSC_vs_HEIGHT","EXTRA_vs_HEIGHT","AGREE_vs_HEIGHT","NEURO_vs_HEIGHT",
  "IQ_vs_HEIGHT",
  "OPEN_vs_CONSC","OPEN_vs_EXTRA","OPEN_vs_AGREE","OPEN_vs_NEURO",
  "CONSC_vs_EXTRA","CONSC_vs_AGREE","CONSC_vs_NEURO",
  "EXTRA_vs_AGREE","EXTRA_vs_NEURO",
  "AGREE_vs_NEURO"
)


read_mixer_pair <- function(path) {
  tryCatch({
    df <- suppressWarnings(
      read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE, trim_ws = TRUE)
    )
    
    # If it collapsed to 1 column, manually split
    if (ncol(df) == 1) {
      x <- df[[1]]
      hdr <- str_split(x[1], "\t", simplify = TRUE)
      dat <- str_split(x[-1], "\t", simplify = TRUE)
      
      df <- as_tibble(dat, .name_repair = "minimal")
      names(df) <- as.character(hdr)
    }
    
    df
  }, error = function(e) {
    message("Skipping unreadable file: ", path)
    NULL
  })
}
all_df <- map_dfr(pairs, \(p) {
  f <- file.path("figures", paste0(p, ".csv"))
  if (!file.exists(f)) return(NULL)
  
  read_mixer_pair(f) %>%
    mutate(pair = p, .before = 1)
})

# (Optional) keep only fit rows if you don't want duplicate fit/test
all_df <- all_df %>% 
              filter(str_detect(fname, "\\.fit\\.json$")) %>%
              drop_na()

# TE I prefer the Jaccard coefficients over the DICE
# Since DICE is what is recommended by the authors it is perhaps best to keep to it
# Also to get "right" standard errors for Jaccard we would need to alter the python script
# instead of using the delta theorem as I have used below. 

all_df <- all_df %>%
  mutate(
    dice    = as.numeric(as.character(`dice (mean)`)),
    dice_se = as.numeric(as.character(`dice (std)`)),
    
    # transform
    jaccard = dice / (2 - dice),
    
    # delta-method SE: g'(D) = 2 / (2-D)^2
    jaccard_se = abs(2 / (2 - dice)^2) * dice_se,
    
    # CLT CI on J directly (not by transforming endpoints)
    jaccard_low  = jaccard - 1.96 * jaccard_se,
    jaccard_high = jaccard + 1.96 * jaccard_se
  )

all_df %>% 
  write_csv("figures/personality_bivariate.csv")


# Use a fixed trait set (your desired rows/cols)
traits <- c("OPEN","CONSC","EXTRA","AGREE","NEURO","IQ","HEIGHT")

# Map from the names in your MiXeR CSVs -> short trait labels
label_map <- c(
  open_schwaba2025   = "OPEN",
  consc_schwaba2025  = "CONSC",
  extra_schwaba2025  = "EXTRA",
  agree_schwaba2025  = "AGREE",
  neurot_schwaba2025 = "NEURO",
  iq_savage2018      = "IQ",
  height_yengo2022   = "HEIGHT"
)

# Pull dice means from FIT rows, label traits, keep only your target set
# Pull dice mean + SE from FIT rows, label traits, keep only your target set
dice_edges <- all_df %>%
  filter(str_detect(fname, "\\.fit\\.json$")) %>%
  transmute(
    t1      = recode(trait1, !!!label_map),
    t2      = recode(trait2, !!!label_map),
    dice    = as.numeric(as.character(`dice (mean)`)),
    dice_se = as.numeric(as.character(`dice (std)`))
  ) %>%
  filter(!is.na(t1), !is.na(t2), !is.na(dice), !is.na(dice_se)) %>%
  group_by(t1, t2) %>%                 # guard duplicates
  summarise(
    dice    = mean(dice),
    dice_se = mean(dice_se),
    .groups = "drop"
  )

# Initialize matrix
dice_mat <- matrix(NA_real_,
                   nrow = length(traits), ncol = length(traits),
                   dimnames = list(traits, traits)
)
diag(dice_mat) <- 1  # keep diagonal as 1s

# Fill: mean in upper triangle, SE in lower triangle
for (k in seq_len(nrow(dice_edges))) {
  a  <- dice_edges$t1[k]
  b  <- dice_edges$t2[k]
  ia <- match(a, traits)
  ib <- match(b, traits)
  if (is.na(ia) || is.na(ib) || ia == ib) next
  
  m  <- dice_edges$dice[k]
  se <- dice_edges$dice_se[k]
  
  if (ia < ib) {
    dice_mat[a, b] <- m   # upper
    dice_mat[b, a] <- se  # lower
  } else {
    dice_mat[b, a] <- m
    dice_mat[a, b] <- se
  }
}

dice_mat

data.frame(
  trait = colnames(dice_mat),
  dice_mat,
  row.names = NULL
) |>
  write_csv("figures/personality_bivariate_cormat.csv")

big5 <- c("OPEN","CONSC","EXTRA","AGREE","NEURO")

avg_dice_big5 <- mean(dice_mat[big5, big5][upper.tri(dice_mat[big5, big5])])
avg_dice_iq_big5 <- mean(dice_mat[big5, "IQ"])

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

out <- data.frame(
  contrast = c("BigFive_within", "IQ_vs_BigFive"),
  mean_dice = c(avg_dice_big5, avg_dice_iq_big5)
)

write.csv(out, file = "figures/dice_iq_is_a_personality_trait.csv", row.names = FALSE)










