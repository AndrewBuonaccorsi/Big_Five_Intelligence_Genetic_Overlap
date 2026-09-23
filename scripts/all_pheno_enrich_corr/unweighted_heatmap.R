library(dplyr)
library(pheatmap)
library(Cairo)
library(grid)

group     <- Sys.getenv("GROUP")
group_dir <- Sys.getenv("GROUP_DIR")
out_file  <- Sys.getenv("OUT_FILE")


# Nice names table
nice_names <- c(
  adhd = "ADHD",
  afb = "Age at First Birth",
  afi = "Age at First Intercourse",
  agemenarche = "Age at Menarche",
  agemenopause = "Age at Menopause",
  agesmokinit = "Age Smoking Initiation",
  agree = "Agreeableness",
  als = "Amyotrophic Lateral Sclerosis",
  alzheimer = "Alzheimer's Disease",
  anorexia = "Anorexia Nervosa",
  anxiety = "Anxiety",
  arthritis = "Rheumatoid Arthritis",
  asd = "Autism Spectrum Disorder",
  asthma = "Asthma",
  atopicdermatitis = "Atopic Dermatitis",
  bipolar = "Bipolar Disorder",
  bmi = "Body Mass Index",
  bodyfat = "Body Fat Percentage",
  cad = "Coronary Artery Disease",
  childless = "Childlessness",
  cigsperday = "Cigarettes Per Day",
  coffeeperday = "Coffee Per Day",
  covid = "Hospitalized with Covid",
  consc = "Conscientiousness",
  crohns = "Crohn's Disease",
  darkskin = "Dark Skin",
  dbp = "Diastolic Blood Pressure",
  depression = "Depression",
  drinksperweek = "Drinks Per Week",
  ea = "Educational Attainment",
  ebb = "Spirituality",
  epilepsy = "Epilepsy",
  extra = "Extraversion",
  E2 = "Verbal Tilt",
  famstat = "Satisfaction with Family",
  fetal = "Birth Weight",
  finsat = "Financial Satisfaction",
  freqfriendvisit = "Frequency of Family/Friend Visits",
  friendsat = "Satisfaction with Friends",
  guiltyfeelings = "Guilty Feelings",
  hdl = "HDL Cholesterol",
  height = "Height",
  income = "Income",
  inflammatoryboweldisease = "Inflammatory Bowel Disease",
  iq = "Cognitive Ability",
  jobsat = "Job Satisfaction",
  ldl = "LDL Cholesterol",
  loneliness = "Loneliness",
  longevity = "Longevity",
  lupus = "Lupus",
  matevalue = "Partner choice index",
  memory = "Memory",
  morningperson = "Morning Person",
  neb = "Number of Children",
  neuro = "Neuroticism",
  nopartners = "Number of Partners",
  occstatus = "Occupational Status",
  ocd = "Obsessive-Compulsive Disorder",
  open = "Openness",
  parkinson = "Parkinson's Disease",
  participationphs = "Phys. Health Survey Participation",
  participationprimary = "Biobank Participation",
  pp = "Pulse Pressure",
  ptsd = "Post-Traumatic Stress Disorder",
  religattend = "Religious Attendance",
  risktaking = "Risk Taking",
  rt = "Reaction Time",
  sbp = "Systolic Blood Pressure",
  schiz = "Schizophrenia",
  selection = "Ancient Selection",
  selfhealth = "Self-Rated Health",
  sleepduration = "Sleep Duration",
  smokecessation = "Smoking Cessation",
  smokinit = "Smoking Initiation",
  ssbf = "Same-Sex Behavior (Female)",
  ssbm = "Same-Sex Behavior (Male)",
  swb = "Subjective Well-Being",
  t2d = "Type 2 Diabetes",
  tiredness = "Tiredness",
  tourette = "Tourette Syndrome",
  townsend = "Townsend Deprivation Index",
  trig = "Triglycerides",
  ulcerativecolitis = "Ulcerative Colitis",
  vitamind = "Vitamin D",
  walkspeed = "Walking Speed",
  whr = "Waist-Hip Ratio",
  neurot = "Neuroticism"
)

#List trait files
trait_files <- list.files(
  path = group_dir,
  pattern = "\\.gsa\\.out$",
  full.names = TRUE
)

stopifnot(length(trait_files) > 0)

short_names <- basename(trait_files)
short_names <- sub("_.*$", "", short_names)

full_names <- ifelse(
  short_names %in% names(nice_names),
  nice_names[short_names],
  short_names
)

#Read and rename each trait file
traits_list <- lapply(seq_along(trait_files), function(k) {
  df <- read.table(trait_files[k],
                   header = TRUE,
                   comment.char = "#",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
  
  if (!"FULL_NAME" %in% colnames(df)) {
    df$FULL_NAME <- df$VARIABLE
  }
  
  df %>%
    dplyr::select(FULL_NAME, BETA_STD, SE) %>%
    dplyr::rename(
      !!paste0("BETA_STD_", short_names[k]) := BETA_STD,
      !!paste0("SE_", short_names[k])       := SE
    )
})

names(traits_list) <- short_names


#Merge all trait tables
merged <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "FULL_NAME"),
  traits_list
)
merged <- as.data.frame(merged)


#Build Beta and SE matrices
beta_cols <- grep("^BETA_STD_", colnames(merged), value = TRUE)
se_cols   <- grep("^SE_", colnames(merged), value = TRUE)

beta_mat <- merged[, c("FULL_NAME", beta_cols)]
se_mat   <- merged[, c("FULL_NAME", se_cols)]

colnames(beta_mat) <- c("FULL_NAME", short_names)
colnames(se_mat)   <- c("FULL_NAME", short_names)

valid_idx <- complete.cases(beta_mat, se_mat)
beta_mat <- beta_mat[valid_idx, ]
se_mat   <- se_mat[valid_idx, ]


# Unweighted Pearson correlations

n_traits <- length(short_names)

cor_mat <- matrix(
  NA,
  n_traits, n_traits,
  dimnames = list(full_names, full_names)
)

for (i in 1:n_traits) {
  for (j in i:n_traits) {
    bi <- beta_mat[[short_names[i]]]
    bj <- beta_mat[[short_names[j]]]
    
    ok <- is.finite(bi) & is.finite(bj)
    
    corr <- if (sum(ok) > 2)
      cor(bi[ok], bj[ok], method = "pearson")
    else NA
    
    cor_mat[i, j] <- corr
    cor_mat[j, i] <- corr
  }
}

cor_mat[!is.finite(cor_mat)] <- 0
if (min(cor_mat) == max(cor_mat))
  cor_mat[1, 1] <- cor_mat[1, 1] + 1e-6

#Save to pdf
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

CairoPDF(
  file = out_file,
  width = 30,
  height = 30
)

ph <- pheatmap(
  cor_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  color = colorRampPalette(c("red3", "white", "blue"))(100),
  breaks = seq(-1, 1, length.out = 101),
  main = paste0(group, " - Unweighted Pearson"),
  fontsize_row = 16,
  fontsize_col = 16,
  angle_col = 315,
  labels_col = paste0(" ", colnames(cor_mat)),
  legend_cex = 1.6,
  legend_width = 3,
  silent = TRUE
)

grid::grid.draw(ph$gtable)
dev.off()

cat("Heatmap saved:", out_file, "\n")