library(ggplot2)
library(dplyr)
library(ggrepel)
library(stringr)

df <- read.table(
  "delta.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

nice_names <- c(
  adhd                  = "ADHD",
  afb                   = "Age at First Birth",
  afi                   = "Age at First Intercourse",  
  agemenarche           = "Age at Menarche",
  agemenopause          = "Age at Menopause",
  agesmokinit           = "Age Smoking Initiation",
  agree                 = "Agreeableness",
  als                   = "Amyotrophic Lateral Sclerosis",
  alzheimer             = "Alzheimer's Disease",
  anorexia              = "Anorexia Nervosa",
  anxiety               = "Anxiety",
  arthritis             = "Rheumatoid Arthritis",
  asd                   = "Autism Spectrum Disorder",
  asthma                = "Asthma",
  atopicdermatitis      = "Atopic Dermatitis",
  bipolar               = "Bipolar Disorder",
  bmi                   = "Body Mass Index",
  bodyfat               = "Body Fat Percentage",
  cad                   = "Coronary Artery Disease",
  childless             = "Childlessness",
  cigsperday            = "Cigarettes Per Day",
  coffeeperday          = "Coffee Per Day",
  covid                 = "Hospitalized with Covid",
  consc                 = "Conscientiousness",
  crohns                = "Crohn's Disease",
  darkskin              = "Dark Skin",
  dbp                   = "Diastolic Blood Pressure",
  depression            = "Depression",
  drinksperweek         = "Drinks Per Week",
  ea                    = "Educational Attainment",
  'private_sumstats/ebb_spiritual_person.sumstats.gz'  = "Spirituality",
  ebb                   = "Spirituality",
  epilepsy              = "Epilepsy",
  extra                 = "Extraversion",
  E2                    = "Verbal Tilt",
  famstat               = "Satisfaction with Family",
  fetal                 = "Birth Weight",  
  finsat                = "Financial Satisfaction",
  freqfriendvisit       = "Frequency of Family/Friend Visits",
  friendsat             = "Satisfaction with Friends",
  guiltyfeelings        = "Guilty Feelings",
  hdl                   = "HDL Cholesterol",
  height                = "Height",
  income                = "Income",
  inflammatoryboweldisease = "Inflammatory Bowel Disease",
  iq                    = "Cognitive Performance",
  jobsat                = "Job Satisfaction",
  ldl                   = "LDL Cholesterol",
  loneliness            = "Loneliness",
  longevity             = "Longevity",
  lupus                 = "Lupus",
  matevalue             = "Partner choice index",
  memory                = "Memory",
  morningperson         = "Morning Person",
  neb                   = "Number of Children",  
  neuro                 = "Neuroticism",
  nopartners            = "Number of Partners",
  occstatus             = "Occupational Status",
  ocd                   = "Obsessive-Compulsive Disorder",
  open                  = "Openness",
  parkinson             = "Parkinson's Disease",
  pp                    = "Pulse Pressure",
  ptsd                  = "Post-Traumatic Stress Disorder",
  religattend           = "Religious Attendance",
  risktaking            = "Risk Taking",
  rt                    = "Reaction Time",
  sbp                   = "Systolic Blood Pressure",
  schiz                 = "Schizophrenia",
  selection             = "Ancient Selection",  
  selfhealth            = "Self-Rated Health",
  sleepduration         = "Sleep Duration",
  smokecessation        = "Smoking Cessation",
  smokinit              = "Smoking Initiation",
  ssbf                  = "Same-Sex Behavior (Female)", 
  ssbm                  = "Same-Sex Behavior (Male)",  
  swb                   = "Subjective Well-Being",
  t2d                   = "Type 2 Diabetes",
  tiredness             = "Tiredness",
  tourette              = "Tourette Syndrome",
  townsend              = "Townsend Deprivation Index",
  trig                  = "Triglycerides",
  ulcerativecolitis     = "Ulcerative Colitis",
  vitamind              = "Vitamin D",
  walkspeed             = "Walking Speed",
  whr                   = "Waist-Hip Ratio",
  open                  = "Openness",
  consc                 = "Conscientiousness",
  extra                 = "Extraversion",
  agree                 = "Agreeableness",
  neurot                = "Neuroticism"
)

reverse_code <- c(
  "finsat", 
  "famstat",
  "friendsat",
  "jobsat",
  "morningperson",
  "selfhealth",
  "freqfriendvisit",
  "childless"
)


df <- df %>%
  mutate(
    q_delta = p.adjust(p_delta, method = "fdr")
  )

df <- df %>%
  mutate(
    phenotype_key = str_extract(PHENOTYPE, "^[a-zA-Z0-9]+"),
    nice_pheno = nice_names[phenotype_key]
  )

df <- df %>%
  mutate(
    reverse = phenotype_key %in% reverse_code,
    DIRECT   = ifelse(reverse, -DIRECT, DIRECT),
    MARGINAL = ifelse(reverse, -MARGINAL, MARGINAL),
    DELTA    = ifelse(reverse, -DELTA, DELTA)
  )

df <- df %>%
  mutate(
    Trait = case_when(
      grepl("cog|iq|intell", PREDICTOR, ignore.case = TRUE) ~ "Cognitive ability",
      grepl("agree", PREDICTOR, ignore.case = TRUE) ~ "Agreeableness",
      grepl("consc", PREDICTOR, ignore.case = TRUE) ~ "Conscientiousness",
      grepl("neuro", PREDICTOR, ignore.case = TRUE) ~ "Neuroticism",
      grepl("open", PREDICTOR, ignore.case = TRUE) ~ "Openness",
      grepl("extra", PREDICTOR, ignore.case = TRUE) ~ "Extraversion",
      TRUE ~ NA_character_
    )
  )

trait_colors <- c(
  "Agreeableness"       = "darkmagenta",
  "Conscientiousness"   = "brown",
  "Neuroticism"         = "red",
  "Openness"            = "darkgreen",
  "Extraversion"        = "darkorange"
)


ggplot(df, aes(x = MARGINAL, y = DIRECT, color = Trait)) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted",
    color = "black"
  ) +
  
  geom_vline(xintercept = 0, color = "black", linewidth = 1.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1.5) +
  
  geom_point(size = 1.5, alpha = 0.8) +
  
  geom_text_repel(
    data = subset(df, q_delta < 0.01 & abs(DELTA) > 0.1),
    aes(label = nice_pheno),
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.size = 0.5,
    force = 2,
    min.segment.length = 0
  ) +

  
  scale_color_manual(values = trait_colors) +
  
  labs(
    x = "Marginal Effect",
    y = "Direct Effect",
    color = "Predictor"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black")
  ) +
coord_cartesian(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6))



