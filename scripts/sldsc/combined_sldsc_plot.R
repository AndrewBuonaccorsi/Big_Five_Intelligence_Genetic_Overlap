

library(tidyverse)
library(stringr)
library(patchwork)


trait_order <- c(
  "Intelligence",
  "Extraversion",
  "Agreeableness",
  "Conscientiousness",
  "Neuroticism",
  "Openness/intellect",
  "Mean"
)

legend_breaks <- c(
  "Intelligence",
  "Extraversion",
  "Agreeableness",
  "Conscientiousness",
  "Neuroticism",
  "Openness/intellect",
  "Mean"
)

trait_colors <- c(
  "Intelligence"      = "blue",
  "Extraversion"      = "darkorange",
  "Agreeableness"     = "hotpink",
  "Conscientiousness" = "brown4",
  "Neuroticism"       = "red",
  "Openness/intellect"          = "darkgreen",
  "Mean"              = "black"
)

pd <- position_dodge(width = 0.6)

master <- read_csv(
  "C:/Users/buona008/Downloads/sldsc_results/master_sldsc_results.csv"
)

brain_keywords <- "BRAIN|CORTEX|CEREBELL|HIPPOCAMP|AMYGDALA|HYPOTHALAM|CAUDATE|PUTAMEN|ACCUMBENS|SPINAL_CORD|SUBSTANTIA_NIGRA"

brain_df <- master %>%
  filter(str_starts(Tissue, "GTEX_")) %>%
  mutate(region_raw = str_remove(Tissue, "^GTEX_")) %>%
  filter(str_detect(region_raw, brain_keywords))

cat("Panel A:", nrow(brain_df), "GTEx brain-region rows matched.\n")

brain_df <- brain_df %>%
  mutate(
    region_clean = str_remove(region_raw, "^BRAIN_"),
    region_clean = str_replace_all(region_clean, "_", " "),
    region_clean = str_to_title(region_clean),
    
    BrainRegion = case_when(
      str_detect(region_raw, "FRONTAL_CORTEX")      ~ "Frontal Cortex",
      region_raw %in% c("RES_CORTEX")    ~ "Frontal Cortex (resample)",
      str_detect(region_raw, "ANTERIOR_CINGULATE")   ~ "Anterior cingulate cortex",
      str_detect(region_raw, "NUCLEUS_ACCUMBENS")    ~ "Nucleus accumbens",
      str_detect(region_raw, "CAUDATE")              ~ "Caudate nucleus",
      str_detect(region_raw, "PUTAMEN")              ~ "Putamen",
      str_detect(region_raw, "SPINAL_CORD")          ~ "Cervical spinal cord",
      str_detect(region_raw, "CEREBELLUM")~ "Cerebellum",
      region_raw %in% c("CEREBELLAR") ~ "Cerebellum (resample)",
      str_detect(region_raw, "HIPPOCAMPUS")          ~ "Hippocampus",
      str_detect(region_raw, "AMYGDALA")             ~ "Amygdala",
      str_detect(region_raw, "HYPOTHALAMUS")         ~ "Hypothalamus",
      str_detect(region_raw, "SUBSTANTIA_NIGRA")     ~ "Substantia nigra",
      TRUE ~ region_clean
    ),
    
    SE = SE
  )

mean_df_brain <- brain_df %>%
  filter(Trait != "Height") %>%
  group_by(BrainRegion) %>%
  summarise(
    mean_enrichment = mean(Enrichment, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_enrichment)) %>%
  mutate(Trait = "Mean")

region_order <- mean_df_brain$BrainRegion

brain_df$BrainRegion <- factor(brain_df$BrainRegion, levels = region_order)
mean_df_brain$BrainRegion <- factor(mean_df_brain$BrainRegion, levels = region_order)

brain_df$Trait <- factor(brain_df$Trait, levels = trait_order)
mean_df_brain$Trait <- factor(mean_df_brain$Trait, levels = trait_order)

brain_plot_df <- brain_df %>%
  filter(Trait != "Height")

p_brain <- ggplot(
  brain_plot_df,
  aes(x = BrainRegion, y = Enrichment, color = Trait)
) +
  
  geom_line(
    data = mean_df_brain,
    aes(x = BrainRegion, y = mean_enrichment, group = 1, color = Trait),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  geom_point(position = pd, size = 1.5, alpha = 0.9) +
  
  geom_errorbar(
    aes(ymin = Enrichment - 1.96 * SE, ymax = Enrichment + 1.96 * SE),
    position = pd,
    width = 0.15,
    alpha = 0.6,
    linewidth = 0.5
  ) +
  
  geom_hline(yintercept = 1, color = "black", linewidth = 0.6) +
  
  scale_color_manual(values = trait_colors, guide = guide_legend(override.aes = list(size = 3, alpha = 1, linewidth = 1)), breaks = legend_breaks) +
  
  labs(x = "Brain Region", y = "Enrichment") +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey60", linewidth = 0.7),
    panel.grid.minor.y = element_blank()
  )

dev_df <- master %>%
  filter(str_starts(Tissue, "BRAINSPAN_")) %>%
  mutate(
    age_raw   = str_remove(Tissue, "^BRAINSPAN_"),
    age_value = as.numeric(str_extract(age_raw, "[0-9]+")),
    
    age_unit = case_when(
      str_detect(age_raw, "PCW") ~ "pcw",
      str_detect(age_raw, "MOS") ~ "months",
      str_detect(age_raw, "YRS") ~ "years",
      TRUE ~ NA_character_
    ),
    
    age_weeks = case_when(
      age_unit == "pcw"    ~ age_value,
      age_unit == "months" ~ (age_value * 4.345) + 40,
      age_unit == "years"  ~ (age_value * 52.18) + 40,
      TRUE ~ NA_real_
    ),
    
    stage_label = paste(age_value, age_unit)
  ) %>%
  filter(!is.na(age_weeks))

cat("Panel B:", nrow(dev_df), "BrainSpan age rows matched.\n")

dev_df <- dev_df %>%
  arrange(age_weeks) %>%
  mutate(stage = factor(stage_label, levels = unique(stage_label)))

mean_df_dev <- dev_df %>%
  filter(Trait != "Height") %>%
  group_by(stage, age_weeks) %>%
  summarise(
    mean_enrichment = mean(Enrichment, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Trait = "Mean")

dev_df$Trait <- factor(dev_df$Trait, levels = trait_order)
mean_df_dev$Trait <- factor(mean_df_dev$Trait, levels = trait_order)

p_dev <- ggplot(
  dev_df,
  aes(x = stage, y = Enrichment, color = Trait)
) +
  
  geom_line(
    data = mean_df_dev,
    aes(x = stage, y = mean_enrichment, group = 1, color = Trait),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  
  geom_point(position = pd, size = 1.2, alpha = 0.9) +
  
  geom_errorbar(
    aes(ymin = Enrichment - 1.96 * SE, ymax = Enrichment + 1.96 * SE),
    position = pd,
    width = 0.2,
    alpha = 0.6,
    linewidth = 0.4
  ) +
  
  geom_hline(yintercept = 1, color = "black", linewidth = 0.6) +
  
  scale_color_manual(values = trait_colors, guide = guide_legend(override.aes = list(size = 3, alpha = 1, linewidth = 1)), breaks = legend_breaks) +
  
  labs(x = "Developmental Stage", y = "Enrichment") +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey60", linewidth = 0.7),
    panel.grid.minor.y = element_blank()
  )

p_dev <- p_dev + theme(legend.position = "none")

figure1 <-
  (p_brain / p_dev) +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  )

figure1 <- figure1 & theme(
  legend.position = "right",
  plot.tag = element_text(size = 20, face = "bold")
)

ggsave(
  filename = "Figure_BrainRegions_DevelopmentalStages.pdf",
  plot = figure1,
  width = 12,
  height = 10,
  units = "in"
)

