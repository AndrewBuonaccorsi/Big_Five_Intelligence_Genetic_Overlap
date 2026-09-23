
library(eulerr)
library(patchwork)
library(ggplotify)

df <- data.frame(
  Trait = c("Extraversion", "Agreeableness", "Conscientiousness",
            "Neuroticism", "Openness", "Height"),
  I_only = c(0.5, 1.1, 0.6, 0.2, 3.0, 9.0),
  Both   = c(10.0, 9.5, 10.0, 10.4, 7.6, 1.5),
  P_only = c(1.1, 0.8, 0.6, 3.1, 0.3, 3.0),
  I_SE   = c(0.4, 0.9, 0.6, 0.1, 1.2, 0.7),
  Both_SE= c(0.7, 1.2, 0.6, 0.7, 1.1, 0.3),
  P_SE   = c(0.9, 0.8, 0.6, 0.9, 0.9, 0.3)
)


intelligence_color <- "blue"
personality_colors <- c(
  Extraversion      = "darkorange",
  Agreeableness     = "hotpink",
  Conscientiousness = "brown4",
  Neuroticism       = "red",
  Openness          = "green4",
  Height            = "gray70"
)


make_venn <- function(row) {
  
  fit <- euler(c(
    Intelligence = row$I_only,
    Personality  = row$P_only,
    "Intelligence&Personality" = row$Both
  ))
  
  quantity_labels <- c(
    sprintf("%.1f (%.1f)", row$I_only, row$I_SE),
    sprintf("%.1f (%.1f)", row$P_only, row$P_SE),
    sprintf("%.1f (%.1f)", row$Both,   row$Both_SE)
  )
  
  p <- plot(
    fit,
    fills = list(
      fill = c(intelligence_color, personality_colors[row$Trait]),
      alpha = 0.6
    ),
    edges = TRUE,
    edge_color = "black",
    labels = FALSE,
    quantities = list(
      labels = quantity_labels,
      cex = 1.1,
      offset = 0.05
    ),
    main = row$Trait,
    main.pos = c(0.05, 0.5),  # left-center
    legend = FALSE,
    mar = c(4, 4, 2, 2)
  )
  
  as.ggplot(p)
}

plots <- lapply(df$Trait, function(trait) {
  row <- df[df$Trait == trait, ]
  make_venn(row)
})


top_row <- wrap_plots(plots[1:3], ncol = 3)
bottom_row <- wrap_plots(plots[4:6], ncol = 3)

final_plot <- top_row / bottom_row + plot_layout(heights = c(1,1))


final_plot
