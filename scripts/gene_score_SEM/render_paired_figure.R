#!/usr/bin/env Rscript
required_packages <- c(
  "data.table",
  "ggplot2",
  "ggrepel",
  "ggthemes",
  "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) == 1L) {
  script_path <- normalizePath(
    sub("^--file=", "", script_argument),
    mustWork = TRUE
  )
} else {
  script_path <- normalizePath(
    "scripts/render_paired_figure.R",
    mustWork = TRUE
  )
}
project_directory <- dirname(dirname(script_path))
input_path <- file.path(
  project_directory,
  "output",
  "regression",
  "collection_point_uncertainty.csv"
)
figure_directory <- file.path(project_directory, "output", "figures")
if (!file.exists(input_path)) {
  stop(
    "Precomputed plot data not found: ",
    input_path,
    "\nRun scripts/analysis.R once to create it.",
    call. = FALSE
  )
}
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
theme_publication <- function(base_size = 17, base_family = "sans") {
  ggthemes::theme_foundation(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold"),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(
        fill = scales::alpha("white", 0.88),
        colour = NA
      ),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.position = "inside",
      legend.position.inside = c(0.03, 0.97),
      legend.justification.inside = c(0, 1),
      legend.direction = "vertical",
      legend.title.position = "top",
      legend.margin = margin(5, 7, 5, 7)
    )
}
collection_label_map <- c(
  GO_cell_component = "GO cellular component",
  akingbuwa_gs = "Akingbuwa",
  brainspan = "BrainSpan",
  computational_perturbation_signatures = "Computational perturbations",
  gtex = "GTEx tissues",
  hallmark_processes = "Hallmark processes",
  kegg_medicus_pathways = "KEGG Medicus",
  microRNA_targets = "microRNA targets"
)
collection_plot_data <- fread(input_path)
two_series_labels <- c(
  iq_big5 = "IQ - Big Five",
  big5_big5 = "Big Five - Big Five"
)
two_series_plot_data <- rbindlist(list(
  collection_plot_data[
    ,
    .(
      avg_reliability,
      correlation = avg_iq_correlation,
      series_label = unname(two_series_labels[["iq_big5"]])
    )
  ],
  collection_plot_data[
    ,
    .(
      avg_reliability,
      correlation = avg_big5_correlation,
      series_label = unname(two_series_labels[["big5_big5"]])
    )
  ]
))
two_series_plot_data[
  ,
  series_label := factor(
    series_label,
    levels = c("Big Five - Big Five", "IQ - Big Five")
  )
]
two_series_colors <- c(
  "Big Five - Big Five" = "#2C7FB8",
  "IQ - Big Five" = "#D95F02"
)
midpoint_obstacles <- collection_plot_data[
  ,
  .(
    collection,
    avg_reliability,
    paired_correlation_midpoint = (
      avg_iq_correlation + avg_big5_correlation
    ) / 2,
    annotation_label = fifelse(
      collection %chin% names(collection_label_map),
      unname(collection_label_map[collection]),
      ""
    )
  )
]
endpoint_obstacles <- rbindlist(list(
  collection_plot_data[
    ,
    .(
      collection,
      avg_reliability,
      paired_correlation_midpoint = avg_iq_correlation,
      annotation_label = ""
    )
  ],
  collection_plot_data[
    ,
    .(
      collection,
      avg_reliability,
      paired_correlation_midpoint = avg_big5_correlation,
      annotation_label = ""
    )
  ]
))
repel_data <- rbindlist(list(midpoint_obstacles, endpoint_obstacles))
paired_plot <- ggplot() +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "#A3A3A3",
    linetype = 2,
    linewidth = 0.65
  ) +
  geom_segment(
    data = collection_plot_data,
    aes(
      x = avg_reliability,
      xend = avg_reliability,
      y = avg_iq_correlation,
      yend = avg_big5_correlation
    ),
    color = "#A3A3A3",
    linewidth = 0.65
  ) +
  geom_point(
    data = two_series_plot_data,
    aes(
      x = avg_reliability,
      y = correlation,
      color = series_label
    ),
    size = 4.4
  ) +
  ggrepel::geom_label_repel(
    data = repel_data,
    aes(
      x = avg_reliability,
      y = paired_correlation_midpoint,
      label = annotation_label
    ),
    size = 4.4,
    seed = 20260729,
    box.padding = 0.7,
    label.padding = grid::unit(0.18, "lines"),
    point.padding = 0,
    label.r = grid::unit(0.2, "lines"),
    linewidth = 0.22,
    min.segment.length = 0,
    force = 4,
    force_pull = 0.25,
    max.time = 8,
    max.iter = 100000,
    max.overlaps = Inf,
    direction = "both",
    fill = scales::alpha("#EAF7FC", 0.92),
    color = "#30414D",
    segment.color = scales::alpha("#71838F", 0.82),
    segment.size = 0.45,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = NULL,
    values = two_series_colors,
    drop = FALSE
  ) +
  labs(
    x = "Average same-trait gene set correlation",
    y = "Average cross-trait gene set correlation"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.06))) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.06))) +
  coord_fixed(ratio = 1) +
  theme_publication()
png_path <- file.path(
  figure_directory,
  "geneset_reliability.png"
)
pdf_path <- file.path(
  figure_directory,
  "geneset_reliability.pdf"
)
ggsave(
  png_path,
  paired_plot,
  width = 8.2,
  height = 7.4,
  dpi = 240
)
ggsave(
  pdf_path,
  paired_plot,
  width = 8.2,
  height = 7.4
)
message("Wrote retained paired-series figure to ", png_path, " and ", pdf_path)
