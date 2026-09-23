# ---------------------------------------------------------------------------
# Participation-bias-adjusted genetic correlations of intelligence with the
# Big Five, plotted against the phenotypic mean shift from population to sample.
#
# File naming: the number in the file name is INTELLIGENCE's mean shift, which
# is fixed at +0.44 in every file used here. The `mean_shift` column inside each
# file is the BIG FIVE TRAIT's mean shift, and is what the x axis varies.
# Neuroticism is the one trait whose own shift is negative, so its curve is
# marked at -0.44 while the other four are marked at +0.44.
#
# Expects the five tab-delimited input files (columns: trait1, trait2,
# mean_shift, adj_rg, se) and this script to live in the same folder.
#
# Requires ggplot2 >= 3.4 (for `linewidth`).
# ---------------------------------------------------------------------------

library(ggplot2)

## ---- settings -------------------------------------------------------------

DIR       <- "C:/Users/buona008/Downloads"
IQ_SHIFT  <- 0.44     # intelligence's mean shift; fixed, and encoded in the file names
X_RANGE   <- c(-1, 1) # range of Big Five mean shifts to plot
CI_MULT   <- 1.96     # ribbon half-width, in SEs. 1.96 = 95% CI; set to 1 for +/-1 SE
SAVE      <- TRUE     # write PNG + PDF to DIR

ZERO_LW   <- 0.7      # weight of the x = 0 and y = 0 reference lines
ZERO_COL  <- "grey25"
GUIDE_LW  <- 0.4      # weight of the dashed lines at each trait's own mean shift
GUIDE_COL <- "grey45"

# One row per trait. `mark` is the trait's own mean shift, highlighted on its curve.
spec <- data.frame(
  file   = c("extra_iq_adj_0.44.txt",
             "agree_iq_adj_0.44.txt",
             "consc_iq_adj_0.44.txt",
             "neurot_iq_adj_0.44.txt",
             "open_iq_adj_0.44.txt"),
  trait  = c("Extraversion",
             "Agreeableness",
             "Conscientiousness",
             "Neuroticism",
             "Openness/Intellect"),
  mark   = c(0.44, 0.44, 0.44, -0.44, 0.44),
  colour = c("#E69F00",   # extraversion      - orange
             "#E377A8",   # agreeableness     - pink
             "#8C2D19",   # conscientiousness - dark red
             "#D7301F",   # neuroticism       - red
             "#1B7837"),  # openness          - green
  stringsAsFactors = FALSE
)

## ---- read -----------------------------------------------------------------

read_one <- function(i) {
  path <- file.path(DIR, spec$file[i])
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  d <- read.delim(path, stringsAsFactors = FALSE)
  need <- c("mean_shift", "adj_rg", "se")
  missing <- setdiff(need, names(d))
  if (length(missing)) {
    stop("Missing column(s) ", paste(missing, collapse = ", "),
         " in ", spec$file[i], call. = FALSE)
  }
  d$trait <- spec$trait[i]
  d[, c("trait", "mean_shift", "adj_rg", "se")]
}

dat <- do.call(rbind, lapply(seq_len(nrow(spec)), read_one))
dat$trait <- factor(dat$trait, levels = spec$trait)

dat$lo <- dat$adj_rg - CI_MULT * dat$se
dat$hi <- dat$adj_rg + CI_MULT * dat$se

dat <- dat[dat$mean_shift >= X_RANGE[1] & dat$mean_shift <= X_RANGE[2], ]

## ---- the highlighted points ------------------------------------------------

marks <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i) {
  d <- dat[dat$trait == spec$trait[i], ]
  d[which.min(abs(d$mean_shift - spec$mark[i])), ]
}))

# per-trait guide lines (facet plot) and distinct guide lines (combined plot)
vlines   <- data.frame(trait = factor(spec$trait, levels = spec$trait),
                       mark  = spec$mark)
mark_lab <- data.frame(mark = sort(unique(spec$mark)))
mark_lab$label <- sprintf("%.2f", mark_lab$mark)

cat("\nAdjusted rg at each trait's own mean shift (intelligence fixed at ",
    IQ_SHIFT, "):\n\n", sep = "")
print(data.frame(trait      = marks$trait,
                 mean_shift = marks$mean_shift,
                 adj_rg     = round(marks$adj_rg, 3),
                 se         = round(marks$se, 3)),
      row.names = FALSE)
cat("\n")

pal <- setNames(spec$colour, spec$trait)

x_lab <- "Mean shift of the Big Five trait"
y_lab <- "Adjusted genetic correlation with intelligence"
cap   <- sprintf(paste0("Intelligence's mean shift fixed at %.2f. Bands: estimate +/- %.2f SE."),
                 IQ_SHIFT, CI_MULT)

## ---- plot 1: all five traits in one panel ---------------------------------

p_all <- ggplot(dat, aes(mean_shift, adj_rg, colour = trait, fill = trait)) +
  geom_hline(yintercept = 0, linewidth = ZERO_LW, colour = ZERO_COL) +
  geom_vline(xintercept = 0, linewidth = ZERO_LW, colour = ZERO_COL) +
  geom_vline(data = mark_lab, aes(xintercept = mark), linetype = "dashed",
             linewidth = GUIDE_LW, colour = GUIDE_COL) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(data = marks, size = 2.3) +
  geom_text(data = mark_lab, aes(x = mark, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.2, vjust = 1.6,
            size = 3.2, colour = GUIDE_COL) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = X_RANGE, breaks = seq(-1, 1, 0.25),
                     expand = c(0.01, 0)) +
  labs(x = x_lab, y = y_lab, caption = cap) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        legend.key.width = unit(1.4, "lines"),
        plot.caption = element_text(colour = "grey40", size = 8))

## ---- plot 2: one panel per trait ------------------------------------------

p_facet <- ggplot(dat, aes(mean_shift, adj_rg)) +
  geom_hline(yintercept = 0, linewidth = ZERO_LW, colour = ZERO_COL) +
  geom_vline(xintercept = 0, linewidth = ZERO_LW, colour = ZERO_COL) +
  geom_vline(data = vlines, aes(xintercept = mark), linetype = "dashed",
             linewidth = GUIDE_LW, colour = GUIDE_COL) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = trait), alpha = 0.25,
              colour = NA) +
  geom_line(aes(colour = trait), linewidth = 0.8) +
  geom_point(data = marks, aes(colour = trait), size = 2.3) +
  geom_text(data = marks,
            aes(label = sprintf("%.2f (%.2f)", adj_rg, se)),
            hjust = -0.15, vjust = -0.9, size = 3, colour = "grey20") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = X_RANGE, breaks = seq(-1, 1, 0.5),
                     expand = c(0.02, 0)) +
  facet_wrap(~ trait, nrow = 2) +
  labs(x = x_lab, y = y_lab, caption = cap) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey40", size = 8))

print(p_all)

## ---- save -----------------------------------------------------------------

if (SAVE) {
  ggsave(file.path(DIR, "pb_adjusted_rg_combined.png"), p_all,
         width = 6.5, height = 5.2, dpi = 400)
  ggsave(file.path(DIR, "pb_adjusted_rg_combined.pdf"), p_all,
         width = 6.5, height = 5.2)
  ggsave(file.path(DIR, "pb_adjusted_rg_faceted.png"), p_facet,
         width = 9, height = 5.6, dpi = 400)
  ggsave(file.path(DIR, "pb_adjusted_rg_faceted.pdf"), p_facet,
         width = 9, height = 5.6)
  write.csv(data.frame(trait      = marks$trait,
                       mean_shift = marks$mean_shift,
                       adj_rg     = marks$adj_rg,
                       se         = marks$se),
            file.path(DIR, "pb_adjusted_rg_at_mark.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Worth checking against Supplementary Table 10. From the four files I have
# read (all with intelligence fixed at 0.44), at trait shift = +0.44:
#
#       extraversion       -0.084 (0.030)      Table 10: -0.10 (0.04)
#       agreeableness       0.223 (0.032)      Table 10:  0.20 (0.04)
#       conscientiousness  -0.101 (0.034)      Table 10: -0.08 (0.04)
#       openness            0.393 (0.030)      Table 10:  0.38 (0.03)
#
# The extraversion and conscientiousness rows in Table 10 look like each
# other's values, and several SEs disagree by a rounding step. Regenerating
# the table from the CSV this script writes is safer than patching it.
# Neuroticism is not in that list because neurot_iq_adj_0.44.txt was not
# among the files I saw; check its value at -0.44 against Table 10's -0.24.
# ---------------------------------------------------------------------------