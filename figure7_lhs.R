#
# Figure 7 LHS: Forest plot of prevalence estimates by study, with pooled and
# predicted prevalence (from RBesT gMAP). Self-contained; no CSV required.
# Requires: RBesT, interaction_forest_plot.R in same folder.
#

library("RBesT")
library(ggplot2)
library(dplyr)

# Study labels and counts (B.n = subgroup B, A.n = subgroup A; matches paper Figure 7)
labels <- c("DEFINE", "FREEDOMS", "CONFIRM-B", "CONFIRM-C", "ADVANCE", "AFFIRM", "TEMSO", "COP1 PILOT")
B.n <- c(73, 62, 60, 60, 89, 79, 164, 12)
A.n <- c(336, 363, 299, 290, 423, 548, 557, 13)

set.seed(1)
map01 <- gMAP(cbind(B.n, A.n) ~ 1 | labels,
              family = binomial,
              tau.dist = "HalfNormal",
              tau.prior = 1.0,
              iter = 2000 + 10000,
              beta.prior = cbind(0.0, 1.814))

mplots <- plot(map01)
fp <- mplots[[3]]

# Same filter as viz_map: keep study rows (stratified) and MAP/Mean rows
res <- fp@data |>
  dplyr::filter(study %in% c("MAP", "Mean") + model %in% "stratified" > 0)

# subgroup_intervals: Pooled (row 10) and Predicted (row 9); interaction_forest_plot expects yi, ci.lb, ci.ub
subgroup_intervals <- as.data.frame(res[c(10, 9), c(1, 4, 5)])
colnames(subgroup_intervals) <- c("yi", "ci.lb", "ci.ub")
rownames(subgroup_intervals) <- c("Pooled prevalence\n", "Predicted prevalence\n")

source("interaction_forest_plot.R")

prob <- 0.95
# LHS styling: individual studies = black unfilled square (shape 22); Pooled/Predicted = solid blue diamond
fp1 <- interaction_forest_plot(
  yi = res$mean[1:8],
  vi = res$sem[1:8]^2,
  labels = labels,
  exponentiate = FALSE,
  subgroup_intervals = subgroup_intervals,
  subgroups = c("", ""),
  pointsize = 5,
  interaction_intervals = rbind(c(NA, NA, NA), c(NA, NA, NA)),
  xlim = c(0.07, 0.8),
  xticks = c(0.25, 0.5, 0.75),
  col_right_heading = list(paste0("Prevalence \n(", 100 * prob, "% interval)"), ""),
  model_names = c(""),
  study_fill = "white",
  study_colour = "black",
  summary_fill = "#0000FF",
  summary_colour = "#0000FF",
  which_panel = 1L
)

# which_panel = 1 built only panel 1 (prevalence forest plot); use fp1$plot directly
fp2 <- fp1$plot +
  theme(plot.margin = margin(t = 40, r = 150, b = 8, l = 40))

print(fp2)
# ggsave("figure7_lhs.pdf", plot = fp2, width = 10, height = 7, units = "in")
