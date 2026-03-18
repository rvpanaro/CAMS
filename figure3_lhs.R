#
# Figure 1 left panel (conditional effects: RR by Age vs information fraction).
# Requires: age_data.csv. Displays in RStudio Plots tab.
#

if (!file.exists("age_data.csv")) stop("age_data.csv not found. Set wd to the project folder.")

set.seed(1)
library(brms)
library(ggplot2)
library(dplyr)

# Read data directly from age_data.csv (columns: endpoint, study.acronym, treatment, subgroup, n, yi, vi, ifrac, ifrac2)
long_data <- read.csv("age_data.csv", stringsAsFactors = FALSE)
if (!"ifrac2" %in% names(long_data)) {
  long_data$ifrac2 <- rep(c(0, 1), nrow(long_data) / 2) - long_data$ifrac
}
# Only add se and study (required by brms); use yi, vi, etc. from the data as-is
long_data$se <- sqrt(long_data$vi)
long_data$study <- paste0(long_data$study.acronym, long_data$treatment)

# Use cached fit if present (saved in project folder); otherwise fit and cache
cache_fit <- "fit_age.rds"
if (file.exists(cache_fit)) {
  fit_subgroup_adj <- readRDS(cache_fit)
} else {
  fit_subgroup_adj <- brms::brm(
    yi | se(se) ~ 1 + ifrac + subgroup + (1 + ifrac2 || study),
    data = long_data,
    prior = prior(normal(0, 1), class = sd, coef = "Intercept", group = "study") +
      prior(normal(0, 0.5), class = sd, coef = "ifrac2", group = "study"),
    chains = 4, iter = 12000, warmup = 4000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
  saveRDS(fit_subgroup_adj, cache_fit)
}

# Conditional effects (ifrac x subgroup)
ce <- conditional_effects(
  fit_subgroup_adj,
  effects = "ifrac:subgroup",
  re_formula = NA,
  prob = 0.95,
  resolution = 200
)
ce_df <- as.data.frame(ce[[1L]])
sub_col <- setdiff(names(ce_df), c("ifrac", "estimate__", "se__", "lower__", "upper__", "cond__"))
sub_col <- sub_col[grepl("subgroup", sub_col, fixed = TRUE)][1L]
if (is.na(sub_col)) sub_col <- names(ce_df)[2L]
ce_df$subgroup_f <- ce_df[[sub_col]]
age_labs <- c("<40", "≥40")
lvls <- unique(as.character(ce_df$subgroup_f))
if (length(lvls) >= 2L) {
  ce_df$Age <- factor(as.character(ce_df$subgroup_f), levels = lvls, labels = age_labs[seq_along(lvls)])
} else {
  ce_df$Age <- factor(ce_df$subgroup_f)
}

overall_df <- ce_df %>%
  group_by(.data$ifrac) %>%
  summarise(estimate__ = mean(.data$estimate__), .groups = "drop")

raw <- as.data.frame(fit_subgroup_adj$data)
raw$subgroup_f <- raw$subgroup
if (length(lvls) >= 2L) {
  raw$Age <- factor(as.character(raw$subgroup), levels = lvls, labels = age_labs[seq_along(lvls)])
} else {
  raw$Age <- factor(raw$subgroup)
}
raw$precision <- 1 / pmax(raw$se, 1e-6)

# One segment per trial: gray line connecting the two yi values (same trial, subgroups A and B)
trial_segments <- raw %>%
  group_by(study) %>%
  summarise(ifrac = ifrac[1], yi_1 = yi[1], yi_2 = yi[2], .groups = "drop")

# Paper colours (classic red and blue): <40 = red, ≥40 = blue
paper_clrs <- c("<40" = "#FF0000", "≥40" = "#0000FF")
lv <- levels(ce_df$Age)
clrs <- setNames(paper_clrs[lv], lv)
if (length(clrs) == 0L) clrs <- c("#FF0000", "#0000FF")

p_out <- ggplot(ce_df, aes(x = ifrac, y = estimate__, color = Age, fill = Age)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_line(data = overall_df, aes(x = ifrac, y = estimate__), colour = "black", linewidth = 1.2, linetype = "dashed", inherit.aes = FALSE) +
  geom_segment(data = trial_segments, aes(x = ifrac, xend = ifrac, y = yi_1, yend = yi_2), colour = "grey75", linewidth = 1, inherit.aes = FALSE) +
  geom_point(data = raw, aes(x = ifrac, y = yi, color = Age, size = precision), alpha = 0.35, inherit.aes = FALSE) +
  scale_color_manual(values = clrs, name = "Age", na.value = "grey75") +
  scale_fill_manual(values = clrs, name = "Age", na.value = "grey75") +
  scale_size_continuous(range = c(2.5, 6.5), guide = "none") +
  labs(x = "Information fraction", y = "Relapse rate (RR) by Age") +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6)) +
  scale_y_continuous(breaks = c(0, -0.5, -1.0, -1.5, -2.0, -2.5)) +
  coord_cartesian(xlim = c(0.25, 0.55), ylim = c(-1.6, 0.2)) +
  # Theme and size from aggregation_bias_plot1 (viz_theme.R + viz_condeffect.R): 7 x 7 in
  theme_classic(base_size = 20, base_family = "sans") +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.box.background = element_blank(),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )

print(p_out)
# To save at paper size (same as img/aggregation_bias_plot1): uncomment:
# ggsave("figure1_lhs.pdf", plot = p_out, width = 7, height = 7, units = "in")
