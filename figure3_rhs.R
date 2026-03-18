#
# Figure 1 right panel (conditional effects: RR by EDSS vs information fraction).
# Same as figure1_lhs but for EDSS data. Requires: edss_data.csv. Displays in RStudio Plots tab.
#

if (!file.exists("edss_data.csv")) stop("edss_data.csv not found. Set wd to the project folder.")

set.seed(1)
library(brms)
library(ggplot2)
library(dplyr)

# Read data directly from edss_data.csv (columns: endpoint, study.acronym, treatment, subgroup, n, yi, vi, ifrac, ifrac2)
long_data <- read.csv("edss_data.csv", stringsAsFactors = FALSE)
if (!"ifrac2" %in% names(long_data)) {
  long_data$ifrac2 <- rep(c(0, 1), nrow(long_data) / 2) - long_data$ifrac
}
# Only add se and study (required by brms); use yi, vi, etc. from the data as-is
long_data$se <- sqrt(long_data$vi)
long_data$study <- paste0(long_data$study.acronym, long_data$treatment)

# Use cached fit if present (saved in project folder); otherwise fit and cache
cache_fit <- "fit_edss.rds"
if (file.exists(cache_fit)) {
  fit_subgroup_adj <- readRDS(cache_fit)
} else {
  fit_subgroup_adj <- brms::brm(
    yi | se(se) ~ 1 + ifrac + subgroup + (1 + ifrac2 || study),
    data = long_data,
    prior = prior(normal(0, 1), class = sd, coef = "Intercept", group = "study") +
      prior(normal(0, 0.5), class = sd, coef = "ifrac2", group = "study"),
    chains = 4, iter = 4000, warmup = 2000, seed = 1,
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
subgroup_labs <- c("Lower EDSS", "Higher EDSS")
lvls <- unique(as.character(ce_df$subgroup_f))
if (length(lvls) >= 2L) {
  ce_df$Subgroup <- factor(as.character(ce_df$subgroup_f), levels = lvls, labels = subgroup_labs[seq_along(lvls)])
} else {
  ce_df$Subgroup <- factor(ce_df$subgroup_f)
}

overall_df <- ce_df %>%
  group_by(.data$ifrac) %>%
  summarise(estimate__ = mean(.data$estimate__), .groups = "drop")

raw <- as.data.frame(fit_subgroup_adj$data)
raw$subgroup_f <- raw$subgroup
if (length(lvls) >= 2L) {
  raw$Subgroup <- factor(as.character(raw$subgroup), levels = lvls, labels = subgroup_labs[seq_along(lvls)])
} else {
  raw$Subgroup <- factor(raw$subgroup)
}
raw$precision <- 1 / pmax(raw$se, 1e-6)

# One segment per trial: gray line connecting the two yi values (same trial, subgroups A and B)
trial_segments <- raw %>%
  group_by(study) %>%
  summarise(ifrac = ifrac[1], yi_1 = yi[1], yi_2 = yi[2], .groups = "drop")

# Paper colours (same as LHS): red and blue
paper_clrs <- c("Lower EDSS" = "#FF0000", "Higher EDSS" = "#0000FF")
lv <- levels(ce_df$Subgroup)
clrs <- setNames(paper_clrs[lv], lv)
if (length(clrs) == 0L) clrs <- c("#FF0000", "#0000FF")

p_out <- ggplot(ce_df, aes(x = ifrac, y = estimate__, color = Subgroup, fill = Subgroup)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_line(data = overall_df, aes(x = ifrac, y = estimate__), colour = "black", linewidth = 1.2, linetype = "dashed", inherit.aes = FALSE) +
  geom_segment(data = trial_segments, aes(x = ifrac, xend = ifrac, y = yi_1, yend = yi_2), colour = "grey75", linewidth = 1, inherit.aes = FALSE) +
  geom_point(data = raw, aes(x = ifrac, y = yi, color = Subgroup, size = precision), alpha = 0.35, inherit.aes = FALSE) +
  scale_color_manual(values = clrs, name = "EDSS", na.value = "grey75") +
  scale_fill_manual(values = clrs, name = "EDSS", na.value = "grey75") +
  scale_size_continuous(range = c(2.5, 6.5), guide = "none") +
  labs(x = "Information fraction", y = "Relapse rate (RR) by baseline EDSS") +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6)) +
  scale_y_continuous(breaks = c(0, -0.5, -1.0, -1.5, -2.0, -2.5)) +
  coord_cartesian(xlim = c(0.1, 0.8), ylim = c(-2.5, 0.2)) +
  # Theme and size: same as figure1_lhs
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
# To save at paper size: uncomment:
# ggsave("figure1_rhs.pdf", plot = p_out, width = 7, height = 7, units = "in")
