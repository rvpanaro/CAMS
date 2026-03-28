#
# Figure 4 right: Relapse rate (RR) by baseline EDSS — same as figure3_rhs (conditional
# effects with 95% CrI). Error bars show subgroup estimates at the optimal prevalence
# (IF that minimizes total credible interval width). Requires: edss_data.csv. Reuses fit_edss.rds.
#

if (!file.exists("edss_data.csv")) stop("edss_data.csv not found. Set wd to the project folder.")

set.seed(1)
library(brms)
library(ggplot2)
library(dplyr)

# ---- Data ----
long_data <- read.csv("edss_data.csv", stringsAsFactors = FALSE)
long_data$se <- sqrt(long_data$vi)
long_data$study <- paste0(long_data$study.acronym, long_data$treatment)
if (!"ifrac2" %in% names(long_data)) long_data$ifrac2 <- rep(c(0, 1), nrow(long_data) / 2) - long_data$ifrac

# ---- Fit (EDSS CAMS) ----
if (file.exists("fit_edss.rds")) {
  fit_subgroup_adj <- readRDS("fit_edss.rds")
} else {
  fit_subgroup_adj <- brms::brm(
    yi | se(se) ~ 1 + ifrac + subgroup + (1 + ifrac2 || study),
    data = long_data,
    prior = prior(normal(0, 1), class = sd, coef = "Intercept", group = "study") +
      prior(normal(0, 0.5), class = sd, coef = "ifrac2", group = "study"),
    chains = 4, iter = 20000, warmup = 6000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
  saveRDS(fit_subgroup_adj, "fit_edss.rds")
}

# ---- Optimal prevalence (same logic as figure4_lhs) ----
draws <- as.data.frame(fit_subgroup_adj)
nms   <- names(draws)
b_Intercept <- draws[[which(nms == "b_Intercept")]]
b_ifrac     <- draws[[which(nms == "b_ifrac")]]
sub_par     <- nms[grepl("subgroup", nms, fixed = TRUE) & startsWith(nms, "b_")][1]
b_subgroupB <- draws[[sub_par]]
ifrac_grid  <- seq(0.1, 0.7, length.out = 150)
total_width <- numeric(length(ifrac_grid))
for (i in seq_along(ifrac_grid)) {
  pi_val <- ifrac_grid[i]
  mean_A <- b_Intercept + b_ifrac * pi_val
  mean_B <- b_Intercept + b_ifrac * pi_val + b_subgroupB
  total_width[i] <- diff(quantile(mean_A, c(0.025, 0.975))) + diff(quantile(mean_B, c(0.025, 0.975)))
}
optimal_ifrac <- ifrac_grid[which.min(total_width)]

# ---- Estimates at optimal prevalence (for error bars) ----
mean_A_opt <- b_Intercept + b_ifrac * optimal_ifrac
mean_B_opt <- b_Intercept + b_ifrac * optimal_ifrac + b_subgroupB
ci_A <- quantile(mean_A_opt, c(0.025, 0.5, 0.975))
ci_B <- quantile(mean_B_opt, c(0.025, 0.5, 0.975))

# ---- Conditional effects (ifrac x subgroup) — same as figure3_rhs ----
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
lvls <- unique(as.character(ce_df$subgroup_f))
# Internal levels match figure3_rhs (brms uses A/B; labels map to threshold wording)
subgroup_labs <- c("<3.5", "ge35")
if (length(lvls) >= 2L) {
  ce_df$Subgroup <- factor(as.character(ce_df$subgroup_f), levels = lvls, labels = subgroup_labs[seq_along(lvls)])
} else {
  ce_df$Subgroup <- factor(ce_df$subgroup_f)
}
lv <- levels(ce_df$Subgroup)

opt_est <- data.frame(
  ifrac = optimal_ifrac,
  Subgroup = factor(subgroup_labs[seq_len(2)], levels = lv),
  estimate__ = c(ci_A[2], ci_B[2]),
  lower__ = c(ci_A[1], ci_B[1]),
  upper__ = c(ci_A[3], ci_B[3])
)

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

trial_segments <- raw %>%
  group_by(study) %>%
  summarise(ifrac = ifrac[1], yi_1 = yi[1], yi_2 = yi[2], .groups = "drop")

paper_clrs <- c("<3.5" = "#FF0000", ge35 = "#0000FF")
clrs <- setNames(paper_clrs[lv], lv)
if (length(clrs) == 0L) clrs <- c("#FF0000", "#0000FF")
edss_legend_labels <- lapply(lv, function(lvl) {
  switch(as.character(lvl), ge35 = expression(phantom() >= 3.5), `<3.5` = expression(phantom() < 3.5), lvl)
})

opt_dodge <- position_dodge(width = 0.025)

# ---- Plot: same as figure3_rhs + error bars at optimal prevalence ----
p <- ggplot(ce_df, aes(x = ifrac, y = estimate__, color = Subgroup, fill = Subgroup)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_errorbar(
    data = opt_est,
    aes(x = ifrac, ymin = lower__, ymax = upper__, color = Subgroup),
    linewidth = 0.8, width = 0.02, inherit.aes = FALSE, position = opt_dodge
  ) +
  geom_point(
    data = opt_est,
    aes(x = ifrac, y = estimate__, color = Subgroup, fill = Subgroup),
    size = 3.2, shape = 21, stroke = 0.9, inherit.aes = FALSE, position = opt_dodge
  ) +
  scale_color_manual(values = clrs, name = "EDSS", na.value = "grey75", labels = edss_legend_labels) +
  scale_fill_manual(values = clrs, name = "EDSS", na.value = "grey75", labels = edss_legend_labels) +
  scale_size_continuous(range = c(2.5, 6.5), guide = "none") +
  labs(
    x = expression(Information~fraction~("EDSS" >= 3.5)),
    y = "Relapse rate (RR) by baseline EDSS"
  ) +
  scale_x_continuous(breaks = c(0.2, 0.3, 0.4, 0.5, 0.6)) +
  scale_y_continuous(breaks = c(0, -0.5, -1.0, -1.5, -2.0, -2.5)) +
  coord_cartesian(xlim = c(0.1, 0.8), ylim = c(-2.5, 0.2)) +
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

print(p)
ggsave("figure4_rhs.pdf", plot = p, width = 7, height = 7, units = "in", device = grDevices::pdf)
