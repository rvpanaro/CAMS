#
# Figure 4 left: Minimizing credible interval width — total subgroup CrI width vs
# information fraction; red dot at optimal prevalence. EDSS example.
# Self-contained. Requires: edss_data.csv. Reuses fit_edss.rds.
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

# ---- Posterior of fixed effects (population-level) ----
draws <- as.data.frame(fit_subgroup_adj)
nms   <- names(draws)
b_Intercept <- draws[[which(nms == "b_Intercept")]]
b_ifrac     <- draws[[which(nms == "b_ifrac")]]
sub_par     <- nms[grepl("subgroup", nms, fixed = TRUE) & startsWith(nms, "b_")][1]
b_subgroupB <- draws[[sub_par]]

# ---- Total credible interval width vs IF ----
# At each IF (prevalence): mean_A = Intercept + ifrac*b_ifrac, mean_B = Intercept + ifrac*b_ifrac + subgroupB
# 95% CrI width for A and B, then total width = width_A + width_B
ifrac_grid <- seq(0.1, 0.7, length.out = 150)
total_width <- numeric(length(ifrac_grid))
for (i in seq_along(ifrac_grid)) {
  pi_val <- ifrac_grid[i]
  mean_A <- b_Intercept + b_ifrac * pi_val
  mean_B <- b_Intercept + b_ifrac * pi_val + b_subgroupB
  ci_A <- quantile(mean_A, c(0.025, 0.975))
  ci_B <- quantile(mean_B, c(0.025, 0.975))
  total_width[i] <- diff(ci_A) + diff(ci_B)
}

# Optimal prevalence (IF that minimizes total width)
idx_opt <- which.min(total_width)
optimal_ifrac <- ifrac_grid[idx_opt]
optimal_width <- total_width[idx_opt]

df_curve <- data.frame(ifrac = ifrac_grid, width = total_width)
y_upper <- max(df_curve$width, na.rm = TRUE) * 1.08

# ---- Plot: total CrI width vs information fraction (y from 0, ribbon under curve) ----
p <- ggplot(df_curve, aes(x = ifrac, y = width)) +
  geom_ribbon(aes(ymin = 0, ymax = width), fill = "#4169E1", alpha = 0.18) +
  geom_line(colour = "#0000FF", linewidth = 1) +
  geom_point(data = data.frame(ifrac = optimal_ifrac, width = optimal_width), aes(x = ifrac, y = width),
             colour = "#FF0000", size = 4) +
  geom_segment(data = data.frame(x = optimal_ifrac, y = optimal_width),
               aes(x = x, xend = x, y = y, yend = y_upper),
               colour = "#FF0000", linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE) +
  labs(
    title = "Minimizing credible interval width",
    x = expression(Information~fraction~("Higher EDSS")),
    y = "Total credible interval width"
  ) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  coord_cartesian(xlim = c(0.1, 0.7), ylim = c(0, y_upper), clip = "on") +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

print(p)
ggsave("figures/figure4_lhs.pdf", plot = p, width = 7, height = 7, units = "in", device = grDevices::pdf)
