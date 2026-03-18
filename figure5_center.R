#
# Figure 5 center: Ecological effect / aggregation bias assessment.
# Self-contained. Requires: edss_data.csv. Reuses fit_overall_edss, fit_subgroup_edss, fit_interaction_edss, fit_edss.
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
long_data$study.name <- long_data$study
long_data$esti <- long_data$yi
long_data$subgroup12 <- as.numeric(factor(long_data$subgroup)) - 1.5
if (!"ifrac2" %in% names(long_data)) long_data$ifrac2 <- rep(c(0, 1), nrow(long_data) / 2) - long_data$ifrac

idx_a <- seq(1L, nrow(long_data), by = 2L)
idx_b <- seq(2L, nrow(long_data), by = 2L)
w_va <- long_data$vi[idx_a]
w_vb <- long_data$vi[idx_b]
w_ya <- long_data$yi[idx_a]
w_yb <- long_data$yi[idx_b]
wide_data <- data.frame(
  study.name = long_data$study.name[idx_a],
  contrast.esti = w_yb - w_ya,
  contrast.se = sqrt(w_va + w_vb),
  overall.esti = (w_ya / w_va + w_yb / w_vb) / (1 / w_va + 1 / w_vb),
  overall.se = sqrt(1 / (1 / w_va + 1 / w_vb))
)

# ---- Fits (EDSS) ----
if (file.exists("fit_subgroup_edss.rds")) fit_subgroup <- readRDS("fit_subgroup_edss.rds") else {
  fit_subgroup <- brms::brm(esti | se(se) ~ 1 + subgroup + (-1 + subgroup12 | study.name), data = long_data,
    prior = prior(normal(0, 0.5), class = sd), chains = 4, iter = 4000, warmup = 2000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)); saveRDS(fit_subgroup, "fit_subgroup_edss.rds")
}
if (file.exists("fit_interaction_edss.rds")) fit_interaction <- readRDS("fit_interaction_edss.rds") else {
  fit_interaction <- brms::brm(contrast.esti | se(contrast.se) ~ 1 + (1 | study.name), data = wide_data,
    prior = prior(normal(0, 0.5), class = sd), chains = 4, iter = 4000, warmup = 2000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)); saveRDS(fit_interaction, "fit_interaction_edss.rds")
}
if (file.exists("fit_edss.rds")) fit_subgroup_adj <- readRDS("fit_edss.rds") else {
  fit_subgroup_adj <- brms::brm(yi | se(se) ~ 1 + ifrac + subgroup + (1 + ifrac2 || study), data = long_data,
    prior = prior(normal(0, 1), class = sd, coef = "Intercept", group = "study") +
      prior(normal(0, 0.5), class = sd, coef = "ifrac2", group = "study"),
    chains = 4, iter = 4000, warmup = 2000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15))
  saveRDS(fit_subgroup_adj, "fit_edss.rds")
}
if (file.exists("fit_overall_edss.rds")) fit_overall <- readRDS("fit_overall_edss.rds") else {
  fit_overall <- brms::brm(overall.esti | se(overall.se) ~ 1 + (1 | study.name), data = wide_data,
    prior = prior(normal(0, 1), class = sd), chains = 4, iter = 4000, warmup = 2000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)); saveRDS(fit_overall, "fit_overall_edss.rds")
}

# ---- overall_ifrac and delta (project fit uses group "study") ----
d <- fit_subgroup_adj$data
re_nm <- names(summary(fit_subgroup_adj)$random)[1L]
if (is.null(re_nm)) re_nm <- "study"
tau <- brms::hypothesis(fit_subgroup_adj, "study__Intercept = 0", class = "sd")$samples |> unlist()
grp <- (seq_len(nrow(d)) + 1) %/% 2
pf <- d$ifrac - rep(c(0, 1), times = nrow(d) / 2)
den_se <- tapply(d$se^2 * pf^2, grp, sum)
wi <- sapply(tau, function(t) 1 / (t^2 + den_se))
pi <- tapply(d$ifrac, grp, `[`, 1)
overall_ifrac <- apply(wi, 2, FUN = function(w) crossprod(pi, w) / sum(w))

delta <- brms::hypothesis(fit_subgroup_adj, paste0(mean(overall_ifrac), " *(ifrac + subgroupB)< 0"))$samples |> unlist()
density_delta <- density(delta, adjust = 3)
df_delta <- data.frame(x = density_delta$x, y = density_delta$y)
prob_positive <- round(100 * mean(delta < 0))

# ---- Plot: Ecological effect / aggregation bias ----
ylim_range <- range(df_delta$y)
p1 <- ggplot(df_delta, aes(x = x, y = y)) +
  geom_line(linewidth = 1) +
  geom_area(data = df_delta |> filter(x < 0), aes(x = x, y = y), fill = "blue", alpha = 0.3) +
  labs(
    title = "Ecological effect assessment",
    x = expression("Overall prevalence x ecological slope" ~ (pi[overall] ~ (delta + gamma))),
    y = "Posterior density"
  ) +
  xlim(-2, 1) +
  annotate("text", x = 0.57, y = max(ylim_range) * 0.9, size = 5,
    label = sprintf("P(pi[overall] * (delta + gamma) < 0) == %d*'%%'", prob_positive), parse = TRUE) +
  coord_cartesian(ylim = ylim_range) +
  theme_classic(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    text = element_text(size = 22),
    panel.grid = element_blank()
  )

print(p1)
ggsave("figure5_center.pdf", plot = p1, width = 8, height = 7, units = "in")
