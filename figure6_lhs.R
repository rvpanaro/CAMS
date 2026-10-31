#
# Figure 6 left: Posterior interaction inference comparison (Interaction-MA, Subgroup-MA, CAMS).
# Self-contained. Requires: edss_data.csv. Reuses fit_subgroup_edss, fit_interaction_edss, fit_edss.
#

if (!file.exists("edss_data.csv")) stop("edss_data.csv not found. Set wd to the project folder.")

set.seed(1)
library(brms)
library(ggplot2)
library(dplyr)
library(gghalves)

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
wide_data <- data.frame(
  study.name = long_data$study.name[idx_a],
  contrast.esti = long_data$yi[idx_b] - long_data$yi[idx_a],
  contrast.se = sqrt(w_va + w_vb)
)

# ---- Fits (EDSS) ----
if (file.exists("fit_subgroup_edss.rds")) fit_subgroup <- readRDS("fit_subgroup_edss.rds") else {
  fit_subgroup <- brms::brm(esti | se(se) ~ 1 + subgroup + (-1 + subgroup12 | study.name), data = long_data,
    prior = prior(normal(0, 0.5), class = sd), chains = 4, iter = 20000, warmup = 6000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)); saveRDS(fit_subgroup, "fit_subgroup_edss.rds")
}
if (file.exists("fit_interaction_edss.rds")) fit_interaction <- readRDS("fit_interaction_edss.rds") else {
  fit_interaction <- brms::brm(contrast.esti | se(contrast.se) ~ 1 + (1 | study.name), data = wide_data,
    prior = prior(normal(0, 0.5), class = sd), chains = 4, iter = 20000, warmup = 6000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)); saveRDS(fit_interaction, "fit_interaction_edss.rds")
}
if (file.exists("fit_edss.rds")) fit_subgroup_adj <- readRDS("fit_edss.rds") else {
  fit_subgroup_adj <- brms::brm(yi | se(se) ~ 1 + ifrac + subgroup + (1 + ifrac2 || study), data = long_data,
    prior = prior(normal(0, 1), class = sd, coef = "Intercept", group = "study") +
      prior(normal(0, 0.5), class = sd, coef = "ifrac2", group = "study"),
    chains = 4, iter = 20000, warmup = 6000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15))
  saveRDS(fit_subgroup_adj, "fit_edss.rds")
}

# ---- Posterior samples (exponentiated interaction) ----
samples_interaction <- brms::hypothesis(fit_interaction, "Intercept = 0")$samples |> unlist() |> exp()
samples_subgroup <- brms::hypothesis(fit_subgroup, "subgroupB = 0")$samples |> unlist() |> exp()
samples_subgroup_adj <- brms::hypothesis(fit_subgroup_adj, "subgroupB = 0")$samples |> unlist() |> exp()

n <- min(length(samples_interaction), length(samples_subgroup), length(samples_subgroup_adj))
samples_interaction <- samples_interaction[seq_len(n)]
samples_subgroup <- samples_subgroup[seq_len(n)]
samples_subgroup_adj <- samples_subgroup_adj[seq_len(n)]
df_samples <- data.frame(
  value = c(samples_interaction, samples_subgroup, samples_subgroup_adj),
  group = rep(c("Interaction-MA", "Subgroup-MA", "CAMS"), each = n)
)
df_gold <- df_samples |> filter(group == "Interaction-MA")
df_bms  <- df_samples |> filter(group == "Subgroup-MA")
df_cams <- df_samples |> filter(group == "CAMS")

df_before <- bind_rows(
  df_gold |> mutate(panel = "Before", model = "Interaction-MA (interaction-data)"),
  df_bms  |> mutate(panel = "Before", model = "Subgroup-MA (subgroup-data)")
)
df_after <- bind_rows(
  df_gold |> mutate(panel = "After", model = "Interaction-MA (interaction-data)"),
  df_cams |> mutate(panel = "After", model = "CAMS (subgroup-data + IF-adjustment)")
)
df_plot <- bind_rows(df_before, df_after) |>
  mutate(
    panel = factor(panel, levels = c("After", "Before")),
    model = factor(model, levels = c("Interaction-MA (interaction-data)", "Subgroup-MA (subgroup-data)", "CAMS (subgroup-data + IF-adjustment)"))
  )

# ---- Plot: half violins ----
p2 <- ggplot(df_plot, aes(y = value)) +
  geom_half_violin(
    data = df_plot |> filter(model == "Interaction-MA (interaction-data)"),
    aes(x = 0, fill = model), side = "l", alpha = 0.5, draw_quantiles = c(0.05, 0.95), adjust = 3
  ) +
  geom_half_violin(
    data = df_plot |> filter(model == "Subgroup-MA (subgroup-data)"),
    aes(x = 0, fill = model), side = "r", alpha = 0.5, draw_quantiles = c(0.05, 0.95), adjust = 3
  ) +
  geom_half_violin(
    data = df_plot |> filter(model == "CAMS (subgroup-data + IF-adjustment)"),
    aes(x = 0, fill = model), side = "r", alpha = 0.5, draw_quantiles = c(0.05, 0.95), adjust = 3
  ) +
  stat_summary(
    data = df_plot |> filter(model == "Interaction-MA (interaction-data)"),
    aes(x = 0), fun = median, fun.min = ~ quantile(.x, 0.025), fun.max = ~ quantile(.x, 0.975),
    geom = "pointrange", color = "black", size = 0.5, position = position_nudge(x = -0.4)
  ) +
  stat_summary(
    data = df_plot |> filter(model == "Subgroup-MA (subgroup-data)"),
    aes(x = 0), fun = median, fun.min = ~ quantile(.x, 0.025), fun.max = ~ quantile(.x, 0.975),
    geom = "pointrange", color = "#F8766D", size = 0.5, position = position_nudge(x = -0.425)
  ) +
  stat_summary(
    data = df_plot |> filter(model == "CAMS (subgroup-data + IF-adjustment)"),
    aes(x = 0), fun = median, fun.min = ~ quantile(.x, 0.025), fun.max = ~ quantile(.x, 0.975),
    geom = "pointrange", color = "#00BFC4", size = 0.5, position = position_nudge(x = -0.45)
  ) +
  scale_fill_manual(
    name = "Model",
    values = c(
      "Interaction-MA (interaction-data)" = "#FFFFFF",
      "Subgroup-MA (subgroup-data)" = "#F8766D",
      "CAMS (subgroup-data + IF-adjustment)" = "#00BFC4"
    )
  ) +
  labs(x = " ", y = expression("Interaction effect" ~ (gamma))) +
  theme_classic(base_size = 22) +
  theme(
    legend.position = c(0.05, 0.975),
    legend.justification = c("left", "top"),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_flip() +
  xlim(-0.5, 0.65) +
  ylim(0, 3)

print(p2)
ggsave("figures/figure6_lhs.pdf", plot = p2, width = 8, height = 7, units = "in", device = grDevices::pdf)

