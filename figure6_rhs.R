#
# Figure 6 right: Conditional posterior (gamma | tau) trace plot + tau density.
# Self-contained. Requires: edss_data.csv. Reuses fit_subgroup_edss, fit_interaction_edss, fit_edss.
#

if (!file.exists("edss_data.csv")) stop("edss_data.csv not found. Set wd to the project folder.")

set.seed(1)
library(brms)
library(ggplot2)
library(dplyr)
library(patchwork)

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

# ---- Tau posterior for bottom panel ----
tau_samples <- brms::hypothesis(fit_subgroup_adj, "study__Intercept = 0", class = "sd")$samples |> unlist()
density_tau <- density(tau_samples, adjust = 3)
df_tau <- data.frame(x = density_tau$x, y = density_tau$y)

# ---- Conditional gamma | tau (Vectorized) ----
gamma_interaction <- Vectorize(function(tau) {
  yi <- fit_interaction$data$contrast.esti
  wi <- 1 / (fit_interaction$data$contrast.se^2 + tau^2)
  mu <- (wi %*% yi) / sum(wi)
  sigma <- 1 / sqrt(sum(wi))
  ci.lb <- mu - qnorm(1 - 0.50 / 2) * sigma
  ci.ub <- mu + qnorm(1 - 0.50 / 2) * sigma
  exp(c(mu, ci.lb, ci.ub))
})

gamma_subgroup <- Vectorize(function(tau) {
  y <- fit_subgroup$data$esti
  se <- fit_subgroup$data$se
  n <- length(y) / 2
  sigma2 <- se^2
  tau2_4 <- tau^2 / 4
  V <- diag(sigma2) + kronecker(diag(n), matrix(c(tau2_4, -tau2_4, -tau2_4, tau2_4), 2))
  W <- solve(V)
  X <- cbind(rep(c(1, 0), n), rep(c(0, 1), n))
  mu <- solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% y)
  Sigma <- solve(t(X) %*% W %*% X)
  ci.lb <- diff(mu) - qnorm(1 - 0.50 / 2) * sqrt(sum(diag(Sigma)) - 2 * Sigma[1, 2])
  ci.ub <- diff(mu) + qnorm(1 - 0.50 / 2) * sqrt(sum(diag(Sigma)) - 2 * Sigma[1, 2])
  exp(c(diff(mu), ci.lb, ci.ub))
})

gamma_subgroup_adj <- Vectorize(function(tau) {
  dat <- as.data.frame(fit_subgroup_adj$data)
  y   <- as.numeric(dat$yi)
  se  <- as.numeric(dat$se)
  ifrac_vec <- as.numeric(dat$ifrac)
  n   <- length(y) / 2
  if (n < 1) return(exp(c(NA, NA, NA)))
  sigma2 <- se^2
  tau2 <- tau^2
  heterogeneity <- tau2 * matrix(c(1, 1, 1, 1), 2)
  c_vec <- rep(c(0, 1), n) - ifrac_vec
  C <- diag(c_vec)
  K <- kronecker(diag(n), heterogeneity)
  V <- diag(sigma2) + C %*% K %*% t(C)
  W <- solve(V)
  X <- cbind(rep(c(1, 1), n), ifrac_vec, rep(c(0, 1), n))
  mu <- solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% y)
  Sigma <- solve(t(X) %*% W %*% X)
  ci.lb <- mu[3, 1] - qnorm(1 - 0.50 / 2) * sqrt(Sigma[3, 3])
  ci.ub <- mu[3, 1] + qnorm(1 - 0.50 / 2) * sqrt(Sigma[3, 3])
  exp(c(mu[3, 1], ci.lb, ci.ub))
})

# ---- Grid and curves ----
model_levels <- c("Interaction-MA", "Subgroup-MA", "CAMS")
tau <- seq(0, 1.0, length.out = 100)
data_plot <- data.frame(
  tau = rep(tau, 3),
  gamma = c(gamma_interaction(tau)[1, ], gamma_subgroup(tau)[1, ], gamma_subgroup_adj(tau)[1, ]),
  lower = c(gamma_interaction(tau)[2, ], gamma_subgroup(tau)[2, ], gamma_subgroup_adj(tau)[2, ]),
  upper = c(gamma_interaction(tau)[3, ], gamma_subgroup(tau)[3, ], gamma_subgroup_adj(tau)[3, ]),
  type = factor(rep(model_levels, each = length(tau)), levels = model_levels)
)

# ---- Plots ----
x_min <- 0
x_max <- 1.0
common_x <- scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(mult = 0))

p_trace <- ggplot(data_plot, aes(x = tau, y = gamma, color = type, linetype = type, linewidth = type)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.05, linetype = "blank") +
  labs(
    x = NULL,
    y = expression("Conditional interaction  (" ~ gamma ~ "|" ~ tau[gamma] ~ ")"),
    title = "Interaction trace plot"
  ) +
  scale_color_manual(
    values = c("Interaction-MA" = "black", "Subgroup-MA" = "#F8766D", "CAMS" = "#00BFC4"),
    limits = model_levels, name = "Model"
  ) +
  scale_fill_manual(
    values = c("Interaction-MA" = "black", "Subgroup-MA" = "#F8766D", "CAMS" = "#00BFC4"),
    limits = model_levels, name = "Model"
  ) +
  scale_linetype_discrete(name = "Model") +
  scale_linewidth_manual(
    values = c("Interaction-MA" = 3, "Subgroup-MA" = 1.5, "CAMS" = 1.5),
    limits = model_levels, name = "Model"
  ) +
  scale_y_continuous(breaks = c(1, 1.25, 1.5, 1.75)) +
  coord_cartesian(xlim = c(x_min, x_max), ylim = c(1, 1.75)) +
  common_x +
  theme_classic(base_size = 22) +
  theme(
    legend.position = c(0.05, 0.975),
    legend.justification = c("left", "top"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.grid = element_blank()
  )

p_dens <- ggplot(df_tau, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = 0, ymax = y), fill = "grey50", alpha = 0.7) +
  geom_line(linewidth = 0.8, color = "black", alpha = 0.7) +
  labs(x = expression("Heterogeneity (" ~tau[gamma]~")"), y = "") +
  common_x +
  coord_cartesian(xlim = c(x_min, x_max)) +
  theme_classic(base_size = 22) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
    panel.grid = element_blank()
  )

p <- p_trace / p_dens + plot_layout(heights = c(2, 1))
print(p)
ggsave("figures/figure6_rhs.pdf", plot = p, width = 8, height = 7, units = "in", device = grDevices::pdf)
