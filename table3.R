#
# Table 3: Comparison of methods for choosing information fraction (IF) — EDSS example.
# Matches paper: Overall-IF, Optimal/optimistic IF, Closer to Subgroups, Average (naive),
# Trial-weighted, Clinical Expertise, MAP prevalence. Reproducible (set.seed(1), seed=1 in brm).
#

if (!file.exists("edss_data.csv")) stop("edss_data.csv not found. Set wd to the project folder.")

set.seed(1)

library(brms)
library(dplyr)
library(xtable)

# ---- Data ----
long_data <- read.csv("edss_data.csv", stringsAsFactors = FALSE)
long_data$se <- sqrt(long_data$vi)
long_data$study <- paste0(long_data$study.acronym, long_data$treatment)
long_data$study.name <- long_data$study
long_data$esti <- long_data$yi
long_data$subgroup12 <- as.numeric(factor(long_data$subgroup)) - 1.5
if (!"ifrac2" %in% names(long_data)) {
  long_data$ifrac2 <- rep(c(0, 1), nrow(long_data) / 2) - long_data$ifrac
}

idx_a <- seq(1L, nrow(long_data), by = 2L)
study_ifrac <- long_data$ifrac[idx_a]
study_n     <- long_data$n[idx_a] + long_data$n[idx_a + 1L]

# Shared grid for methods that choose "closest" IF (match conditional effects resolution)
ifrac_grid <- seq(0.10, 0.70, length.out = 201L)

# ---- Load or fit CAMS (fit_edss) ----
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

# ---- Load or fit Subgroup-MA and interaction (for table and Closer to Subgroups) ----
wide_data <- data.frame(
  study.name = long_data$study.name[idx_a],
  contrast.esti = long_data$yi[idx_a + 1L] - long_data$yi[idx_a],
  contrast.se   = sqrt(long_data$vi[idx_a] + long_data$vi[idx_a + 1L])
)

if (file.exists("fit_subgroup_edss.rds")) {
  fit_subgroup <- readRDS("fit_subgroup_edss.rds")
} else {
  fit_subgroup <- brms::brm(
    esti | se(se) ~ 1 + subgroup + (-1 + subgroup12 | study.name),
    data = long_data,
    prior = prior(normal(0, 0.5), class = sd),
    chains = 4, iter = 20000, warmup = 6000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
  saveRDS(fit_subgroup, "fit_subgroup_edss.rds")
}

if (file.exists("fit_interaction_edss.rds")) {
  fit_interaction <- readRDS("fit_interaction_edss.rds")
} else {
  fit_interaction <- brms::brm(
    contrast.esti | se(contrast.se) ~ 1 + (1 | study.name),
    data = wide_data,
    prior = prior(normal(0, 0.5), class = sd),
    chains = 4, iter = 20000, warmup = 6000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
  saveRDS(fit_interaction, "fit_interaction_edss.rds")
}

# ---- Posterior draws ----
draws <- as.data.frame(fit_subgroup_adj)
nms   <- names(draws)
b_Intercept <- draws[[which(nms == "b_Intercept")]]
b_ifrac     <- draws[[which(nms == "b_ifrac")]]
sub_par     <- nms[grepl("subgroup", nms, fixed = TRUE) & startsWith(nms, "b_")][1]
b_subgroupB <- draws[[sub_par]]

# ---- MAP prevalence: Beta(5, 21) prior ----
# Report in table: mean and SD of Beta(5,21).
# RRs: two options:
#   (A) Point: evaluate CAMS at IF = E[Beta(5,21)] = 5/26. Simple, comparable to other rows.
#   (B) Principled: IF is uncertain — sample one IF per CAMS draw from Beta(5,21), compute RR
#       at that IF for each draw, then take quantiles. This propagates IF uncertainty into the
#       RR posterior (wider CrI). We use (B) for the MAP row.
map_ifrac  <- 5 / (5 + 21)
map_sd     <- sqrt((5 * 21) / ((5 + 21)^2 * (5 + 21 + 1)))
n_cams     <- length(b_Intercept)
map_samp   <- rbeta(n_cams, 5, 21)  # one IF per CAMS draw for joint posterior

# CAMS posterior means (for Closer to Subgroups: distance to Subgroup-MA)
m_Intercept <- mean(b_Intercept)
m_ifrac     <- mean(b_ifrac)
m_subgroupB <- mean(b_subgroupB)

# Subgroup-MA (fit_subgroup): posterior medians for Low and High EDSS (log RR)
draws_sg <- as.data.frame(fit_subgroup)
sg_Intercept <- draws_sg[[which(names(draws_sg) == "b_Intercept")]]
sg_par <- names(draws_sg)[grepl("subgroup", names(draws_sg), fixed = TRUE) & startsWith(names(draws_sg), "b_")][1]
sg_subgroupB <- draws_sg[[sg_par]]
subgroup_low_log  <- median(sg_Intercept)
subgroup_high_log <- median(sg_Intercept + sg_subgroupB)

# ---- overall_ifrac (posterior) for Overall-IF row ----
d <- fit_subgroup_adj$data
tau_samp <- brms::hypothesis(fit_subgroup_adj, "study__Intercept = 0", class = "sd")$samples |> unlist()
grp <- (seq_len(nrow(d)) + 1L) %/% 2L
pf  <- d$ifrac - rep(c(0, 1), times = nrow(d) / 2)
den_se <- tapply(d$se^2 * pf^2, grp, sum)
wi  <- sapply(tau_samp, function(t) 1 / (t^2 + den_se))
pi  <- tapply(d$ifrac, grp, `[`, 1)
overall_ifrac_full <- apply(wi, 2, function(w) crossprod(pi, w) / sum(w))

# Overall-IF: report mean (SD) for table; RRs propagate IF uncertainty (one IF per CAMS draw)
overall_ifrac_mean   <- mean(overall_ifrac_full)
overall_ifrac_sd     <- sd(overall_ifrac_full)
# Ensure same length as CAMS draws for joint computation
n_draws <- length(b_Intercept)
if (length(overall_ifrac_full) != n_draws) stop("overall_ifrac_full length must match CAMS draws")

# ---- IF for each method (paper logic) ----
# 1) Optimal: IF in grid that minimizes total CrI width
total_width <- numeric(length(ifrac_grid))
for (i in seq_along(ifrac_grid)) {
  pv <- ifrac_grid[i]
  mean_A <- b_Intercept + b_ifrac * pv
  mean_B <- b_Intercept + b_ifrac * pv + b_subgroupB
  total_width[i] <- diff(quantile(mean_A, c(0.025, 0.975))) + diff(quantile(mean_B, c(0.025, 0.975)))
}
optimal_ifrac <- ifrac_grid[which.min(total_width)]

# 2) Closer to Subgroups: grid IF where CAMS subgroup means (posterior mean) are closest to Subgroup-MA medians (log scale)
cam_low_at_g  <- m_Intercept + m_ifrac * ifrac_grid
cam_high_at_g <- cam_low_at_g + m_subgroupB
dist_to_subgroup <- (cam_low_at_g - subgroup_low_log)^2 + (cam_high_at_g - subgroup_high_log)^2
closer_to_subgroups_ifrac <- ifrac_grid[which.min(dist_to_subgroup)]

# 3) Average (naive): grid IF closest to mean(study-level prevalence)
avg_study_ifrac <- mean(study_ifrac)
average_naive_ifrac <- ifrac_grid[which.min(abs(ifrac_grid - avg_study_ifrac))]

# 4) Trial-weighted: grid IF closest to sum(n*ifrac)/sum(n)
trial_weighted_raw <- sum(long_data$n * long_data$ifrac) / sum(long_data$n)
trial_weighted_ifrac <- ifrac_grid[which.min(abs(ifrac_grid - trial_weighted_raw))]

# 5) Clinical Expertise: grid IF closest to 0.32 (32%)
clinical_ifrac <- ifrac_grid[which.min(abs(ifrac_grid - 0.32))]

# 6) MAP prevalence: Beta(5, 21) mean (SD) for display; RRs integrate over Beta(5,21)

cat("\n[DEBUG] IFs (proportions before formatting):\n")
print(list(
  overall_ifrac_mean   = overall_ifrac_mean,
  overall_ifrac_sd     = overall_ifrac_sd,
  optimal_ifrac        = optimal_ifrac,
  closer_ifrac         = closer_to_subgroups_ifrac,
  avg_study_ifrac      = avg_study_ifrac,
  average_naive_ifrac  = average_naive_ifrac,
  trial_weighted_raw   = trial_weighted_raw,
  trial_weighted_ifrac = trial_weighted_ifrac,
  clinical_ifrac       = clinical_ifrac,
  map_ifrac            = map_ifrac,
  map_sd               = map_sd  # Beta(5,21) SD
))

# ---- Function: at a given IF (scalar), get Low-EDSS, High-EDSS, Overall (RR, exp scale) ----
at_ifrac <- function(ifrac_val) {
  mean_A <- b_Intercept + b_ifrac * ifrac_val
  mean_B <- b_Intercept + b_ifrac * ifrac_val + b_subgroupB
  overall_lin <- b_Intercept + ifrac_val * (b_ifrac + b_subgroupB)
  list(
    low_edss  = quantile(exp(mean_A), c(0.5, 0.025, 0.975)),
    high_edss = quantile(exp(mean_B), c(0.5, 0.025, 0.975)),
    overall   = quantile(exp(overall_lin), c(0.5, 0.025, 0.975))
  )
}

fmt_rr <- function(q) sprintf("%.2f (%.2f; %.2f)", q[1], q[2], q[3])

# ---- Interaction (RRR): same for all methods (from fit_interaction) ----
h_int <- brms::hypothesis(fit_interaction, "Intercept > 0")$hypothesis
int_est <- as.numeric(h_int[c(2, 4, 5)])
interaction_str <- fmt_rr(exp(int_est))

# ---- Build each row ----
# Overall-IF: RR posterior integrates over overall-IF uncertainty (one IF per draw)
mean_A_overall   <- b_Intercept + b_ifrac * overall_ifrac_full
mean_B_overall   <- b_Intercept + b_ifrac * overall_ifrac_full + b_subgroupB
overall_lin_full <- b_Intercept + overall_ifrac_full * (b_ifrac + b_subgroupB)
r_overall        <- list(
  low_edss  = quantile(exp(mean_A_overall),   c(0.5, 0.025, 0.975)),
  high_edss = quantile(exp(mean_B_overall),   c(0.5, 0.025, 0.975)),
  overall   = quantile(exp(overall_lin_full), c(0.5, 0.025, 0.975))
)
r_optimal   <- at_ifrac(optimal_ifrac)
r_closer    <- at_ifrac(closer_to_subgroups_ifrac)
r_naive     <- at_ifrac(average_naive_ifrac)
r_trial     <- at_ifrac(trial_weighted_ifrac)
r_clinical  <- at_ifrac(clinical_ifrac)
# MAP row (principled): RR posterior integrates over Beta(5,21) IF uncertainty
mean_A_map  <- b_Intercept + b_ifrac * map_samp
mean_B_map  <- b_Intercept + b_ifrac * map_samp + b_subgroupB
overall_lin_map <- b_Intercept + map_samp * (b_ifrac + b_subgroupB)
r_map       <- list(
  low_edss  = quantile(exp(mean_A_map),  c(0.5, 0.025, 0.975)),
  high_edss = quantile(exp(mean_B_map),  c(0.5, 0.025, 0.975)),
  overall   = quantile(exp(overall_lin_map), c(0.5, 0.025, 0.975))
)

# ---- Table: Method, IF (%), Low-EDSS, High-EDSS, Overall, Interaction ----
# Footnote *: overall-IF mean (SD); **: Beta(5, 21) mean (SD)
tab <- data.frame(
  Method = c(
    "Overall-IF",
    "Optimal/optimistic IF",
    "Closer to Subgroups",
    "Average (naive)",
    "Trial-weighted",
    "Clinical Expertise",
    "MAP prevalence"
  ),
  IF_pct = c(
    sprintf("%.1f%%*", 100 * overall_ifrac_mean),
    sprintf("%.1f%%", 100 * optimal_ifrac),
    sprintf("%.1f%%", 100 * closer_to_subgroups_ifrac),
    sprintf("%.1f%%", 100 * average_naive_ifrac),
    sprintf("%.1f%%", 100 * trial_weighted_ifrac),
    sprintf("%.1f%%", 100 * clinical_ifrac),
    sprintf("%.1f%%**", 100 * map_ifrac)
  ),
  IF_extra = c(
    sprintf("(%.2f)", 100 * overall_ifrac_sd),
    "",
    "",
    "",
    "",
    "",
    sprintf("(%.1f)", 100 * map_sd)
  ),
  Low_EDSS = c(
    fmt_rr(r_overall$low_edss),
    fmt_rr(r_optimal$low_edss),
    fmt_rr(r_closer$low_edss),
    fmt_rr(r_naive$low_edss),
    fmt_rr(r_trial$low_edss),
    fmt_rr(r_clinical$low_edss),
    fmt_rr(r_map$low_edss)
  ),
  High_EDSS = c(
    fmt_rr(r_overall$high_edss),
    fmt_rr(r_optimal$high_edss),
    fmt_rr(r_closer$high_edss),
    fmt_rr(r_naive$high_edss),
    fmt_rr(r_trial$high_edss),
    fmt_rr(r_clinical$high_edss),
    fmt_rr(r_map$high_edss)
  ),
  Overall = c(
    fmt_rr(r_overall$overall),
    fmt_rr(r_optimal$overall),
    fmt_rr(r_closer$overall),
    fmt_rr(r_naive$overall),
    fmt_rr(r_trial$overall),
    fmt_rr(r_clinical$overall),
    fmt_rr(r_map$overall)
  ),
  Interaction = rep(interaction_str, 7),
  stringsAsFactors = FALSE
)

# ---- Write CSV ----
write.csv(tab, "tables/table3.csv", row.names = FALSE)

# ---- Write LaTeX ----
xtable(tab)

# ---- Print ----
cat("\n-------- Table 3: Methods for choosing information fraction (IF), EDSS example --------\n\n")
cat("Reproducible: set.seed(1), seed=1 in brm. Grid length =", length(ifrac_grid), "\n")
cat("* overall-IF: mean (SD in parentheses); ** Beta(5, 21) mean (SD).\n\n")
print(tab)
cat("\n* overall-IF (mean and standard deviation)\n")
cat("** MAP prevalence: Beta(5, 21) mean and standard deviation.\n")
cat("\nTable 3 saved to table2.csv\n")
