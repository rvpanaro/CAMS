# Compare Bayesian (brms) vs Frequentist (metafor) CAMS - EDSS
#
# Focus: Compare EDSS fixed-effect coefficients (brms vs metafor) and the interaction (RRR).
#
# Usage (from this folder):
#   Rscript compare_bayesian_frequentist_edss_overall_if.R
#
# Notes:
# - The script caches model fits as local .rds files on first run.
# - The repository .gitignore excludes *.rds from version control.

if (!file.exists("edss_data.csv")) {
  stop("edss_data.csv not found. Run this script from the CAMS repository folder.")
}

set.seed(1)

require_brms <- requireNamespace("brms", quietly = TRUE)
require_metafor <- requireNamespace("metafor", quietly = TRUE)

if (!require_brms) stop("Missing package: brms")
if (!require_metafor) stop("Missing package: metafor")

# -----------------------------
# Data
# -----------------------------
long_data <- read.csv("edss_data.csv", stringsAsFactors = FALSE)
long_data$se <- sqrt(long_data$vi)
long_data$study <- paste0(long_data$study.acronym, long_data$treatment)

# Ensure subgroup coding matches brms' default A-reference / B coefficient.
long_data$subgroup <- factor(long_data$subgroup, levels = c("A", "B"))

# -----------------------------
# Bayesian: brms CAMS fit (fit_edss)
# -----------------------------
fit_file_bayes <- "fit_edss.rds"

if (file.exists(fit_file_bayes)) {
  fit_bayes <- readRDS(fit_file_bayes)
} else {
  # Match CAMS/figure3_rhs.R model settings.
  fit_bayes <- brms::brm(
    yi | se(se) ~ 1 + ifrac + subgroup + (1 + ifrac2 || study),
    data = long_data,
    prior = brms::prior(normal(0, 1), class ="sd", coef = "Intercept", group = "study") +
      brms::prior(normal(0, 0.5), class = "sd", coef = "ifrac2", group = "study"),
    chains = 4,
    iter = 20000,
    warmup = 6000,
    seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
  saveRDS(fit_bayes, fit_file_bayes)
}

# -----------------------------
# Frequentist: metafor rma.mv CAMS analogue
# -----------------------------
fit_file_freq <- "fit_edss_metafor.rds"

if (file.exists(fit_file_freq)) {
  fit_freq <- readRDS(fit_file_freq)
} else {
  # Random effects emulate (1 + ifrac2 || study) from brms:
  # - (||) in brms is diagonal covariance.
  fit_freq <- metafor::rma.mv(
    yi,
    V = vi,
    mods = ~ ifrac + subgroup,
    random = list(~ 1 | study, ~ 0 + ifrac2 | study),
    data = long_data,
    struct = "GEN",
    method = "REML"
  )
  saveRDS(fit_freq, fit_file_freq)
}

# brms: posterior intervals for fixed effects
bayes <- brms::fixef(fit_bayes, probs = c(0.025, 0.975)) |>
  exp() |>
  round(2) |>
  (\(x) x[, -2])()

# metafor: Wald CIs for coefficients
freq <- confint(fit_freq, fixed = TRUE)$fixed |>
  exp() |> 
  round(2)

print("intercept")
print(rbind(freq = freq[1,], bayes  = bayes[1,]))
print("ifrac")
print(rbind(freq = freq[2,], bayes  = bayes[2,]))
print("subgroup")
print(rbind(freq = freq[3,], bayes  = bayes[3,]))

