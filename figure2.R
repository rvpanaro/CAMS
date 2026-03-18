#
# Figure 2: Interaction forest plot (edss_data).
# Same as figure1 but for EDSS. Uses only fit_subgroup and fit_interaction.
# Caches: fit_subgroup_edss.rds, fit_interaction_edss.rds.
#

if (!file.exists("edss_data.csv")) stop("edss_data.csv not found. Set wd to the project folder.")

set.seed(1)
library(brms)
library(ggplot2)
library(dplyr)
library(gtable)
library(grid)

# Read data
long_data <- read.csv("edss_data.csv", stringsAsFactors = FALSE)
long_data$se <- sqrt(long_data$vi)
long_data$study <- paste0(long_data$study.acronym, long_data$treatment)
long_data$study.name <- long_data$study
long_data$esti <- long_data$yi
long_data$subgroup12 <- as.numeric(factor(long_data$subgroup)) - 1.5

# Build wide_data (one row per study) for fit_interaction: contrast = B - A
idx_a <- seq(1L, nrow(long_data), by = 2L)
idx_b <- seq(2L, nrow(long_data), by = 2L)
wide_data <- data.frame(
  study.name = long_data$study.name[idx_a],
  contrast.esti = long_data$yi[idx_b] - long_data$yi[idx_a],
  contrast.se = sqrt(long_data$vi[idx_a] + long_data$vi[idx_b])
)

# Fit or load fit_subgroup (subgroup-MA)
cache_subgroup <- "fit_subgroup_edss.rds"
if (file.exists(cache_subgroup)) {
  fit_subgroup <- readRDS(cache_subgroup)
} else {
  fit_subgroup <- brms::brm(
    esti | se(se) ~ 1 + subgroup + (-1 + subgroup12 | study.name),
    data = long_data,
    prior = prior(normal(0, 0.5), class = sd),
    chains = 4, iter = 12000, warmup = 4000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
  saveRDS(fit_subgroup, cache_subgroup)
}

# Fit or load fit_interaction (interaction-MA)
cache_interaction <- "fit_interaction_edss.rds"
if (file.exists(cache_interaction)) {
  fit_interaction <- readRDS(cache_interaction)
} else {
  fit_interaction <- brms::brm(
    contrast.esti | se(contrast.se) ~ 1 + (1 | study.name),
    data = wide_data,
    prior = prior(normal(0, 0.5), class = sd),
    chains = 4, iter = 12000, warmup = 4000, seed = 1,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
  saveRDS(fit_interaction, cache_interaction)
}

# Extract hypothesis (Estimate, CI.Lower, CI.Upper)
extract_hypothesis <- function(fit, hypothesis_str) {
  brms::hypothesis(fit, hypothesis_str)$hypothesis[c(2, 4, 5)]
}

# Subgroup intervals: 2 rows (Lower EDSS, Higher EDSS) from fit_subgroup
h_a <- extract_hypothesis(fit_subgroup, "Intercept > 0")
h_b <- extract_hypothesis(fit_subgroup, "Intercept + subgroupB > 0")
subgroup_intervals <- matrix(
  c(as.numeric(h_a), as.numeric(h_b)),
  nrow = 2, ncol = 3, byrow = TRUE,
  dimnames = list(c("Lower EDSS", "Higher EDSS"), c("yi", "ci.lb", "ci.ub"))
)

# Interaction interval: from fit_interaction
h_int <- extract_hypothesis(fit_interaction, "Intercept > 0")
interaction_intervals <- rbind(
  rep(NA_real_, 3),
  as.numeric(h_int)
)
colnames(interaction_intervals) <- c("yi", "ci.lb", "ci.ub")

# Labels and subgroups (match example_analysis for EDSS)
ld <- long_data
ld$subgroup_label <- sapply(ld$subgroup, function(i) switch(as.character(i), A = "Lower EDSS", B = "Higher EDSS", i))
labels <- ld$study.acronym
labels[5:6] <- "CONFIRM-B"
labels[7:8] <- "CONFIRM-C"
labels[15:16] <- "COP1 PILOT"
subgroups <- ld$subgroup_label
subgroup.total <- ld$n
total <- with(ld, rep(n[c(TRUE, FALSE)] + n[c(FALSE, TRUE)], each = 2))

# Source forest plot function (in same folder)
source("interaction_forest_plot.R")

# Theme and xlab (RR for relapse)
theme_set(
  theme_classic(base_size = 15, base_family = "sans") +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 20, face = "bold"),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.key = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray85"),
      panel.grid.minor = element_blank()
    )
)
xlab <- c("Relapse rate ratio (RR)", "Ratio of relapse rate ratios (RRR)")

# Build forest plot (fp1 style: one model row, prevalence)
fp1 <- interaction_forest_plot(
  yi = ld$yi,
  vi = ld$vi,
  labels = labels,
  scalepoints = TRUE,
  subgroups = subgroups,
  xlim   = list(c(0.2, 1.2), c(0.75, 2.2)),
  xticks = list(c(0.3, 0.6, 1), c(1, 1.5, 2)),
  xlab = xlab,
  subgroup_intervals = subgroup_intervals,
  interaction_intervals = interaction_intervals,
  model_names = c(""),
  subgroup.total = subgroup.total,
  total = total,
  row.labels.space = c(2, 2, 2, 2)
)

fp1 <- gtable::gtable_add_padding(fp1$figure, padding = unit(c(10, 0, 0, 0), "mm"))
plot(fp1)

# Optional: save
# ggsave("figure2.pdf", plot = fp1, width = 17, height = 10, units = "in")
