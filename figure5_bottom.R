#
# Figure 5 (bottom): EDSS-based three-panel forest plot (Interaction-MA, Subgroup-MA, CAMS + Overall).
# Self-contained. Requires: edss_data.csv. Uses fit_overall_*.rds and
# reuses fit_subgroup*.rds, fit_interaction*.rds, fit_*.rds (from other figure scripts).
# Needs: brms, ggplot2, dplyr, patchwork, gtable, grid, ckbplotr, metafor, ggtext.
#

# ---- Helper: run one Figure 5 (three-panel plot) for a given dataset ----
run_figure5 <- function(csv_file, prefix, subgroup_names) {
  if (!file.exists(csv_file)) stop(csv_file, " not found. Set wd to the project folder.")
  set.seed(1)
  library(brms)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(gtable)
  library(grid)

  long_data <- read.csv(csv_file, stringsAsFactors = FALSE)
  long_data$se <- sqrt(long_data$vi)
  long_data$study <- paste0(long_data$study.acronym, long_data$treatment)
  long_data$study.name <- long_data$study
  long_data$esti <- long_data$yi
  long_data$subgroup12 <- as.numeric(factor(long_data$subgroup)) - 1.5
  if (!"ifrac2" %in% names(long_data)) {
    long_data$ifrac2 <- rep(c(0, 1), nrow(long_data) / 2) - long_data$ifrac
  }

  idx_a <- seq(1L, nrow(long_data), by = 2L)
  idx_b <- seq(2L, nrow(long_data), by = 2L)
  # Wide: contrast (for fit_interaction) and overall (for fit_overall)
  w_va <- long_data$vi[idx_a]
  w_vb <- long_data$vi[idx_b]
  w_ya <- long_data$yi[idx_a]
  w_yb <- long_data$yi[idx_b]
  w_inv <- 1 / w_va + 1 / w_vb
  wide_data <- data.frame(
    study.name = long_data$study.name[idx_a],
    contrast.esti = w_yb - w_ya,
    contrast.se = sqrt(w_va + w_vb),
    overall.esti = (w_ya / w_va + w_yb / w_vb) / (1 / w_va + 1 / w_vb),
    overall.se = sqrt(1 / w_inv)
  )

  # ---- Fits ----
  cache_subgroup   <- paste0("fit_subgroup", if (prefix != "age") "_edss" else "", ".rds")
  cache_interaction <- paste0("fit_interaction", if (prefix != "age") "_edss" else "", ".rds")
  cache_adj        <- if (prefix == "age") "fit_age.rds" else "fit_edss.rds"
  cache_overall    <- paste0("fit_overall_", prefix, ".rds")

  if (file.exists(cache_subgroup)) {
    fit_subgroup <- readRDS(cache_subgroup)
  } else {
    fit_subgroup <- brms::brm(
      esti | se(se) ~ 1 + subgroup + (-1 + subgroup12 | study.name),
      data = long_data,
      prior = prior(normal(0, 0.5), class = sd),
      chains = 4, iter = 20000, warmup = 6000, seed = 1,
      control = list(adapt_delta = 0.999, max_treedepth = 15)
    )
    saveRDS(fit_subgroup, cache_subgroup)
  }

  if (file.exists(cache_interaction)) {
    fit_interaction <- readRDS(cache_interaction)
  } else {
    fit_interaction <- brms::brm(
      contrast.esti | se(contrast.se) ~ 1 + (1 | study.name),
      data = wide_data,
      prior = prior(normal(0, 0.5), class = sd),
      chains = 4, iter = 20000, warmup = 6000, seed = 1,
      control = list(adapt_delta = 0.999, max_treedepth = 15)
    )
    saveRDS(fit_interaction, cache_interaction)
  }

  if (file.exists(cache_adj)) {
    fit_subgroup_adj <- readRDS(cache_adj)
  } else {
    fit_subgroup_adj <- brms::brm(
      yi | se(se) ~ 1 + ifrac + subgroup + (1 + ifrac2 || study),
      data = long_data,
      prior = prior(normal(0, 1), class = sd, coef = "Intercept", group = "study") +
        prior(normal(0, 0.5), class = sd, coef = "ifrac2", group = "study"),
      chains = 4, iter = 20000, warmup = 6000, seed = 1,
      control = list(adapt_delta = 0.999, max_treedepth = 15)
    )
    saveRDS(fit_subgroup_adj, cache_adj)
  }

  if (file.exists(cache_overall)) {
    fit_overall <- readRDS(cache_overall)
  } else {
    fit_overall <- brms::brm(
      overall.esti | se(overall.se) ~ 1 + (1 | study.name),
      data = wide_data,
      prior = prior(normal(0, 1), class = sd),
      chains = 4, iter = 20000, warmup = 6000, seed = 1,
      control = list(adapt_delta = 0.999, max_treedepth = 15)
    )
    saveRDS(fit_overall, cache_overall)
  }

  # ---- overall_ifrac (from fit_subgroup_adj; group is "study" in project fits) ----
  d <- fit_subgroup_adj$data
  tau <- brms::hypothesis(
    fit_subgroup_adj,
    hypothesis = "study__Intercept = 0",
    class = "sd"
  )$samples |> unlist()
  grp <- (seq_len(nrow(d)) + 1) %/% 2
  pf <- d$ifrac - rep(c(0, 1), times = nrow(d) / 2)
  den_se <- tapply(d$se^2 * pf^2, grp, sum)
  wi <- sapply(tau, function(t) 1 / (t^2 + den_se))
  pi <- tapply(d$ifrac, grp, `[`, 1)
  overall_ifrac <- apply(wi, 2, FUN = function(w) crossprod(pi, w) / sum(w))

  # ---- Summaries (extract_summary logic) ----
  extract_hypothesis <- function(fit, hypothesis_str) {
    brms::hypothesis(fit, hypothesis_str)$hypothesis[c(2, 4, 5)]
  }

  # Subsample overall_ifrac so hypothesis() is not called thousands of times
  n_ifrac <- length(overall_ifrac)
  if (n_ifrac > 4000L) {
    set.seed(1)
    overall_ifrac <- overall_ifrac[sample.int(n_ifrac, 4000L)]
  }

  summary_interaction <- rbind(
    rep(NA_real_, 3),
    as.numeric(extract_hypothesis(fit_interaction, "Intercept > 0")),
    rep(NA_real_, 3),
    as.numeric(extract_hypothesis(fit_subgroup, "subgroupB > 0")),
    rep(NA_real_, 3),
    as.numeric(extract_hypothesis(fit_subgroup_adj, "subgroupB > 0"))
  )
  colnames(summary_interaction) <- c("yi", "ci.lb", "ci.ub")

  subgroupA_samples <- sapply(overall_ifrac, function(x) {
    brms::hypothesis(fit_subgroup_adj, paste0("Intercept + ", x, "* ifrac > 0"))$samples |> unlist()
  })
  subgroupB_samples <- sapply(overall_ifrac, function(x) {
    brms::hypothesis(fit_subgroup_adj, paste0("Intercept + subgroupB + ", x, "* ifrac > 0"))$samples |> unlist()
  })

  summary_subgroup <- rbind(
    rep(NA_real_, 3),
    rep(NA_real_, 3),
    as.numeric(extract_hypothesis(fit_subgroup, "Intercept > 0")),
    as.numeric(extract_hypothesis(fit_subgroup, "Intercept + subgroupB > 0")),
    quantile(subgroupA_samples, probs = c(0.5, 0.025, 0.975)),
    quantile(subgroupB_samples, probs = c(0.5, 0.025, 0.975))
  )
  colnames(summary_subgroup) <- c("yi", "ci.lb", "ci.ub")
  rownames(summary_subgroup) <- rep(subgroup_names, 3)

  overall_samples <- sapply(overall_ifrac, function(x) {
    brms::hypothesis(fit_subgroup_adj, paste0("Intercept + ", x, "* (ifrac + subgroupB) > 0"))$samples |> unlist()
  })
  summary_overall <- rbind(
    rep(NA_real_, 3),
    as.numeric(extract_hypothesis(fit_overall, "Intercept > 0")),
    rep(NA_real_, 3),
    rep(NA_real_, 3),
    rep(NA_real_, 3),
    quantile(overall_samples, probs = c(0.5, 0.025, 0.975))
  )
  colnames(summary_overall) <- c("yi", "ci.lb", "ci.ub")

  re_nm <- names(summary(fit_subgroup_adj)$random)[1L]
  if (is.null(re_nm)) re_nm <- "study"
  summary_tau <- rbind(
    summary(fit_interaction)$random$study.name[c(1, 3, 4)],
    summary(fit_subgroup)$random$study.name[c(1, 3, 4)],
    summary(fit_subgroup_adj)$random[[re_nm]][2, c(1, 3, 4)]
  )
  rownames(summary_tau) <- c("Interaction-MA", "Subgroup-MA", "CAMS")

  # ---- Relabel subgroup for forest (study.name as label) ----
  ld <- long_data
  ld$subgroup <- sapply(ld$subgroup, function(i) {
    switch(as.character(i), A = subgroup_names[1], B = subgroup_names[2], i)
  })

  xlab <- c("Relapse rate ratio (RR)", "Ratio of relapse rate ratios (RRR)")
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

  source("interaction_forest_plot.R", local = TRUE)

  # fp2: Interaction-MA/Subgroup-MA/CAMS subgroup + interaction
  fp2 <- interaction_forest_plot(
    yi = ld$esti,
    vi = ld$se^2,
    labels = ld$study.name,
    pointsize = 10,
    scalepoints = TRUE,
    estimates_only = TRUE,
    base_size = 20,
    xlim   = list(c(0.3, 0.8), c(0.9, 2.1)),
    xticks = list(c(0.3, 0.7, 1), c(1, 1.5, 2)),
    subgroups = ld$subgroup,
    xlab = xlab,
    subgroup_intervals = summary_subgroup[-c(7, 8), , drop = FALSE],
    interaction_intervals = summary_interaction[-c(7, 8), , drop = FALSE],
    model_names = c("Interaction-MA", "Subgroup-MA", "CAMS")
  )

  fp1_sub <- fp2$plots[[1]]$plot +
    theme(plot.margin = margin(t = 40, r = 150, b = 40, l = 40))
  fp2_p2 <- fp2$plots[[2]]$plot +
    theme(plot.margin = margin(t = 40, r = 150, b = 40, l = 40)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  # fp3: overall column
  fp3 <- interaction_forest_plot(
    yi = ld$esti,
    vi = ld$se^2,
    labels = ld$study.name,
    col_right_heading = list("", paste0("Overall \n(95% interval)")),
    pointsize = 10,
    scalepoints = TRUE,
    estimates_only = TRUE,
    base_size = 20,
    xlim   = list(c(0.3, 0.9), c(0.3, 0.9)),
    xticks = list(c(0.3, 0.7, 1), c(0.3, 0.7, 1)),
    subgroups = ld$subgroup,
    xlab = c("", "Relapse rate ratio (RR)"),
    subgroup_intervals = summary_subgroup[-c(7, 8), , drop = FALSE],
    interaction_intervals = summary_overall[-c(7, 8), , drop = FALSE],
    model_names = c("Interaction-MA", "Subgroup-MA", "CAMS"),
    row.labels.space = c(1, 1, 1, 1)
  )
  fp3_p2 <- fp3$plots[[2]]$plot +
    theme(plot.margin = margin(t = 40, r = 150, b = 40, l = 40)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  addtext <- "\n"
  for (i in 1:3) {
    addtext <- paste0(
      addtext,
      sprintf(
        "%s %4.2f(%4.2f; %4.2f)    ",
        rownames(summary_tau)[i],
        summary_tau[i, 1], summary_tau[i, 2], summary_tau[i, 3]
      )
    )
  }
  qs <- round(100 * quantile(overall_ifrac, probs = c(0.5, 0.025, 0.975)), 1)
  addtext <- paste0(
    addtext,
    sprintf(
      "%s %4.1f(%4.1f; %4.1f)    ",
      "Overall-IF (%): ", qs[1], qs[2], qs[3]
    )
  )

  final_plot <- (fp1_sub + fp2_p2 + fp3_p2) +
    plot_layout(guides = "collect", widths = c(1, 1, 1)) +
    labs(caption = bquote(Interaction~heterogeneity(tau[gamma]):~ .(addtext))) +
    theme(
      plot.caption = element_text(hjust = 0.85, vjust = -10, size = 15, face = "italic"),
      plot.caption.position = "plot"
    )

  list(final_plot = final_plot, summary_tau = summary_tau)
}

# ---- Run for EDSS (bottom panel) ----
if (!file.exists("edss_data.csv")) stop("edss_data.csv not found. Set wd to the project folder.")

out_edss <- run_figure5("edss_data.csv", "edss", c("Lower EDSS", "Higher EDSS"))

plot(out_edss$final_plot)

ggsave(
  filename = "figures/figure5_bottom.pdf",
  plot = out_edss$final_plot,
  width = 17,
  height = 5,
  units = "in",
  device = grDevices::pdf
)
