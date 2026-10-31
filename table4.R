# Frequentist counterparts of Interaction-MA, Subgroup-MA, CAMS, and
# the within-trial (Stata metafloat) framework, on the Age and EDSS datasets.
# Estimates on RR / RRR scale with 95% CIs.
#
# CAMS follows CAMS/compare_bayesian_frequentist_edss_overall_if.R
# (metafor rma.mv with struct = "GEN", V = vi; optional .rds cache).
# Within-trial uses Stata metafloat REML (wt_metafloat_export.do).
#
# Run from this folder:  Rscript table4.R

if (!file.exists("age_data.csv")) stop("age_data.csv not found. Set wd to the project folder.")

suppressPackageStartupMessages({
  library(metafor)
})

zcrit <- 1.96
fmt <- function(est, lcl, ucl, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f; %.", digits, "f)"), est, lcl, ucl)
}
wald_exp <- function(est, se, z = zcrit) {
  c(est = exp(est), lcl = exp(est - z * se), ucl = exp(est + z * se))
}
lincom_exp <- function(L, b, V, z = zcrit) {
  est <- drop(L %*% b)
  se  <- sqrt(drop(L %*% V %*% L))
  wald_exp(est, se, z = z)
}

prep_data <- function(csv) {
  d <- read.csv(csv, stringsAsFactors = FALSE)
  d$study <- paste0(d$study.acronym, d$treatment)
  d$subgroup <- factor(d$subgroup, levels = c("A", "B"))
  d$x <- as.integer(d$subgroup == "B")
  d$subgroup12 <- as.numeric(d$subgroup) - 1.5
  if (!"ifrac2" %in% names(d)) {
    d$ifrac2 <- d$x - d$ifrac
  }
  d <- d[order(d$study, d$subgroup), ]
  rownames(d) <- NULL
  d
}

write_wt_input <- function(d, file) {
  wt_in <- data.frame(
    study    = as.integer(factor(d$study, levels = unique(d$study))),
    subgroup = as.integer(d$subgroup),
    ES       = d$yi,
    seES     = sqrt(d$vi)
  )
  wt_in <- wt_in[order(wt_in$study, wt_in$subgroup), ]
  write.csv(wt_in, file, row.names = FALSE)
}

read_wt_results <- function(file) {
  res <- read.csv(file, stringsAsFactors = FALSE)
  rownames(res) <- res$quantity
  b1 <- res["beta1", "estimate"]
  b2 <- res["beta2", "estimate"]
  v11 <- res["beta1", "se"]^2
  v22 <- res["beta2", "se"]^2
  vg  <- res["gamma", "se"]^2
  # gamma = beta2 - beta1 => Cov(b1,b2) = (v11 + v22 - Var(gamma)) / 2
  v12 <- (v11 + v22 - vg) / 2
  V <- matrix(c(v11, v12, v12, v22), 2, 2)
  ones <- c(1, 1)
  Vinv <- solve(V)
  th <- drop(ones %*% Vinv %*% c(b1, b2) / drop(ones %*% Vinv %*% ones))
  se_th <- sqrt(1 / drop(ones %*% Vinv %*% ones))
  list(
    A   = c(est = exp(res["beta1", "estimate"]),
            lcl = exp(res["beta1", "ci_lb"]),
            ucl = exp(res["beta1", "ci_ub"])),
    B   = c(est = exp(res["beta2", "estimate"]),
            lcl = exp(res["beta2", "ci_lb"]),
            ucl = exp(res["beta2", "ci_ub"])),
    int = c(est = exp(res["gamma", "estimate"]),
            lcl = exp(res["gamma", "ci_lb"]),
            ucl = exp(res["gamma", "ci_ub"])),
    overall = wald_exp(th, se_th)
  )
}

run_stata_metafloat <- function() {
  needed <- c("wt_stata_age_results.csv", "wt_stata_edss_results.csv")
  if (all(file.exists(needed))) {
    cat("Using existing Stata metafloat results.\n")
    return(invisible(TRUE))
  }
  stata <- find_stata()
  if (is.na(stata)) {
    stop("Stata not found. Install Stata or add it to PATH, then re-run.")
  }
  unlink("wt_metafloat_export.log")
  status <- system2(stata, args = c("-b", "do", "wt_metafloat_export.do"))
  ok <- is.null(status) || identical(as.integer(status), 0L)
  if (!ok) {
    stop("Stata metafloat failed (exit ", status, "). See wt_metafloat_export.log")
  }
  needed <- c("wt_stata_age_results.csv", "wt_stata_edss_results.csv")
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    stop("Stata ran but did not write: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

fit_cams_metafor <- function(d, cache_file = NULL) {
  # Same specification as CAMS/compare_bayesian_frequentist_edss_overall_if.R
  if (!is.null(cache_file) && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  fit <- metafor::rma.mv(
    yi,
    V = vi,
    mods = ~ ifrac + subgroup,
    random = list(~ 1 | study, ~ 0 + ifrac2 | study),
    data = d,
    struct = "GEN",
    method = "REML"
  )
  if (!is.null(cache_file)) saveRDS(fit, cache_file)
  fit
}

analyse <- function(d, label, levels_A_B, cams_cache = NULL, wt, latex_group = label) {
  idx_a <- which(d$subgroup == "A")
  idx_b <- which(d$subgroup == "B")
  stopifnot(length(idx_a) == length(idx_b))
  stopifnot(identical(d$study[idx_a], d$study[idx_b]))
  
  cat("\n====================\n", label, "\n====================\n")
  cat("k studies =", length(idx_a), "\n")
  
  # ---- 1. Interaction-MA: RE MA of within-trial contrasts ----
  contrast <- d$yi[idx_b] - d$yi[idx_a]
  sei_c    <- sqrt(d$vi[idx_a] + d$vi[idx_b])
  fit_int  <- rma(yi = contrast, sei = sei_c, method = "REML", test = "z")
  int_ma   <- wald_exp(as.numeric(fit_int$b), fit_int$se)
  cat("\nInteraction-MA (REML)\n")
  print(round(int_ma, 4))
  cat("  tau2_gamma =", round(fit_int$tau2, 4), "\n")
  
  # ---- 2. Subgroup-MA: intercept + subgroup, RE on centered contrast ----
  # Matches BMS: ~ 1 + subgroup + (-1 + subgroup12 | study)
  # struct = "GEN" treats subgroup12 as a quantitative slope (rho = -1),
  # same as brms; default CS would treat the two +/- 0.5 levels as a factor.
  fit_sg <- rma.mv(
    yi, vi,
    mods = ~ subgroup,
    random = ~ 0 + subgroup12 | study,
    struct = "GEN",
    data = d,
    method = "REML",
    test = "z"
  )
  b_sg <- coef(fit_sg)
  V_sg <- vcov(fit_sg)
  sg_A   <- lincom_exp(c(1, 0), b_sg, V_sg)
  sg_B   <- lincom_exp(c(1, 1), b_sg, V_sg)
  sg_int <- lincom_exp(c(0, 1), b_sg, V_sg)
  cat("\nSubgroup-MA (REML, struct = GEN)\n")
  print(round(rbind(A = sg_A, B = sg_B, interaction = sg_int), 4))
  cat("  tau2_slope =", round(as.numeric(fit_sg$G), 4), "\n")
  
  # ---- 3. CAMS: same frequentist analogue as compare_bayesian_frequentist_edss_overall_if.R
  # brms: ~ 1 + ifrac + subgroup + (1 + ifrac2 || study)
  fit_cams <- fit_cams_metafor(d, cache_file = cams_cache)
  b_c <- coef(fit_cams)
  V_c <- vcov(fit_cams)
  ci_fixed <- confint(fit_cams, fixed = TRUE)$fixed
  z_cams <- qnorm(0.975)
  tau2_m <- fit_cams$sigma2[1]
  pi_j   <- d$ifrac[idx_a]
  den_se <- (1 - pi_j)^2 * d$vi[idx_a] + pi_j^2 * d$vi[idx_b]
  w_j    <- 1 / (tau2_m + den_se)
  pi_star <- drop(crossprod(w_j, pi_j) / sum(w_j))
  pi_trial_raw <- sum(d$n * d$ifrac) / sum(d$n)
  # Same IF-selection logic as CAMS/table2.R (grid snap where that script snaps)
  ifrac_grid <- seq(0.10, 0.70, length.out = 201L)
  snap <- function(x) ifrac_grid[which.min(abs(ifrac_grid - x))]
  se_lin <- function(L) sqrt(drop(L %*% V_c %*% L))
  total_width <- vapply(ifrac_grid, function(pv) {
    2 * z_cams * (se_lin(c(1, pv, 0)) + se_lin(c(1, pv, 1)))
  }, numeric(1))
  pi_optimal <- ifrac_grid[which.min(total_width)]
  cam_A_g <- b_c[[1]] + ifrac_grid * b_c[[2]]
  cam_B_g <- cam_A_g + b_c[[3]]
  dist_sg <- (cam_A_g - b_sg[[1]])^2 + (cam_B_g - (b_sg[[1]] + b_sg[[2]]))^2
  pi_closer  <- ifrac_grid[which.min(dist_sg)]
  pi_naive   <- snap(mean(pi_j))
  pi_trial   <- snap(pi_trial_raw)
  pi_clinical <- snap(0.32)
  pi_map     <- 5 / (5 + 21)
  if_methods <- list(
    "Overall-IF"              = pi_star,
    "Optimal/optimistic IF"   = pi_optimal,
    "Closer to Subgroups"     = pi_closer,
    "Average (naive)"         = pi_naive,
    "Trial-weighted"          = pi_trial,
    "Clinical Expertise"      = pi_clinical,
    "MAP prevalence"          = pi_map
  )
  
  # Remove EDSS-specific assumptions from Age only
  if (label == "Age (<40 vs >=40)") {
    if_methods[c("Clinical Expertise", "MAP prevalence")] <- NULL
  }
  
  cat("\nCAMS (REML, struct = GEN)  tau2_m =", round(tau2_m, 4), "\n")
  cat("  IFs (%):", paste(sprintf("%s=%.1f", names(if_methods), 100 * unlist(if_methods)), collapse = "; "), "\n")
  cat("  fixed-effect CIs (confint, then exp), as in compare script:\n")
  print(round(exp(as.matrix(ci_fixed)), 4))
  cams_at <- function(pi) {
    list(
      A  = lincom_exp(c(1, pi, 0), b_c, V_c, z = z_cams),
      B  = lincom_exp(c(1, pi, 1), b_c, V_c, z = z_cams),
      ov = lincom_exp(c(1, pi, pi), b_c, V_c, z = z_cams)
    )
  }
  # Interaction is the subgroup coefficient (does not depend on IF)
  sg_row <- grep("subgroup", rownames(ci_fixed), ignore.case = TRUE)[1]
  cams_int <- c(
    est = exp(ci_fixed[sg_row, "estimate"]),
    lcl = exp(ci_fixed[sg_row, "ci.lb"]),
    ucl = exp(ci_fixed[sg_row, "ci.ub"])
  )
  cams_by_if <- lapply(if_methods, function(pi) c(cams_at(pi), list(int = cams_int, pi = pi)))
  for (nm in names(cams_by_if)) {
    cat("  CAMS at", nm, sprintf("(%.1f%%)", 100 * cams_by_if[[nm]]$pi), "\n")
    x <- cams_by_if[[nm]]
    print(round(rbind(A = x$A, B = x$B, overall = x$ov, interaction = x$int), 4))
  }
  cat("  CAMS interaction vs Interaction-MA |diff| =",
      round(abs(unname(cams_int["est"]) - unname(int_ma["est"])), 6), "\n")
  
  # ---- 4. Within-trial framework: Stata metafloat (REML) ----
  cat("\nWithin-trial (Stata metafloat, REML)\n")
  print(round(rbind(A = wt$A, B = wt$B, interaction = wt$int), 4))
  
  out <- list(
    label = label,
    levels = levels_A_B,
    latex_group = latex_group,
    pi_star = pi_star,
    interaction_ma = int_ma,
    subgroup_ma = list(A = sg_A, B = sg_B, int = sg_int),
    cams_by_if = cams_by_if,
    wt = list(A = wt$A, B = wt$B, int = wt$int)
  )
  invisible(out)
}

age_d  <- prep_data("age_data.csv")
edss_d <- prep_data("edss_data.csv")
write_wt_input(age_d,  "wt_metafloat_age_input.csv")
write_wt_input(edss_d, "wt_metafloat_edss_input.csv")
run_stata_metafloat()
wt_age  <- read_wt_results("wt_stata_age_results.csv")
wt_edss <- read_wt_results("wt_stata_edss_results.csv")

age  <- analyse(
  age_d, "Age (<40 vs >=40)", c("Younger", "Older"),
  cams_cache = "fit_age_metafor.rds",
  wt = wt_age,
  latex_group = "Age (subgroup B: $\\geq 40$)"
)
edss <- analyse(
  edss_d, "EDSS (<=3.5 vs >3.5)", c("Lower EDSS", "Higher EDSS"),
  cams_cache = "fit_edss_metafor.rds",
  wt = wt_edss,
  latex_group = "EDSS (subgroup B: $>3.5$)"
)

row_int <- function(x) fmt(x["est"], x["lcl"], x["ucl"])
fmt_if <- function(p) sprintf("%.1f\\\\%%", 100 * p)
dash <- "--"

cat("\n\n===== LaTeX caption =====\n")
cat(paste0(
  "\\caption{Frequentist analyses for the age and EDSS subgroup comparisons.\n",
  "Estimates are relapse-rate ratios (RR) for subgroups and the overall effect, and ratios of RRs (RRR) for the interaction, with 95\\% Wald confidence intervals.\n",
  "Between-study heterogeneity was estimated by restricted maximum likelihood (REML).\n",
  "The overall effect is available only from CAMS, evaluated at each information fraction (IF).\n",
  "Subgroup A is age ${<}40$ or EDSS ${\\leq}3.5$.}\n"
))

cat("\n===== LaTeX rows =====\n")
for (res in list(age, edss)) {
  cat(res$latex_group, "\n", sep = "")
  cat("% A=", res$levels[1], "  B=", res$levels[2], "\n", sep = "")
  cat(" & Interaction-MA         & ", dash, " & ", dash, " & ", dash, " & ",
      row_int(res$interaction_ma), " \\\\\n", sep = "")
  cat(" & Subgroup-MA            & ",
      row_int(res$subgroup_ma$A), " & ",
      row_int(res$subgroup_ma$B), " & ", dash, " & ",
      row_int(res$subgroup_ma$int), " \\\\\n", sep = "")
  for (nm in names(res$cams_by_if)) {
    x <- res$cams_by_if[[nm]]
    cat(" & CAMS (", nm, " ", fmt_if(x$pi), ") & ",
        row_int(x$A), " & ",
        row_int(x$B), " & ",
        row_int(x$ov), " & ",
        row_int(x$int), " \\\\\n", sep = "")
  }
  cat(" & Within-trial framework & ",
      row_int(res$wt$A), " & ",
      row_int(res$wt$B), " & ", dash, " & ",
      row_int(res$wt$int), " \\\\\n", sep = "")
}


# ---- Save frequentist results as CSV ----

make_csv_rows <- function(res) {
  
  add_row <- function(method, IF = NA_real_,
                      A = NULL, B = NULL,
                      overall = NULL, interaction = NULL) {
    
    value <- function(x) {
      if (is.null(x)) NA_character_ else row_int(x)
    }
    
    data.frame(
      Dataset = res$label,
      Method = method,
      IF_percent = IF,
      Subgroup_A = value(A),
      Subgroup_B = value(B),
      Overall = value(overall),
      Interaction = value(interaction)
    )
  }
  
  rows <- list(
    add_row("Interaction-MA",
            interaction = res$interaction_ma),
    
    add_row("Subgroup-MA",
            A = res$subgroup_ma$A,
            B = res$subgroup_ma$B,
            interaction = res$subgroup_ma$int)
  )
  
  for (nm in names(res$cams_by_if)) {
    x <- res$cams_by_if[[nm]]
    
    rows[[length(rows) + 1L]] <- add_row(
      paste0("CAMS (", nm, ")"),
      IF = 100 * x$pi,
      A = x$A,
      B = x$B,
      overall = x$ov,
      interaction = x$int
    )
  }
  
  rows[[length(rows) + 1L]] <- add_row(
    "Within-trial framework",
    A = res$wt$A,
    B = res$wt$B,
    interaction = res$wt$int
  )
  
  do.call(rbind, rows)
}

csv_results <- do.call(
  rbind,
  lapply(list(age, edss), make_csv_rows)
)

write.csv(
  csv_results,
  "tables/table4.csv",
  row.names = FALSE,
  na = ""
)

cat("\nSaved: table4.csv\n")
