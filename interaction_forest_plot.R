# interaction_forest_plot.R
#
#    Copyright (C) 2024  Renato Panaro
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
# R code to generate interaction forest plots for subgroup analyses,
# visualizing effect size estimates and confidence intervals.
# Project copy: same as R/interaction_forest_plot.R for self-contained reproduction.
#
#   Authors: Renato Panaro
#   Last Updated: December 2024
#
#' Generate Interaction Forest Plots
#'
#' This function creates forest plots to visualize effect size estimates and confidence intervals
#' for subgroup analyses, including interaction effects.
#' The function relies on the `ckbplotr` package for forest plot visualization and the `metafor`
#' package for effect size and variance calculations.
#'
#' @param yi Numeric vector of effect size estimates.
#' @param vi Numeric vector of variances for the effect size estimates.
#' @param labels Character vector of labels for each subgroup.
#' @param subgroups Character vector indicating subgroup membership.
#' @param subgroup_intervals Data frame containing confidence intervals for subgroups.
#' @param interaction_intervals Data frame containing confidence intervals for interaction effects.
#' @param model_names Character vector of model names.
#' @param xlab Character string for x-axis label.
#' @param exponentiate Logical; whether to exponentiate effect sizes.
#' @param col_right_heading Character vector for right column headings.
#' @param xlim List or vector defining x-axis limits per panel.
#' @param xticks List or vector of x-axis tick positions.
#' @param ... Additional arguments passed to `ckbplotr::forest_plot`.
#'
#' @return A forest plot object (from `ckbplotr`).
#'
#' @import ckbplotr
#' @import metafor
#' @import dplyr
#' @import magrittr
#' @importFrom grid unit
#' @importFrom ggplot2 margin
#'

interaction_forest_plot <-
  function(yi, vi, labels, subgroups,
           subgroup_intervals, interaction_intervals,
           model_names = "Model 1",
           treat.events = "",
           treat.total = "",
           cont.events = "",
           cont.total = "",
           subgroup.total = "",
           total = "",
           reference_group = 1,
           xlab, exponentiate = TRUE, col_right_heading,
           stroke = 3,
           digits = 2,
           fill = "white",
           study_fill = NULL,
           study_colour = NULL,
           summary_fill = NULL,
           summary_colour = NULL,
           pointsize = 13, base_size = 20, prob = 0.95,
           left.space = unit(30, "mm"),
           xlim = c(0.125, 8),
           xticks = c(0.125, 1, 8),
           addaes = c("identity", "log"),
           grouping,
           row.labels.space = c(0,1,0,0),
           estimates_only = FALSE,
           which_panel = NULL,
           ...) {

    group_labels = rep("", length(labels))

    library(ckbplotr)
    library(metafor)
    library(dplyr)
    library(magrittr)
    library(grid)
    library(ggplot2)
    library(ggtext)
    addaes <- match.arg(addaes)

    addaes <- switch(
      addaes,
      log = list(
        paste0(
          "stroke = stroke / 10 * ((2 * 1.96) / ",
          "(log(uci_transformed) - log(lci_transformed)))^(1/2)"
        )
      ),
      identity = list(
        paste0(
          "stroke = stroke / 10 * ((2 * 1.96) / ",
          "(uci_transformed - lci_transformed))^(1/2)"
        )
      )
    )

    stopifnot(!is.null(rownames(subgroup_intervals)))
    num_subgroups <- length(unique(rownames(subgroup_intervals)))
    k <- nrow(subgroup_intervals)

    row_labels <- data.frame(
      subheading = labels,
      heading = group_labels,
      label = subgroups
    )

    n <- nrow(row_labels)
    additional_rows <- num_subgroups * length(model_names)

    row_labels[(n + 1):(n + additional_rows), ] <- cbind(
      rep(model_names, each = num_subgroups),
      sub("\\|.*", "", rownames(subgroup_intervals)),
      sub(".*\\|", "", rownames(subgroup_intervals))
    )

    resultsA <- metafor::escalc(yi = yi, vi = vi) %>%
      summary(level = 100 * prob, transf = identity)
    resultsA[resultsA$ci.ub > 195, ] <- NA
    resultsA[(n + 1):(n + additional_rows), ] <- subgroup_intervals

    resultsA$fill <- resultsA$colour <- NA
    study_f <- if (!is.null(study_fill)) study_fill else "black"
    study_c <- if (!is.null(study_colour)) study_colour else "black"
    sum_f   <- if (!is.null(summary_fill)) summary_fill else "white"
    sum_c   <- if (!is.null(summary_colour)) summary_colour else "royalblue3"
    resultsA$fill[row_labels$subheading %in% labels] <- study_f
    resultsA$colour[row_labels$subheading %in% labels] <- study_c
    resultsA$fill[!(row_labels$subheading %in% labels)] <- sum_f
    resultsA$colour[!(row_labels$subheading %in% labels)] <- sum_c
    resultsA$bold <- resultsA$colour != "black"
    resultsA$colour[resultsA$bold] <- sum_c

    resultsA$variable <- row_labels$variable <- apply(row_labels, 1, paste0, collapse = "")
    resultsB <- resultsA

    compute_interactions <- function(grouped, reference = 1, point = TRUE) {
      if (!all(c("endpoints", "yi", "vi") %in% colnames(grouped))) {
        stop("Input data 'grouped' must contain 'endpoints', 'yi', and 'vi' columns.")
      }
      unique_subgroups <- unique(subgroups)
      if (reference < 1 || reference > length(unique_subgroups)) {
        stop("Reference group index is out of bounds.")
      }
      C <- diag(length(unique_subgroups))
      C[, reference] <- C[, reference] - 1
      idx <- match(grouped$endpoints, unique_subgroups)
      if (any(is.na(idx))) {
        stop("Some endpoints in 'grouped' do not match the subgroups.")
      }
      if (!point) {
        C <- abs(C)
        result <- C[idx, idx] %*% grouped$vi
      } else {
        result <- C[idx, idx] %*% grouped$yi
      }
      return(result)
    }

    data_subset <- data.frame(
      yi = yi,
      vi = vi,
      labels = labels,
      endpoints = subgroups,
      groups = group_labels
    )
    grouped_data <- data_subset %>%
      group_by(groups, labels)
    interaction_results_yi <- grouped_data %>%
      do(data.frame(yi = compute_interactions(., reference = reference_group)))
    interaction_results_vi <- grouped_data %>%
      do(data.frame(vi = compute_interactions(., point = FALSE, reference = reference_group)))
    matched_yi <- interaction_results_yi %>%
      .[pmatch(labels, .$labels), ] %>%
      .$yi
    matched_vi <- interaction_results_vi %>%
      .[pmatch(labels, .$labels), ] %>%
      .$vi

    resultsB$yi[!resultsB$bold] <- matched_yi
    resultsB$vi[!resultsB$bold] <- matched_vi

    resultsB[, c("ci.lb", "ci.ub")] <-
      summary(escalc(yi = resultsB$yi, vi = resultsB$vi),
              level = 100 * prob)[, c("ci.lb", "ci.ub")]

    rows_to_update <- which(row_labels$subheading %in% model_names, c("yi", "ci.lb", "ci.ub"))
    include_indices <- seq(rows_to_update[1], rows_to_update[k], num_subgroups) + reference_group - 1

    resultsB[include_indices, c("yi", "ci.lb", "ci.ub")] <- interaction_intervals[seq(2, 2 * length(include_indices), 2), ]
    resultsB[resultsB$ci.ub > 195 & !is.na(resultsB$ci.ub), ] <- NA
    resultsB[which(resultsB$vi == 0), c("yi", "ci.lb", "ci.ub")] <- NA
    resultsB[is.na(resultsB$ci.lb), c("yi", "ci.lb", "ci.ub")] <- NA
    resultsB$colour[!is.na(resultsB$yi)] <- "black"
    resultsB$colour[resultsB$bold] <- "chartreuse4"

    if (missing(xlab)) {
      xlab <- c("", "")
    }
    if (missing(col_right_heading)) {
      col_right_heading <- list(paste0("Subgroup \n(", 100 * prob, "% interval)"), paste0("Interaction \n(", 100 * prob, "% interval)"))
    }

    row_labels$heading <- row_labels$subheading
    row_labels$subheading <- NA_character_

    resultsA$shape <- 22
    resultsB$shape <- 21
    resultsA$shape[resultsA$bold] <- 18
    resultsB$shape[resultsB$bold] <- 18

    col.left <- NULL
    left.heading <- ""

    if (all(treat.total != "")) {
      resultsA$treat <- c(paste0(treat.events, "/", treat.total),
                          rep("", nrow(resultsA) - length(treat.total)))
      resultsA$cont <- c(paste0(cont.events, "/", cont.total),
                         rep("", nrow(resultsA) - length(treat.total)))
      left.heading_list1 <- list(c("Control \n (n/N)", "Treatment \n (n/N)"), c("", ""))
    }
    if (!is.null(resultsA$treat)) {
      resultsA$treat[resultsA$treat == "0/0"] <- ""
      resultsA$cont[resultsA$cont == "0/0"] <- ""
    }

    if (all(total != "")) {
      resultsA$prev <- c(paste0("       ", subgroup.total, " (", round(subgroup.total / total * 100, digits = 0), "%)"),
                         rep("", nrow(resultsA) - length(total)))
      left.heading_list2 <- list(paste0("Prevalence \n (n/N)"), "")
    }

    if (all(treat.total != "") & !all(total != "")) {
      col.left = c("cont", "treat")
      left.heading <- left.heading_list1
    } else if (!all(treat.total != "") & all(total != "")) {
      col.left = c("prev")
      left.heading <- left.heading_list2
    } else if (all(treat.total != "") & all(total != "")) {
      col.left = c("cont", "treat", "prev")
      left.heading <- c(left.heading_list1[[1]], left.heading_list2)
    } else {
      col.left <- NULL
      left.heading <- ""
    }

    resultsB$treat <- resultsB$cont <- resultsB$prev <- ""

    idxA <- which(resultsA[, "ci.ub"] - resultsA[, "ci.lb"] > 39)
    idxB <- which(resultsB[, "ci.ub"] - resultsB[, "ci.lb"] > 39)
    resultsA[idxA, c("yi", "ci.lb", "ci.ub")] <- c(NA, NA, NA)
    resultsB[idxB, c("yi", "ci.lb", "ci.ub")] <- c(NA, NA, NA)

    if (!is.null(summary_fill)) {
      # Per-row fill/colour already set (e.g. prevalence forest plot); do not overwrite
      resultsB$fill <- resultsA$fill
    } else {
      if (is.list(fill)) {
        if (length(fill) != 2) stop("`fill` must be a list of length 2 when provided as a list")
        fill_list <- fill
      } else if (is.null(fill)) {
        fill_list <- list("white", "white")
      } else {
        fill_list <- rep(fill, 2)
      }
      resultsA$fill <- fill_list[[1]]
      resultsB$fill <- fill_list[[2]]
    }

    if (!missing(grouping)) {
      row_labels$subheading <- grouping
    }
    bold <- is.na(row_labels$heading)

    if (estimates_only) {
      idx <- resultsB$bold
      idx[is.na(idx)] <- FALSE
      row_labels <- row_labels[idx, ]
      resultsA <- resultsA[idx, ]
      resultsB <- resultsB[idx, ]
      bold <- bold[resultsB$bold]
    }

    # Option to build only one panel (e.g. which_panel = 1 for prevalence-only forest plot)
    if (!is.null(which_panel) && which_panel %in% c(1L, 2L)) {
      panels_to_use <- list(resultsA, resultsB)[which_panel]
      col_right_heading <- list(col_right_heading[[which_panel]])
      if (!missing(xlab) && length(xlab) >= which_panel) xlab <- xlab[which_panel]
      xlim <- if (is.list(xlim) && length(xlim) >= which_panel) xlim[[which_panel]] else xlim
      xticks <- if (is.list(xticks) && length(xticks) >= which_panel) xticks[[which_panel]] else xticks
    } else {
      panels_to_use <- list(resultsA, resultsB)
    }

    forestplot <-
      ckbplotr::forest_plot(
        base_size = base_size,
        panels = panels_to_use,
        col.key = "variable", col.estimate = "yi",
        col.lci = "ci.lb", col.uci = "ci.ub",
        right.heading = col_right_heading,
        heading.space = 1,
        row.labels = row_labels,
        row.labels.space = row.labels.space,
        xlab = xlab,
        showcode = FALSE,
        digits = digits,
        fill = "fill",
        shape = "shape",
        colour = "colour",
        col.left = col.left,
        left.heading = left.heading,
        bold.labels = bold,
        pointsize = pointsize,
        envir = environment(),
        left.space = left.space,
        exponentiate = exponentiate,
        xlim = xlim,
        xticks = xticks,
        addaes = list(point = paste0("stroke = stroke/10*((2*1.96)/(log(uci_transformed)-log(lci_transformed)))^(1/2)")), ...
      )

    attr(forestplot, "panels") <- list(resultsA, resultsB)
    return(forestplot)
  }
