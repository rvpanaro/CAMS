#
# Figure 7 RHS: Approximate predictive density (MAP samples) and Beta(5, 21) prior.
# Self-contained; uses RBesT gMAP + automixfit. Requires: RBesT, scales.
#

library("RBesT")
library(ggplot2)
library(scales)

# Same data as figure7_lhs (prevalence by study)
labels <- c("DEFINE", "FREEDOMS", "CONFIRM-B", "CONFIRM-C", "ADVANCE", "AFFIRM", "TEMSO", "COP1 PILOT")
B.n <- c(73, 62, 60, 60, 89, 79, 164, 12)
A.n <- c(336, 363, 299, 290, 423, 548, 557, 13)

set.seed(1)
map01 <- gMAP(cbind(B.n, A.n) ~ 1 | labels,
              family = binomial,
              tau.dist = "HalfNormal",
              tau.prior = 1.0,
              iter = 20000,
              beta.prior = cbind(0.0, 1.814))

map <- RBesT::automixfit(map01)
dens <- plot(map)[[2]] + labs(title = "Approximate predictive density")

map_key   <- 'plain("MAP samples")'
prior_key <- 'pi[j] %~% B(5, 21)'

dens1 <- dens +
  geom_line(
    data = data.frame(x = 0, y = 0),
    aes(x, y, colour = map_key),
    inherit.aes = FALSE,
    linewidth = 1,
    alpha = 0
  ) +
  geom_line(
    data = data.frame(x = 0, y = 0),
    aes(x, y, colour = prior_key),
    inherit.aes = FALSE,
    linewidth = 1,
    alpha = 0
  ) +
  scale_colour_manual(
    values = c(
      `plain("MAP samples")` = "#2C7FB8",
      `pi[j] %~% B(5, 21)`   = "black"
    ),
    breaks = c(map_key, prior_key),
    labels = label_parse()
  ) +
  guides(colour = guide_legend(
    title = NULL,
    override.aes = list(alpha = 1, linewidth = 1)
  )) +
  theme_classic(base_size = 20) +
  theme(
    legend.position      = c(0.6, 0.95),
    legend.justification = c("left", "top"),
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

print(dens1)
ggsave("figure7_rhs.pdf", plot = dens1, width = 6, height = 6, units = "in", device = grDevices::pdf)
