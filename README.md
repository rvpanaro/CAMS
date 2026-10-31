# CAMS

Source code and materials for *Consistent Bayesian meta-analysis on subgroup specific effects and interactions*  
Renato Panaro, Christian Röver, Tim Friede

Preprint: https://arxiv.org/abs/2512.18785  
DOI: https://doi.org/10.48550/arXiv.2512.18785

## Project overview

This repository contains the R scripts and CSV inputs used to reproduce the manuscript figures and tables. Set the working directory to **this folder** (`CAMS`). Model fits (`.rds`), figure PDFs, and exported table CSVs are generated locally and are **not** stored on GitHub.

Use **this folder** as the reproducibility bundle: zip for upload excluding `.git`, `.Rhistory`, and gitignored caches (`.rds`, `figures/*.pdf`, `tables/*.csv`).

## Repository structure

- `figure*.R`: scripts for manuscript Figures 1–7 (PDFs written to `figures/`).
- `table3.R`: Bayesian comparison of IF methods (EDSS) → `tables/table3.csv`.
- `table4.R`: Frequentist Age and EDSS analyses → `tables/table4.csv` (no Stata if `wt_stata_*_results.csv` are present).
- `interaction_forest_plot.R`: shared plotting helper.
- `age_data.csv`, `edss_data.csv`: long-format inputs.
- `wt_stata_age_results.csv`, `wt_stata_edss_results.csv`: cached within-trial (metafloat) results for `table4.R`.

## R packages

Scripts use **R ≥ 4.1** (recommended). Bayesian fits (`brms`) need a **Stan** backend and a working **C++** toolchain (Xcode Command Line Tools on macOS, `r-base-dev` / `g++` on Linux).

### CRAN

Install these from CRAN before running any script:

| Package | Used for |
|---------|----------|
| `brms` | Bayesian CAMS and figure fits |
| `metafor` | Frequentist CAMS (`table4.R`, forest plots) |
| `RBesT` | Figure 7 (`gMAP`, prevalence) |
| `dplyr`, `ggplot2`, `scales`, `patchwork`, `gtable`, `grid`, `magrittr` | Data wrangling and figures |
| `ggtext`, `gghalves` | Figure labels / half-violin panels |
| `xtable` | LaTeX table output in `table3.R` |

One-shot install:

```r
install.packages(c(
  "brms", "metafor", "RBesT", "dplyr", "ggplot2", "scales", "patchwork",
  "gtable", "magrittr", "ggtext", "gghalves", "xtable"
))
```

### Not on CRAN (required for Figures 1, 5, and `interaction_forest_plot.R`)

**`ckbplotr`** — install from the [neilstats R-universe](https://neilstats.r-universe.dev/ckbplotr) (or [GitHub `neilstats/ckbplotr`](https://github.com/neilstats/ckbplotr)):

```r
install.packages(
  "ckbplotr",
  repos = c("https://neilstats.r-universe.dev", getOption("repos"))
)
# alternative:
# remotes::install_github("neilstats/ckbplotr")
```

### Stan backend for `brms` (required for `table3.R` and Bayesian figures)

**Recommended:** [`cmdstanr`](https://mc-stan.org/cmdstanr/) + CmdStan (not on CRAN):

```r
install.packages(
  "cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
)
cmdstanr::install_cmdstan()
```

Point `brms` at CmdStan if needed: `options(brms.backend = "cmdstan")`.

**Alternative:** CRAN package **`rstan`** (older toolchain; `brms` can use it if CmdStan is unavailable).

## Reproduction

1. Set R’s working directory to **`CAMS`**.
2. Install packages as above (including **`ckbplotr`** and **Stan**).
3. Run e.g. `Rscript figure1.R`, `Rscript table3.R`, `Rscript table4.R`.

**Settings:** `set.seed(1)`, `seed = 1` in `brm()`, `chains = 4`, `iter = 20000`, `warmup = 6000`, `control = list(adapt_delta = 0.999, max_treedepth = 15)`.

**Note:** `table4.R` does not need Stata if `wt_stata_age_results.csv` and `wt_stata_edss_results.csv` are present. First runs of Bayesian scripts can take a long time while `.rds` caches are created.

## Manuscript mapping

| Output | Script(s) | Local output path |
|--------|-----------|-------------------|
| Figure 1 (age) | `figure1.R` | `figures/figure1.pdf` |
| Figure 2 (EDSS) | `figure2.R` | `figures/figure2.pdf` |
| Figure 3 | `figure3_lhs.R`, `figure3_rhs.R` | `figures/figure3_*.pdf` |
| Figure 4 | `figure4_lhs.R`, `figure4_rhs.R` | `figures/figure4_*.pdf` |
| Figure 5 | `figure5_top.R`, `figure5_bottom.R` | `figures/figure5_*.pdf` |
| Figure 6 | `figure6_lhs.R`, `figure6_center.R`, `figure6_rhs.R` | `figures/figure6_*.pdf` |
| Figure 7 | `figure7_lhs.R`, `figure7_rhs.R` | `figures/figure7_*.pdf` |
| Table 3 (Bayesian IF, EDSS) | `table3.R` | `tables/table3.csv` |
| Table 4 (Frequentist age + EDSS) | `table4.R` | `tables/table4.csv` |

Optional local caches (gitignored): `fit_*.rds` produced on first run of figure or table scripts.
