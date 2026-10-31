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

## Reproduction

1. Set R’s working directory to **`CAMS`**.
2. Install R packages used in the scripts (`brms`, `metafor`, `RBesT`, etc.).
3. Run e.g. `Rscript figure1.R`, `Rscript table3.R`, `Rscript table4.R`.

**Settings:** `set.seed(1)`, `seed = 1` in `brm()`, `chains = 4`, `iter = 20000`, `warmup = 6000`, `control = list(adapt_delta = 0.999, max_treedepth = 15)`.

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
