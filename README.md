# CAMS

Source code and materials for *Consistent Bayesian meta-analysis on subgroup specific effects and interactions*  
Renato Panaro, Christian Röver, Tim Friede

Preprint: https://arxiv.org/abs/2512.18785  
DOI: https://doi.org/10.48550/arXiv.2512.18785

## Project overview
This repository contains the R scripts and CSV inputs used to reproduce the manuscript figures and Table 2 for the paper.
The scripts are designed to be run with the working directory set to **this folder** (`CAMS`); they read the CSV files here and generate cached model objects (`.rds`) on first run.

## Repository structure
The repository mirrors the original materials layout and keeps files at the project root:
- `figure*.R`: R scripts generating the manuscript figures (interaction forest plots, subgroup effects, conditional effects, and prevalence plots).
- `table2.R`: script generating `table2.csv`, the dataset underlying Table 2.
- `interaction_forest_plot.R`: shared plotting helper used by multiple figure scripts.
- `age_data.csv`, `edss_data.csv`: long-format input data for Age and EDSS.
- `table2.csv`: output dataset for Table 2.

## Zenodo / self-contained reproduction
Use **this folder** (`CAMS`) as the reproducibility bundle: zip it for upload (exclude `.git`, `.Rhistory`, and large optional `.rds` caches if you want a smaller archive).

1. Set the R working directory to **`CAMS`** (the folder containing the scripts and CSVs).
2. Run or source the desired `figure*.R` or `table2.R`. Figures can be written with `ggsave()` where enabled; `brms` fits are cached as `.rds` on first run (delete a `.rds` to force refit).

**Reproducibility settings** (see scripts): `set.seed(1)`, `seed = 1` in `brm()`, `chains = 4`, `iter = 20000`, `warmup = 6000`, `control = list(adapt_delta = 0.999, max_treedepth = 15)`.

**Manuscript figure → script (matches RSM paper Figure 1–7 order)**

| Figure | Script(s) | Output PDF(s) |
|--------|------------|---------------|
| 1 — Interaction forest (age) | `figure1.R` | `figure1.pdf` |
| 2 — Interaction forest (EDSS) | `figure2.R` | `figure2.pdf` |
| 3 — Conditional effects / bubble plots | `figure3_lhs.R` (age), `figure3_rhs.R` (EDSS) | `figure3_lhs.pdf`, `figure3_rhs.pdf` |
| 4 — Optimal IF (CrI width) + conditional RR (EDSS) | `figure4_lhs.R`, `figure4_rhs.R` | `figure4_lhs.pdf`, `figure4_rhs.pdf` |
| 5 — Three-panel forest (+ overall), age / EDSS | `figure5_top.R`, `figure5_bottom.R` | `figure5_top.pdf`, `figure5_bottom.pdf` |
| 6 — Posterior panels (inference / ecological / trace) | `figure6_lhs.R`, `figure6_center.R`, `figure6_rhs.R` | `figure6_lhs.pdf`, `figure6_center.pdf`, `figure6_rhs.pdf` |
| 7 — Prevalence (RBesT) | `figure7_lhs.R`, `figure7_rhs.R` | `figure7_lhs.pdf`, `figure7_rhs.pdf` |
| Table 2 | `table2.R` | `table2.csv` |

**Optional cache files** (speed reruns; omit from minimal uploads if desired): `fit_age.rds`, `fit_edss.rds`, `fit_subgroup.rds`, `fit_interaction.rds`, `fit_subgroup_edss.rds`, `fit_interaction_edss.rds`, `fit_overall_age.rds`, `fit_overall_edss.rds`, and other `fit_*.rds` produced by the figure scripts.
