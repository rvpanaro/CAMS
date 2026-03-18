# CAMS

Source code and materials for *Consistent Bayesian meta-analysis on subgroup specific effects and interactions*  
Renato Panaro, Christian Röver, Tim Friede

Preprint: https://arxiv.org/abs/2512.18785  
DOI: https://doi.org/10.48550/arXiv.2512.18785

## Project overview
This repository contains the R scripts and CSV inputs used to reproduce the manuscript figures and Table 2 for the paper.
The scripts are designed to be run from the repository root; they read the CSV files in this directory and generate cached model objects (`.rds`) on first run.

## Repository structure
The repository mirrors the original materials layout and keeps files at the project root:
- `figure*.R`: R scripts generating the manuscript figures (interaction forest plots, subgroup effects, conditional effects, and prevalence plots).
- `table2.R`: script generating `table2.csv`, the dataset underlying Table 2.
- `interaction_forest_plot.R`: shared plotting helper used by multiple figure scripts.
- `age_data.csv`, `edss_data.csv`: long-format input data for Age and EDSS.
- `table2.csv`: output dataset for Table 2.
