CAMS submission cache (optional)
================================

Companion to https://github.com/rvpanaro/CAMS — not required for GitHub clone.

Contents
--------
- fit_*.rds          : brms / metafor fits from the submitted manuscript run
- tables/table3.csv  : reference output (Bayesian EDSS IF table)
- tables/table4.csv  : reference output (frequentist age + EDSS table)

Install
-------
1. Clone or download the CAMS GitHub repository.
2. Unzip this archive into the CAMS folder (same level as figure1.R).
   Overwrite nothing except adding missing fit_*.rds and tables/*.csv.

Quick check
-----------
  Rscript table3.R
  Rscript table4.R
  diff tables/table3.csv  <(unzip -p ... table3.csv)   # or compare to bundled reference

Figure scripts (figure1.R … figure7_*.R) use the same fit_*.rds and skip long refits.

Generated: manuscript submission cache bundle.
