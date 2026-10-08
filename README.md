# Serum 25-hydroxyvitamin D and bowel movement phenotypes among US adults

Analysis code and derived data for the manuscript:

> **Serum 25-hydroxyvitamin D and bowel movement phenotypes among US adults: a cross-sectional study of NHANES 2005–2010** (manuscript in preparation; target journal: *Scientific Reports*)

## Study summary

- **Design:** cross-sectional analysis of three pooled cycles of the National Health and Nutrition Examination Survey (NHANES 2005–2006, 2007–2008, 2009–2010) — the only cycles that fielded the Bowel Health Questionnaire (BHQ) with the Bristol Stool Form Scale (BSFS).
- **Analytic sample:** 13,266 adults aged ≥20 years (representing approximately 175.8 million US adults).
- **Exposure:** serum 25-hydroxyvitamin D [25(OH)D], modeled continuously, in quartiles, and with a 50 nmol/L clinical threshold; values harmonized across cycles following NCHS analytic guidance on the 2005–2006 assay transition.
- **Outcomes:** bowel movement phenotypes derived from the BSFS — constipation-type (types 1–2), normal (3–5), diarrhea-type (6–7).
- **Analysis:** survey-weighted multivariable logistic regression under NCHS weighting guidance (MEC examination weights divided across pooled cycles, stratification and PSUs); effect modification by sex, age group and BMI category; restricted cubic splines for dose–response; BMI mediation analysis; multiple sensitivity analyses.
- **Supplementary analysis (Part B):** cholecystectomy and hepatic biochemistry in the NHANES 2017–2020 pre-pandemic cycle (n = 8,377; 914 with a history of cholecystectomy).

## Repository structure

```
code/                     Analysis pipeline (R) — run from the repository root
  download_d1_data.r      downloads all required NHANES XPT files from CDC and verifies integrity
  w1_merge_d1.r           merges/cleans cycles, applies 25(OH)D harmonization, builds A_pool / B_pool
  w2_analysis_d1.r        main models (OR tables), RCS dose–response, subgroup forest (Fig2 / Fig3)
  w2b_analysis_d1.r       Part B supplementary analysis (Fig S1)
  w3_mediation_d1.r       BMI and dietary-fat mediation analysis (ACME/ADE); requires the DR1TOT XPT files, which are fetched by the download script
  w2_table1.r             weighted Table 1
  w2_full_tables.r        full model-output tables
  w2_prev.r               weighted prevalence estimates
  w2_pop_ci.r             weighted population totals / prevalence CIs
  w2_normal_ci.r          CI helper (normal-type prevalence)
  make_fig1.r             Fig 1 — participant flow diagram
  check_cols.r / check_pred.r / check_med.r                       QC helper scripts
nhanes_data/              derived pooled analysis datasets (A_pool.RDS, B_pool.RDS)
                          + merge verification summary (W1_verify_summary.csv)
                          raw XPT files are fetched by the download script (not committed)
W2_results/               result tables (CSV), spline point objects (RDS) and figures (PNG, 300 dpi)
```

## Reproduction

1. Install R (≥ 4.3) with the required packages:
   `foreign, dplyr, survey, splines, ggplot2, mediation`
2. From the repository root:

```r
source("code/download_d1_data.r")   # downloads NHANES XPT files into nhanes_data/ and verifies them
source("code/w1_merge_d1.r")        # builds A_pool.RDS / B_pool.RDS (pre-built copies included)
source("code/w2_analysis_d1.r")     # main analysis + Fig2 / Fig3 + result tables
source("code/w2b_analysis_d1.r")    # Part B supplementary analysis + Fig S1
source("code/w3_mediation_d1.r")    # BMI + dietary-fat mediation (uses DR1TOT fetched above)
source("code/w2_table1.r")          # Table 1
source("code/make_fig1.r")          # Fig 1
```

All scripts assume the repository root as the working directory (paths such as `nhanes_data/` and `W2_results/` are relative). Pre-built `A_pool.RDS` / `B_pool.RDS` are included so the analysis steps can be run without re-downloading the raw files.

## Data sources

All data are publicly available from the CDC National Center for Health Statistics:

- NHANES 2005–2006 / 2007–2008 / 2009–2010 — https://wwwn.cdc.gov/nchs/nhanes/
- NHANES 2017–2020 pre-pandemic (Part B) — https://wwwn.cdc.gov/nchs/nhanes/

Serum 25(OH)D harmonization across the 2005–2010 cycles follows the NCHS analytic note on the 2005–2006 assay transition.

## License

MIT License — see `LICENSE`.

## Citation

Please cite the manuscript once published. A code archive DOI will be added upon publication.
