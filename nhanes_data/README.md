# Derived analysis datasets

- `A_pool.RDS` — pooled, cleaned analysis dataset for the main (Part A) analysis: NHANES 2005–2010 adults aged ≥20 with serum 25(OH)D, Bristol stool type and complete covariates (analysis n = 13,266). Built by `code/w1_merge_d1.r`.
- `B_pool.RDS` — cleaned dataset for the Part B supplementary analysis (NHANES 2017–2020 pre-pandemic; cholecystectomy and hepatic biochemistry).
- `W1_verify_summary.csv` — merge/verification summary produced by `code/w1_merge_d1.r`.

Raw NHANES XPT files are **not** committed. Run `code/download_d1_data.r` from the repository root to download and verify them into this folder.
