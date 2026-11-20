# Sleep-duration trajectories and depressive symptom burden in CHARLS

Analysis scripts and software environment records for *Sleep-duration trajectories and depressive symptom burden: complementary trajectory and landmark analyses in CHARLS*.

The full-cohort analysis describes co-occurrence of five-wave sleep and symptom trajectories (2011-2020). The landmark analysis relates sleep trajectories during 2011-2015 to CES-D-10 scores of at least 10 at both the 2018 and 2020 visits. These are explanatory association analyses. The analyses were not preregistered; the paired history comparison was post hoc and exploratory.

## Files

- `01_cohort_derivation.R`: linkage, baseline covariates, sleep and CES-D scoring, vital status and eligibility.
- `02_trajectory_models.R`: full-cohort and landmark latent class mixed models.
- `03_landmark_main_and_mice.R`: covariate imputation, Models 1-3 and landmark L0-L5, clustered inference and standardization.
- `04_sensitivity_analyses.R`: restrictions, item deletion, paired history contrasts, continuous sleep, weighting, pain/ADL adjustment and nine-item history coverage (Table S19).
- `05_descriptive_outputs.R`: descriptive summaries, missingness, classification diagnostics and class-correspondence estimates.
- `06_index_training_sensitivity.R`: the reported 2015-survival-only training restriction and classification agreement (Table S13D).
- `07_exclude_long_sleep.R`: saved-MI checks and the post hoc >8-hour exclusion sensitivity (Table S18); no trajectory refitting or re-imputation.
- `R/`: shared analysis functions (`cluster_variance.R`, `mi_helpers.R`, `standardize_risk.R`).
- `sessionInfo.txt` and `environment/`: software information from the formal analyses.

No participant-level data, completed imputations, fitted participant-level objects, credentials or internal working records are distributed.

## Data access

Obtain CHARLS data independently through the [CHARLS website](https://charls.pku.edu.cn/) and [data portal](https://charls.charlsdata.com/), subject to registration and the data-use agreement. Do not upload these data or locally generated RDS files.

Set `CHARLS_RAW_ROOT` to a private directory with this layout:

```text
2011/demographic_background.dta
2011/household_roster.dta
2011/health_status_and_functioning.dta
2011/biomarkers.dta
2011/Blood_20140429.dta
2011/weight.dta
2013/Health_Status_and_Functioning.dta
2013/Exit_Interview.dta
2013/Verbal_Autopsy.dta
2015/Health_Status_and_Functioning.dta
2015/Sample_Infor.dta
2018/Health_Status_and_Functioning.dta
2018/Cognition.dta
2018/Sample_Infor.dta
2020/Health_Status_and_Functioning.dta
2020/Sample_Infor.dta
2020/Exit_Module.dta
```

If official status modules are stored separately, `CHARLS_STATUS_ROOT` can point to a directory with the same year/subfile layout for `Sample_Infor.dta` and `Exit_Module.dta`. The two 2013 death files remain under `CHARLS_RAW_ROOT`.

The baseline release is *2011 CHARLS Wave1 (Baseline)*, with questionnaire modules dated 12 March 2013; the blood file is `Blood_20140429.dta`.

## Environment and execution

Formal analyses used R 4.5.1 on Linux, mice 3.19.0, lcmm 2.2.2, sandwich 3.1-3 and survey 4.5. Session information records the other versions. Scripts also use haven and dplyr (alongside R recommended packages MASS and nnet, and base parallel and splines). Install these dependencies before execution. The checked imputation methods require mice 3.19.0. The long-sleep module checks R 4.5.1, mice 3.19.0 and sandwich 3.1-3 before accessing saved study inputs.

Run from the repository root in separate R sessions:

```bash
export CHARLS_RAW_ROOT=/path/to/private/charls
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript --vanilla 01_cohort_derivation.R
Rscript --vanilla 02_trajectory_models.R
Rscript --vanilla 03_landmark_main_and_mice.R
Rscript --vanilla 04_sensitivity_analyses.R
Rscript --vanilla 05_descriptive_outputs.R
Rscript --vanilla 06_index_training_sensitivity.R
Rscript --vanilla 07_exclude_long_sleep.R run
```

On Windows, set the same variables in PowerShell, for example `$env:CHARLS_RAW_ROOT = 'data/raw_charls'`, then run the Rscript commands. Local data/model outputs go to `output/`, excluded by `.gitignore`. Sourcing an analysis script loads functions without starting that analysis.

## Analysis settings

CES-D-10 requires at least eight observed items, prorated to ten. CES-D-9 removes restless sleep while retaining the original reverse-scoring identities and requires eight of nine retained items. The complete nine-item subset requires all nine retained items at all five waves, not the excluded sleep item. Sleep, symptom scores, class assignments and outcomes are not imputed.

Trajectory models have class-specific quadratic time terms, random intercepts and linear slopes, and shared unstructured random-effect covariance and residual variance. Searches use seeds 2025, 2026 and 2027, 32 starts, 30 short iterations and up to 500 full iterations with a permitted 500-iteration continuation. Within each class count, the lowest-BIC converged candidate is considered first, then class-size (at least 5%) and average posterior-probability (at least 0.70) criteria are applied. The full-cohort sleep, full-cohort ten-item symptom and landmark sleep searches evaluate one to four classes. The paired ten-/nine-item symptom searches and 2015-survival-only training sensitivity evaluate one and two classes. The PSOCK worker settings are 15, 5, 24, 14 and 14 for the five stage-02 domains; stage 06 uses 24.

Covariate MI uses `m = 100`, `maxit = 20`, `seed = 2025` and `RNGkind("Mersenne-Twister", "Inversion", "Rejection")`. Exposure class, outcome and appropriate symptom-history variables are predictors, not targets. BMI uses predictive mean matching; binary and unordered categorical variables use logistic and multinomial imputation.

Regressions use community-clustered HC1 covariance. Rubin pooling uses Barnard-Rubin degrees of freedom with complete-data degrees of freedom `G - 1`. Standardized differences use coefficient delta-method variances, conditional on fitted classes and empirical target distributions. Paired history comparisons use common participants and imputations, with stacked community-level influence functions retaining cross-model covariance.

The survey sensitivity first builds the design from the complete positive-weight baseline file, then subsets the landmark domain. The survivor-response analysis targets survivors with observed three-wave symptom histories; deaths and unresolved or conflicting survival evidence are not treated as ordinary nonresponse.

## Citation and reuse

Please cite the accompanying manuscript and archived release DOI when reusing these scripts. The analysis code is released under the MIT License (see `LICENSE`). This licence applies solely to the software scripts and does not extend to CHARLS participant data, which remain governed by the CHARLS Data Use Agreement.
