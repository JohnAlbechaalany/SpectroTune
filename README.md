# SpectroTune

SpectroTune is an R package for spectral pretreatment, model benchmarking with nested cross‑validation, plotting, and exporting final models. It supports regression (numeric targets) and classification (factor targets).

## Installation

```r
# Ensure devtools/roxygen2 are available
install.packages(c("devtools","roxygen2"))

# From local folder (replace path if needed)
devtools::document("spectrotune")
devtools::install("spectrotune", dependencies = TRUE)

library(spectrotune)
```

If you are on Windows and see permission DLL warnings, restart R, delete any `00LOCK` folders under your R library, and reinstall. Rtools is recommended for building from source: `https://cran.rstudio.com/bin/windows/Rtools/`.

## Key features

- Pretreatments: RAW, ABS, SNV, SG0/SG1/SG2, MSC (with tunable Savitzky–Golay params)
- Range control: exclude noise bands, keep wavelength ranges
- Models (regression): PLSR, SVR (RBF), GLMNET (Elastic Net), RIDGE, Kernel Ridge (RBF), Random Forest, XGBoost on PCA scores
- Models (classification): SVM (RBF), GLMNET (logistic/multinomial), Random Forest, XGBoost on PCA scores
- Nested CV (outer/inner) with flexible sampling:
  - Regression (numeric y): random or stratified by y bins (quantiles)
  - Independent covariates: random or stratified by provided factors (e.g., breed, sex)
- Metrics
  - Regression: SEC, RSQcal, SEP, RSQval, Bias
  - Classification: Accuracy, Balanced Accuracy, Macro F1
- Plot best pipeline (scatter or confusion matrix)
- Fit and save final model for deployment

## API overview

- `st_preprocess(data, target, preprocess, preprocess_params, noise_range, wl_keep)`
  - Returns a list with X (per-method matrices), y, wl
- `st_model(X_list, y, algos, outer_k, inner_k, outer_rep, seed, dependent_sampling, dependent_bins, independent_sampling, stratify, verbose)`
  - Runs nested CV and returns `list(results, summary, best)`
- `st_run(data, target, preprocess, preprocess_params, noise_range, wl_keep, algos, outer_k, inner_k, outer_rep, seed, verbose)`
  - Convenience wrapper: preprocess + model
- `st_plot_best(fit, task = "auto", select_by = NULL, preprocess_filter = NULL, algo_filter = NULL, return_data = FALSE)`
  - Plots best pipeline (regression scatter or classification confusion matrix)
- `st_fit_final(data, target, preprocess, algo, preprocess_params, inner_k, seed)`
  - Refits a chosen pipeline on full data with inner CV for hyperparameters; returns a serializable model
- `st_save_final(final_model, file)`
  - Saves the final model to an RDS file

## Quick start

### Regression (numeric target)
```r
# df: data.frame with spectral columns named by wavelength (e.g., "400","405",...) and a numeric target
fit_reg <- st_run(
  data  = df,
  target = "Age",
  preprocess = c("ABS","SG1"),
  preprocess_params = list(SG1 = list(w = 17, p = 2, m = 1)),
  algos = c("PLSR","SVR","GLMNET","RIDGE","KRR_RBF","RF","XGB_PCA"),
  outer_k = 5, inner_k = 5, outer_rep = 2, seed = 42,
  dependent_sampling = "stratified_bin", dependent_bins = 5,
  independent_sampling = "random",
  verbose = TRUE
)

fit_reg$summary  # columns: SEC, RSQcal, SEP, RSQval, Bias
p <- st_plot_best(fit_reg)  # truth vs pred scatter
print(p)
```

### Classification (factor target)
```r
df$Class <- as.factor(df$Class)

fit_cls <- st_run(
  data  = df,
  target = "Class",
  preprocess = c("ABS","SG1"),
  preprocess_params = list(SG1 = list(w = 17, p = 2, m = 1)),
  algos = c("SVM","GLMNET","RF","XGB_PCA"),
  outer_k = 5, inner_k = 5, outer_rep = 2, seed = 42,
  independent_sampling = "stratified", stratify = list(breed = df$breed, sex = df$sex),
  verbose = TRUE
)

fit_cls$summary  # columns: Accuracy, Balanced_Accuracy, F1_macro
p <- st_plot_best(fit_cls)  # confusion matrix heatmap
print(p)
```

## Save a final model

```r
final <- st_fit_final(
  data = df,
  target = "Age",
  preprocess = "SG1",
  algo = "SVR",
  preprocess_params = list(w = 17, p = 2, m = 1),
  inner_k = 5, seed = 42
)

st_save_final(final, "spectrotune_final_model.rds")
```

## Sampling controls (simple mental model)

- Dependent (numeric y):
  - `dependent_sampling = "random"` → plain k‑fold
  - `dependent_sampling = "stratified_bin"`, `dependent_bins = 5` → y quantile bins, balanced across folds
- Independent covariates:
  - `independent_sampling = "random"` → ignore covariates for folds
  - `independent_sampling = "stratified"`, `stratify = list(breed=df$breed, sex=df$sex)` → balance these levels across folds

## Notes

- Large grids and many algorithms can be expensive; start with smaller sets, then expand.
- For tree/boosted models on spectra, PCA features (`XGB_PCA`) often generalize better.
- Set `seed` for reproducibility.
