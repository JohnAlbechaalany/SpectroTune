# High-level wrapper to choose preprocessing, QC, folds, target, etc.
spectrotune_benchmark <- function(data,
                                  target,
                                  preprocess = c('RAW','ABS','SNV','SG1','SG2','MSC','SG0'),
                                  algos = c('PLSR','SVR'),
                                  noise_range = NOISE_RANGE,
                                  wl_keep = NULL,
                                  qc = FALSE,
                                  outer_k = OUTER_K,
                                  outer_rep = OUTER_REP,
                                  inner_k = INNER_K,
                                  seed = 42,
                                  preprocess_params = list()) {
  stopifnot(is.data.frame(data))
  stopifnot(target %in% names(data))
  wl_cols <- get_spectral_cols(data)
  wl <- as.numeric(wl_cols)
  X_raw <- as.matrix(data[, wl_cols, drop = FALSE])
  if (!is.null(noise_range)) {
    exc <- exclude_noise_range(X_raw, wl, noise_range)
    X_raw <- exc$X
    wl    <- exc$wl
  }
  if (!is.null(wl_keep)) {
    kept <- keep_wavelength_ranges(X_raw, wl, wl_keep)
    X_raw <- kept$X
    wl    <- kept$wl
  }
  # Define preprocessors
  prep_map <- list(RAW = prep_raw, ABS = prep_abs, SNV = prep_snv_abs, SG1 = prep_sg1_abs, SG2 = prep_sg2_abs,
                   MSC = prep_msc_abs, SG0 = prep_sg0_abs)
  preprocess <- match.arg(preprocess, choices = names(prep_map), several.ok = TRUE)
  algos <- match.arg(algos, choices = c('PLSR','SVR'), several.ok = TRUE)
  # No QC by default; if enabled, use PCA+Mahalanobis
  keep_idx <- seq_len(nrow(X_raw))
  if (isTRUE(qc)) {
    X_abs <- prep_abs(X_raw)
    qc_res <- qc_filter_mahal(X_abs)
    keep_idx <- qc_res$keep
  }
  df <- data[keep_idx, , drop = FALSE]
  X  <- X_raw[keep_idx, , drop = FALSE]
  Y  <- df[[target]]
  # Build pipelines
  # Build param-aware preprocessors
  build_preproc <- function(name) {
    fn <- prep_map[[name]]
    params <- preprocess_params[[name]]
    if (is.null(params)) return(function(X) fn(X))
    function(X) do.call(fn, c(list(X = X), params))
  }
  pipes <- expand_grid(prep_name = preprocess, algo = algos) %>%
    mutate(prep_fn = I(lapply(prep_name, build_preproc)))
  results <- vector('list', nrow(pipes))
  for (i in seq_len(nrow(pipes))) {
    cat(sprintf('\nRunning pipeline %d/%d: %s + %s\n', i, nrow(pipes), pipes$prep_name[i], pipes$algo[i]))
    oof <- run_nested_cv(X, Y, preproc_fn = pipes$prep_fn[[i]], algo = pipes$algo[i],
                         outer_k = outer_k, outer_rep = outer_rep, inner_k = inner_k, seed = seed)
    results[[i]] <- oof %>% mutate(prep = pipes$prep_name[i], algo = pipes$algo[i])
  }
  res_df <- bind_rows(results) %>% left_join(tibble(idx = seq_along(Y), truth = Y), by = 'idx')
  summ <- res_df %>%
    group_by(prep, algo) %>%
    summarize(
      RSQ_val = rsq(truth, pred),
      MAE_val = Metrics::mae(truth, pred),
      RMSE_val = Metrics::rmse(truth, pred),
      INNER_R2_mean = mean(inner_R2, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(RSQ_val), MAE_val)
  list(results = res_df, summary = summ, best = summ %>% slice(1))
}
# =============================================================
# Spectral QC, Preprocessing, and Nested CV (PLSR & SVR)
# -------------------------------------------------------------
# Requirements: install.packages(c(
#   'readxl','dplyr','tidyr','tibble','stringr','purrr',
#   'ggplot2','scales','prospectr','pls','kernlab','Metrics','caret'
# ))
# If you prefer tidymodels, you can translate the modeling block accordingly.
# =============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(prospectr)   # SNV, savitzkyGolay
  library(pls)         # PLSR
  library(kernlab)     # SVR (ksvm)
  library(Metrics)     # mae, rmse, etc.
  library(caret)       # CV utilities
})

# (Removed legacy alias; using exported function prospectr::standardNormalVariate)

# -----------------------------
# User parameters
# -----------------------------
SAMPLE_MODE <- TRUE            # set TRUE to run with synthetic sample data
FILE_PATH   <- 'C:/Users/user/Desktop/df25.xlsx'     # adjust path as needed
SHEET_NAME  <- 'Overall'
TARGET_COL  <- 'PV'            # <-- put your target variable here (e.g., 'PV' or 'PM')
NOISE_RANGE <- c(700, 705)     # wavelengths to exclude [nm]
OUTER_K     <- 5               # outer CV folds
OUTER_REP   <- 2               # outer CV repeats
INNER_K     <- 5               # inner CV folds

# Hyperparameter grids
PLS_K_GRID  <- 4:18
SVR_C_GRID  <- c(0.5, 1, 2, 4)
SVR_SIGMA_GRID <- c(0.01, 0.1, 1, 10) # rbfdot sigma in ksvm

# -----------------------------
# Helpers
# -----------------------------
read_data <- function(path = FILE_PATH, sheet = SHEET_NAME) {
  df <- readxl::read_excel(path, sheet = sheet)
  as_tibble(df)
}

# Sample data generator: creates spectral columns named as wavelengths plus a target column
generate_sample_data <- function(num_samples = 80, num_wavelengths = 200, 
                                 wl_start = 400, wl_step = 5, target_col = TARGET_COL) {
  set.seed(1)
  wavelengths <- seq(from = wl_start, by = wl_step, length.out = num_wavelengths)
  X <- matrix(rnorm(num_samples * num_wavelengths, mean = 0.5, sd = 0.2),
              nrow = num_samples, ncol = num_wavelengths)
  colnames(X) <- as.character(wavelengths)
  # Create target as a smooth function of a few wavelengths plus noise
  wl_idx1 <- which.min(abs(wavelengths - 550))
  wl_idx2 <- which.min(abs(wavelengths - 900))
  wl_idx3 <- which.min(abs(wavelengths - 1200))
  signal <- 2*X[, wl_idx1] - 1.5*X[, wl_idx2] + 1.0*X[, wl_idx3]
  target <- scales::rescale(signal + rnorm(num_samples, sd = 0.3), to = c(0, 100))
  df <- as_tibble(X)
  df[[target_col]] <- target
  df
}

get_spectral_cols <- function(df) {
  cols <- colnames(df)
  num_cols <- cols[ suppressWarnings(!is.na(as.numeric(cols))) ]
  num_cols
}

exclude_noise_range <- function(X, wavelengths, noise_range = NOISE_RANGE) {
  keep <- !(wavelengths >= noise_range[1] & wavelengths <= noise_range[2])
  list(X = X[, keep, drop = FALSE], wl = wavelengths[keep])
}

# Keep-only wavelength range(s). wl_keep can be a numeric vector length-2 or a list of pairs
keep_wavelength_ranges <- function(X, wavelengths, wl_keep = NULL) {
  if (is.null(wl_keep)) return(list(X = X, wl = wavelengths))
  if (is.list(wl_keep)) {
    mask <- rep(FALSE, length(wavelengths))
    for (rng in wl_keep) {
      mask <- mask | (wavelengths >= rng[1] & wavelengths <= rng[2])
    }
  } else {
    mask <- (wavelengths >= wl_keep[1] & wavelengths <= wl_keep[2])
  }
  list(X = X[, mask, drop = FALSE], wl = wavelengths[mask])
}

# Preprocessing functions (expect matrix reflectance X in [n_samples, n_wl])
prep_raw <- function(X) X
prep_abs <- function(X) -log10(pmax(X, .Machine$double.eps))
prep_snv_abs <- function(X) prospectr::standardNormalVariate(prep_abs(X))
prep_sg1_abs <- function(X, w = 11, p = 2, m = 1) savitzkyGolay(prep_abs(X), m = m, p = p, w = w)
prep_sg2_abs <- function(X, w = 15, p = 2, m = 2) savitzkyGolay(prep_abs(X), m = m, p = p, w = w)

# Additional preprocessing examples
prep_msc_abs <- function(X) prospectr::msc(prep_abs(X))
prep_sg0_abs <- function(X, w = 11, p = 2, m = 0) savitzkyGolay(prep_abs(X), m = m, p = p, w = w)  # smoothing only

PREP_LIST <- list(
  RAW = prep_raw,
  ABS = prep_abs,
  SNV = prep_snv_abs,
  SG1 = prep_sg1_abs,
  SG2 = prep_sg2_abs,
  MSC = prep_msc_abs,
  SG0 = prep_sg0_abs
)

# PCA-based QC on spectra (after a baseline preprocess, e.g., ABS)
qc_filter_mahal <- function(X, var_explained = 0.95, alpha = 0.001) {
  # Center columns, PCA
  pca <- prcomp(scale(X, center = TRUE, scale = TRUE), center = FALSE, scale. = FALSE)
  # Determine number of PCs to reach var_explained
  cumvar <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  q <- which(cumvar >= var_explained)[1]
  if (is.na(q) || q < 2) q <- min(5, ncol(pca$x))
  Z <- pca$x[, 1:q, drop = FALSE]
  # Mahalanobis on PC scores
  mu <- colMeans(Z)
  S  <- cov(Z)
  d2 <- mahalanobis(Z, center = mu, cov = S)
  thr <- qchisq(1 - alpha, df = q)
  keep_idx <- which(d2 <= thr)
  list(keep = keep_idx, d2 = d2, thr = thr, q = q)
}

# Utility: fit + predict PLSR
fit_plsr <- function(X, y, ncomp) {
  df <- as.data.frame(X)
  df$y <- y
  model <- pls::plsr(y ~ ., data = df, ncomp = ncomp, validation = 'none', scale = TRUE)
  model
}

predict_plsr <- function(model, Xnew) as.numeric(predict(model, newdata = as.data.frame(Xnew), ncomp = model$ncomp))

# Utility: fit + predict SVR (RBF)
fit_svr <- function(X, y, C, sigma) {
  dat <- as.data.frame(X)
  dat$y <- y
  k <- rbfdot(sigma = sigma)
  model <- ksvm(as.matrix(dat[ , setdiff(colnames(dat), 'y')]), y,
                type = 'eps-svr', kernel = k, C = C, scaled = TRUE)
  model
}

predict_svr <- function(model, Xnew) as.numeric(predict(model, as.matrix(Xnew)))

# Metrics
rsq <- function(truth, estimate) {
  ss_res <- sum((truth - estimate)^2)
  ss_tot <- sum( (truth - mean(truth))^2 )
  1 - ss_res/ss_tot
}

# -----------------------------
# Nested CV runner for one (preprocess, algorithm)
# -----------------------------
run_nested_cv <- function(X, y, preproc_fn, algo = c('PLSR','SVR'),
                          outer_k = OUTER_K, outer_rep = OUTER_REP,
                          inner_k = INNER_K, seed = 42) {
  algo <- match.arg(algo)
  
  # Adjust outer folds to data size (caret requires 1 < k <= length(y))
  n_obs <- length(y)
  if (n_obs < 3) {
    stop(sprintf('Too few samples after QC: %d. Need at least 3 to run CV.', n_obs))
  }
  outer_k_eff <- min(OUTER_K, max(2, n_obs))
  
  # Outer CV index list with repeats
  set.seed(seed)
  outer_splits <- vector('list', outer_rep * outer_k_eff)
  idx <- 1
  for (r in seq_len(outer_rep)) {
    folds <- createFolds(y, k = outer_k_eff, returnTrain = TRUE)
    for (k in seq_len(outer_k_eff)) {
      outer_splits[[idx]] <- list(train = folds[[k]], test = setdiff(seq_along(y), folds[[k]]))
      idx <- idx + 1
    }
  }
  
  oof_pred <- rep(NA_real_, length(y))
  oof_fold <- rep(NA_integer_, length(y))
  oof_inner_r2 <- rep(NA_real_, length(y))
  chosen_hparams <- list()
  
  for (i in seq_along(outer_splits)) {
    tr <- outer_splits[[i]]$train
    te <- outer_splits[[i]]$test
    
    Xtr <- preproc_fn(X[tr, , drop = FALSE])
    Xte <- preproc_fn(X[te, , drop = FALSE])
    ytr <- y[tr]
    
    # Adjust inner folds to training size
    inner_k_eff <- min(inner_k, max(2, length(ytr)))
    inner_folds <- createFolds(ytr, k = inner_k_eff, returnTrain = TRUE)
    
    # Hyperparameter grid, bounded by data when needed
    if (algo == 'PLSR') {
      max_ncomp <- max(1, min(ncol(Xtr), length(ytr) - 1))
      Ks <- intersect(PLS_K_GRID, seq_len(max_ncomp))
      if (length(Ks) == 0) Ks <- max(1, min(PLS_K_GRID))
      grid <- tibble(model = 'PLSR', K = Ks)
    } else {
      grid <- expand_grid(model = 'SVR', C = SVR_C_GRID, sigma = SVR_SIGMA_GRID)
    }
    
    val_scores <- grid %>% mutate(R2 = NA_real_)
    
    for (g in seq_len(nrow(grid))) {
      preds <- rep(NA_real_, length(ytr))
      
      for (f in seq_along(inner_folds)) {
        tr_in <- inner_folds[[f]]
        va_in <- setdiff(seq_along(ytr), tr_in)
        
        if (algo == 'PLSR') {
          m <- fit_plsr(Xtr[tr_in, , drop = FALSE], ytr[tr_in], ncomp = grid$K[g])
          preds[va_in] <- predict_plsr(m, Xtr[va_in, , drop = FALSE])
        } else {
          m <- fit_svr(Xtr[tr_in, , drop = FALSE], ytr[tr_in], C = grid$C[g], sigma = grid$sigma[g])
          preds[va_in] <- predict_svr(m, Xtr[va_in, , drop = FALSE])
        }
      }
      val_scores$R2[g] <- rsq(ytr, preds)
    }
    
    # Best inner R2 across hyperparameters
    best_inner_R2 <- max(val_scores$R2)
    
    # Pick best by highest R2; in case of ties, choose by MAE
    if (algo == 'PLSR') {
      best_idx <- which.max(val_scores$R2)
      best <- val_scores[best_idx,]
      final_model <- fit_plsr(Xtr, ytr, ncomp = best$K)
      oof_pred[te] <- predict_plsr(final_model, Xte)
      oof_inner_r2[te] <- best_inner_R2
      chosen_hparams[[i]] <- list(K = best$K)
    } else {
      # tie-breaker by MAE
      top_R2 <- max(val_scores$R2)
      top <- which(val_scores$R2 >= top_R2 - 1e-12)
      if (length(top) > 1) {
        maes <- map_dbl(top, function(j){
          preds <- rep(NA_real_, length(ytr))
          for (f in seq_along(inner_folds)) {
            tr_in <- inner_folds[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_svr(Xtr[tr_in, , drop = FALSE], ytr[tr_in], C = grid$C[j], sigma = grid$sigma[j])
            preds[va_in] <- predict_svr(m, Xtr[va_in, , drop = FALSE])
          }
          Metrics::mae(ytr, preds)
        })
        best_idx <- top[ which.min(maes) ]
      } else best_idx <- top
      
      best <- val_scores[best_idx,]
      final_model <- fit_svr(Xtr, ytr, C = grid$C[best_idx], sigma = grid$sigma[best_idx])
      oof_pred[te] <- predict_svr(final_model, Xte)
      oof_inner_r2[te] <- best_inner_R2
      chosen_hparams[[i]] <- list(C = grid$C[best_idx], sigma = grid$sigma[best_idx])
    }
    
    oof_fold[te] <- i
  }
  
  tibble(idx = seq_along(y), pred = oof_pred, fold = oof_fold, inner_R2 = oof_inner_r2)
}

# -----------------------------
# Main
# -----------------------------
if (SAMPLE_MODE) {
  message('SAMPLE_MODE is TRUE: generating synthetic dataset...')
  raw_df <- generate_sample_data()
} else {
  raw_df <- read_data(FILE_PATH, SHEET_NAME)
  # Normalize headers and ensure TARGET_COL is present (case-insensitive, trimmed)
  orig_names <- names(raw_df)
  norm_names <- stringr::str_trim(orig_names)
  norm_lower <- tolower(norm_names)
  target_lower <- tolower(TARGET_COL)
  if (target_lower %in% norm_lower) {
    match_idx <- which(norm_lower == target_lower)[1]
    norm_names[match_idx] <- TARGET_COL
    names(raw_df) <- norm_names
  } else {
    stop(sprintf("Target column '%s' not found in sheet '%s'. Available columns: %s",
                 TARGET_COL, SHEET_NAME, paste(orig_names, collapse = ", ")))
  }
}

# Identify spectra
wl_cols <- get_spectral_cols(raw_df)
wl <- as.numeric(wl_cols)
X_raw <- as.matrix(raw_df[, wl_cols, drop = FALSE])

# Remove noisy region (700–705 nm)
exc <- exclude_noise_range(X_raw, wl, NOISE_RANGE)
X_raw <- exc$X
wl    <- exc$wl

# Quick spectra plot (first 5)
if (nrow(X_raw) >= 1) {
  n_show <- min(5, nrow(X_raw))
  plot_df <- tibble(sample = rep(1:n_show, each = length(wl)),
                    wl = rep(wl, times = n_show),
                    value = as.vector(t(X_raw[1:n_show, ])))
  ggplot(plot_df, aes(wl, value, group = sample)) +
    geom_line() +
    labs(x = 'Wavelength (nm)', y = 'Reflectance', title = 'First spectra (noise band excluded)') +
    theme_minimal()
}

# Use all samples (no QC) for the script main flow
keep_idx <- seq_len(nrow(X_raw))
df <- raw_df[keep_idx, , drop = FALSE]
X  <- X_raw[keep_idx, , drop = FALSE]
Y  <- df[[TARGET_COL]]

# Guardrail: ensure enough samples for CV
if (length(Y) < 3) {
  stop(sprintf('Too few samples available: %d. Need at least 3 to run CV.', length(Y)))
}

# Define pipelines (preprocessing x algorithm)
pipes <- tribble(
  ~prep_name, ~prep_fn, ~algo,
  'RAW', prep_raw, 'PLSR',
  'ABS', prep_abs, 'PLSR',
  'SNV', prep_snv_abs, 'PLSR',
  'SG1', prep_sg1_abs, 'PLSR',
  'SG2', prep_sg2_abs, 'PLSR',
  'RAW', prep_raw, 'SVR',
  'ABS', prep_abs, 'SVR',
  'SNV', prep_snv_abs, 'SVR',
  'SG1', prep_sg1_abs, 'SVR',
  'SG2', prep_sg2_abs, 'SVR'
)

results <- vector('list', nrow(pipes))

for (i in seq_len(nrow(pipes))) {
  cat(sprintf('\nRunning pipeline %d/%d: %s + %s\n', i, nrow(pipes), pipes$prep_name[i], pipes$algo[i]))
  oof <- run_nested_cv(X, Y, preproc_fn = pipes$prep_fn[[i]], algo = pipes$algo[i],
                       outer_k = OUTER_K, outer_rep = OUTER_REP, inner_k = INNER_K, seed = 42)
  results[[i]] <- oof %>% mutate(prep = pipes$prep_name[i], algo = pipes$algo[i])
}

res_df <- bind_rows(results) %>% left_join(tibble(idx = seq_along(Y), truth = Y), by = 'idx')

# Aggregate metrics per pipeline
summ <- res_df %>%
  group_by(prep, algo) %>%
  summarize(
    RSQ_val = rsq(truth, pred),
    MAE_val = Metrics::mae(truth, pred),
    RMSE_val = Metrics::rmse(truth, pred),
    INNER_R2_mean = mean(inner_R2, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(RSQ_val), MAE_val)

print(summ)

# Plot: OOF predictions vs truth for the best pipeline
best <- summ %>% slice(1)
best_key <- paste(best$prep, best$algo, sep = '__')
res_df <- res_df %>% mutate(key = paste(prep, algo, sep='__'))

best_oof <- res_df %>% filter(key == best_key)

ggplot(best_oof, aes(truth, pred)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = 'lm', se = FALSE) +
  labs(title = sprintf('OOF predictions: %s + %s (R^2=%.3f, MAE=%.3f)', best$prep, best$algo, best$RSQ_val, best$MAE_val),
       x = 'Truth', y = 'OOF prediction') +
  theme_minimal()
