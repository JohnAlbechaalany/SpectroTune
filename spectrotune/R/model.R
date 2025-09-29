#' Nested CV modeling for PLSR and SVR over pretreatments
#' @param X_list list of matrices as returned by st_preprocess()$X
#' @param y numeric response
#' @param algos character vector among c("PLSR","SVR")
#' @param outer_k,outer_rep,inner_k integers for nested CV
#' @param seed integer
#' @param verbose logical, print progress
#' @return list(results = long tibble, summary = per-pipeline metrics, best = top row)
#' @export
st_model <- function(X_list,
                     y,
                     algos = c("PLSR","SVR"),
                     outer_k = 5,
                     outer_rep = 2,
                     inner_k = 5,
                     seed = 42,
                     verbose = FALSE) {
  stopifnot(is.list(X_list), is.numeric(y))
  algos <- match.arg(algos, choices = c("PLSR","SVR"), several.ok = TRUE)
  fit_plsr <- function(X, y, ncomp) {
    df <- as.data.frame(X); df$y <- y
    pls::plsr(y ~ ., data = df, ncomp = ncomp, validation = "none", scale = TRUE)
  }
  predict_plsr <- function(model, Xnew) as.numeric(stats::predict(model, newdata = as.data.frame(Xnew), ncomp = model$ncomp))
  fit_svr <- function(X, y, C, sigma) {
    k <- kernlab::rbfdot(sigma = sigma)
    kernlab::ksvm(as.matrix(X), y, type = "eps-svr", kernel = k, C = C, scaled = TRUE)
  }
  predict_svr <- function(model, Xnew) as.numeric(kernlab::predict(model, as.matrix(Xnew)))
  rsq <- function(truth, estimate) {
    ss_res <- sum((truth - estimate)^2)
    ss_tot <- sum((truth - mean(truth))^2)
    1 - ss_res/ss_tot
  }
  PLS_K_GRID  <- 4:18
  SVR_C_GRID  <- c(0.5, 1, 2, 4)
  SVR_SIGMA_GRID <- c(0.01, 0.1, 1, 10)
  run_nested_cv <- function(X, y, algo, inner_k, outer_k, outer_rep, seed, verbose) {
    n_obs <- length(y)
    if (n_obs < 3) stop(sprintf("Too few samples: %d", n_obs))
    outer_k_eff <- min(outer_k, max(2, n_obs))
    set.seed(seed)
    outer_splits <- vector("list", outer_rep * outer_k_eff)
    idx <- 1
    for (r in seq_len(outer_rep)) {
      folds <- caret::createFolds(y, k = outer_k_eff, returnTrain = TRUE)
      for (k in seq_len(outer_k_eff)) {
        outer_splits[[idx]] <- list(train = folds[[k]], test = setdiff(seq_along(y), folds[[k]]))
        idx <- idx + 1
      }
    }
    oof_pred <- rep(NA_real_, n_obs)
    cal_pred <- rep(NA_real_, n_obs)
    oof_fold <- rep(NA_integer_, n_obs)
    oof_inner_r2 <- rep(NA_real_, n_obs)
    pb <- NULL
    if (isTRUE(verbose)) pb <- utils::txtProgressBar(min = 0, max = length(outer_splits), style = 3)
    for (i in seq_along(outer_splits)) {
      tr <- outer_splits[[i]]$train
      te <- outer_splits[[i]]$test
      Xtr <- X[tr, , drop = FALSE]
      Xte <- X[te, , drop = FALSE]
      ytr <- y[tr]
      inner_k_eff <- min(inner_k, max(2, length(ytr)))
      inner_folds <- caret::createFolds(ytr, k = inner_k_eff, returnTrain = TRUE)
      if (algo == "PLSR") {
        max_ncomp <- max(1, min(ncol(Xtr), length(ytr) - 1))
        Ks <- intersect(PLS_K_GRID, seq_len(max_ncomp))
        if (length(Ks) == 0) Ks <- max(1, min(PLS_K_GRID))
        grid <- tibble::tibble(model = "PLSR", K = Ks)
      } else {
        grid <- tidyr::expand_grid(model = "SVR", C = SVR_C_GRID, sigma = SVR_SIGMA_GRID)
      }
      val_scores <- dplyr::mutate(grid, R2 = NA_real_)
      for (g in seq_len(nrow(grid))) {
        preds <- rep(NA_real_, length(ytr))
        for (f in seq_along(inner_folds)) {
          tr_in <- inner_folds[[f]]
          va_in <- setdiff(seq_along(ytr), tr_in)
          if (algo == "PLSR") {
            m <- fit_plsr(Xtr[tr_in, , drop = FALSE], ytr[tr_in], ncomp = grid$K[g])
            preds[va_in] <- predict_plsr(m, Xtr[va_in, , drop = FALSE])
          } else {
            m <- fit_svr(Xtr[tr_in, , drop = FALSE], ytr[tr_in], C = grid$C[g], sigma = grid$sigma[g])
            preds[va_in] <- predict_svr(m, Xtr[va_in, , drop = FALSE])
          }
        }
        val_scores$R2[g] <- rsq(ytr, preds)
      }
      best_inner_R2 <- max(val_scores$R2)
      if (algo == "PLSR") {
        best_idx <- which.max(val_scores$R2)
        final_model <- fit_plsr(Xtr, ytr, ncomp = val_scores$K[best_idx])
        oof_pred[te] <- predict_plsr(final_model, Xte)
        cal_pred[tr] <- predict_plsr(final_model, Xtr)
        oof_inner_r2[te] <- best_inner_R2
      } else {
        top_R2 <- max(val_scores$R2)
        top <- which(val_scores$R2 >= top_R2 - 1e-12)
        if (length(top) > 1) {
          maes <- purrr::map_dbl(top, function(j){
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
        final_model <- fit_svr(Xtr, ytr, C = grid$C[best_idx], sigma = grid$sigma[best_idx])
        oof_pred[te] <- predict_svr(final_model, Xte)
        cal_pred[tr] <- predict_svr(final_model, Xtr)
        oof_inner_r2[te] <- best_inner_R2
      }
      oof_fold[te] <- i
      if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    if (!is.null(pb)) close(pb)
    tibble::tibble(idx = seq_along(y), pred = oof_pred, cal_pred = cal_pred, fold = oof_fold, inner_R2 = oof_inner_r2)
  }
  results <- list()
  for (prep_name in names(X_list)) {
    for (algo in algos) {
      if (isTRUE(verbose)) cat(sprintf("\nRunning: %s + %s\n", prep_name, algo))
      oof <- run_nested_cv(X_list[[prep_name]], y, algo, inner_k, outer_k, outer_rep, seed, verbose)
      results[[length(results) + 1]] <- dplyr::mutate(oof, prep = prep_name, algo = algo)
    }
  }
  res_df <- dplyr::bind_rows(results) %>% dplyr::left_join(tibble::tibble(idx = seq_along(y), truth = y), by = "idx")
  rsq <- function(truth, estimate) { ss_res <- sum((truth - estimate)^2); ss_tot <- sum((truth - mean(truth))^2); 1 - ss_res/ss_tot }
  summ <- res_df %>%
    dplyr::group_by(prep, algo) %>%
    dplyr::summarize(
      SEC = Metrics::rmse(truth[!is.na(cal_pred)], cal_pred[!is.na(cal_pred)]),
      RSQcal = rsq(truth[!is.na(cal_pred)], cal_pred[!is.na(cal_pred)]),
      SEP = Metrics::rmse(truth, pred),
      RSQval = rsq(truth, pred),
      Bias = mean(pred - truth, na.rm = TRUE),
      .groups = "drop"
    ) %>% dplyr::arrange(dplyr::desc(RSQval), SEP)
  list(results = res_df, summary = summ, best = summ %>% dplyr::slice(1))
}
