#' Nested CV modeling for PLSR and SVR over pretreatments
#' @param X_list list of matrices as returned by st_preprocess()$X
#' @param y numeric response
#' @param algos character vector among c("PLSR","SVR","GLMNET","RIDGE","KRR_RBF","RF","XGB_PCA","RF_PCA")
#' @param outer_k,outer_rep,inner_k integers for nested CV
#' @param seed integer
#' @param verbose logical, print progress
#' @param dependent_sampling character: "random" or "stratified_bin" (for numeric y)
#' @param dependent_bins integer number of quantile bins when dependent_sampling = "stratified_bin"
#' @param independent_sampling character: "random" or "stratified"
#' @param stratify factor or list/data.frame of factors to balance when independent_sampling = "stratified"
#' @return list(results = long tibble, summary = per-pipeline metrics, best = top row)
#' @export
st_model <- function(X_list,
                     y,
                     algos = c("PLSR","SVR","GLMNET","RIDGE","KRR_RBF","RF","XGB_PCA","RF_PCA"),
                     outer_k = 5,
                     outer_rep = 2,
                     inner_k = 5,
                     seed = 42,
                     verbose = FALSE,
                     dependent_sampling = c("random","stratified_bin"),
                     dependent_bins = 5,
                     independent_sampling = c("random","stratified"),
                     stratify = NULL) {
  stopifnot(is.list(X_list))
  is_classification <- is.factor(y) || is.character(y)
  if (is_classification) {
    y <- as.factor(y)
    return(.st_model_classification(X_list, y, algos, outer_k, outer_rep, inner_k, seed, verbose,
                                    independent_sampling, stratify))
  }
  # -------------- Regression path (existing) --------------
  algos <- match.arg(algos, choices = c("PLSR","SVR","GLMNET","RIDGE","KRR_RBF","RF","XGB_PCA","RF_PCA"), several.ok = TRUE)
  dependent_sampling <- match.arg(dependent_sampling)
  independent_sampling <- match.arg(independent_sampling)
  stopifnot(is.numeric(y))

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
  fit_glmnet <- function(X, y, alpha, lambda = NULL) glmnet::glmnet(x = as.matrix(X), y = y, alpha = alpha, lambda = lambda, standardize = TRUE)
  fit_krr_rbf <- function(X, y, lambda, sigma) {
    K <- kernlab::kernelMatrix(kernlab::rbfdot(sigma = sigma), as.matrix(X))
    alpha <- solve(K + lambda * diag(nrow(K)), y)
    list(alpha = alpha, Xtrain = as.matrix(X), sigma = sigma)
  }
  predict_krr_rbf <- function(model, Xnew) {
    Kstar <- kernlab::kernelMatrix(kernlab::rbfdot(sigma = model$sigma), as.matrix(Xnew), model$Xtrain)
    as.numeric(Kstar %*% model$alpha)
  }
  fit_rf <- function(X, y, mtry, num.trees) {
    df <- as.data.frame(X)
    colnames(df) <- paste0("X", seq_len(ncol(df)))
    df$y <- y
    ranger::ranger(y ~ ., data = df, mtry = mtry, num.trees = num.trees, respect.unordered.factors = TRUE)
  }
  predict_rf <- function(model, Xnew) {
    df_new <- as.data.frame(Xnew)
    colnames(df_new) <- paste0("X", seq_len(ncol(df_new)))
    as.numeric(predict(model, data = df_new)$predictions)
  }
  # RF on PCA scores (regression)
  fit_rf_pca <- function(X, y, ncomp, mtry, num.trees) {
    pr <- stats::prcomp(scale(X, center = TRUE, scale = TRUE), center = FALSE, scale. = FALSE)
    q <- min(ncomp, ncol(pr$x))
    df_pca <- as.data.frame(pr$x[, 1:q, drop = FALSE])
    colnames(df_pca) <- paste0("PC", seq_len(ncol(df_pca)))
    df_pca$y <- y
    # Ensure mtry doesn't exceed number of features
    mtry_adj <- min(mtry, ncol(df_pca) - 1)  # -1 because y is also in the data
    mdl <- ranger::ranger(y ~ ., data = df_pca, mtry = mtry_adj, num.trees = num.trees)
    list(pca = pr, ncomp = q, rf = mdl)
  }
  predict_rf_pca <- function(model, Xnew) {
    scores <- predict(model$pca, newdata = scale(Xnew, center = TRUE, scale = TRUE))[, 1:model$ncomp, drop = FALSE]
    df_scores <- as.data.frame(scores)
    colnames(df_scores) <- paste0("PC", seq_len(ncol(df_scores)))
    as.numeric(predict(model$rf, data = df_scores)$predictions)
  }
  fit_xgb_pca <- function(X, y, ncomp, params) {
    pr <- stats::prcomp(scale(X, center = TRUE, scale = TRUE), center = FALSE, scale. = FALSE)
    q <- min(ncomp, ncol(pr$x)); dtrain <- xgboost::xgb.DMatrix(data = pr$x[, 1:q, drop = FALSE], label = y)
    bst <- xgboost::xgb.train(params = params, data = dtrain, nrounds = params$nrounds, verbose = 0)
    list(bst = bst, pca = pr, ncomp = q)
  }
  predict_xgb_pca <- function(model, Xnew) {
    scores <- predict(model$pca, newdata = scale(Xnew, center = TRUE, scale = TRUE))[, 1:model$ncomp, drop = FALSE]
    as.numeric(predict(model$bst, xgboost::xgb.DMatrix(scores)))
  }
  rsq <- function(truth, estimate) { ss_res <- sum((truth - estimate)^2); ss_tot <- sum((truth - mean(truth))^2); 1 - ss_res/ss_tot }
  to_factor <- function(x) { if (is.factor(x)) x else as.factor(x) }
  combine_stratify <- function(stratify_obj, n) {
    if (is.null(stratify_obj)) return(NULL)
    if (is.data.frame(stratify_obj) || is.list(stratify_obj)) { stopifnot(all(vapply(stratify_obj, length, integer(1)) == n)); f_list <- lapply(stratify_obj, to_factor); return(interaction(f_list, drop = TRUE)) } else { stopifnot(length(stratify_obj) == n); return(to_factor(stratify_obj)) }
  }
  make_y_bins <- function(y, bins) {
    probs <- seq(0, 1, length.out = bins + 1); brks  <- unique(stats::quantile(y, probs = probs, na.rm = TRUE, type = 7))
    if (length(brks) < 3) return(as.factor(cut(y, breaks = pretty(range(y), n = max(2, bins)), include.lowest = TRUE)))
    as.factor(cut(y, breaks = brks, include.lowest = TRUE))
  }
  build_strata <- function(y, dep_samp, dep_bins, indep_samp, stratify_obj) {
    s_indep <- if (indep_samp == "stratified") combine_stratify(stratify_obj, length(y)) else NULL
    s_dep   <- if (dep_samp == "stratified_bin") make_y_bins(y, dep_bins) else NULL
    if (!is.null(s_indep) && !is.null(s_dep)) interaction(s_indep, s_dep, drop = TRUE) else if (!is.null(s_indep)) s_indep else if (!is.null(s_dep)) s_dep else NULL
  }
  create_folds_balanced <- function(n, k, strata = NULL, seed = NULL, returnTrain = TRUE) {
    if (!is.null(seed)) set.seed(seed)
    if (is.null(strata)) return(caret::createFolds(seq_len(n), k = k, returnTrain = returnTrain))
    strata <- as.factor(strata); by_stratum <- split(seq_len(n), strata); folds_te <- vector("list", k); for (i in seq_len(k)) folds_te[[i]] <- integer(0)
    for (s in by_stratum) { s <- sample(s); parts <- split(s, rep(seq_len(k), length.out = length(s))); for (i in seq_len(k)) folds_te[[i]] <- c(folds_te[[i]], parts[[i]]) }
    if (isTRUE(returnTrain)) lapply(folds_te, function(te) setdiff(seq_len(n), te)) else folds_te
  }
  PLS_K_GRID  <- 4:18; SVR_C_GRID  <- c(0.5, 1, 2, 4); SVR_SIGMA_GRID <- c(0.01, 0.1, 1, 10)
  GLMNET_ALPHA_GRID <- c(0, 0.25, 0.5, 0.75, 1); KRR_LAMBDA_GRID <- c(1e-3, 1e-2, 1e-1, 1); KRR_SIGMA_GRID  <- c(0.01, 0.1, 1)
  RF_MTRY_GRID    <- c(2, 4, 8); RF_TREES_GRID   <- c(300, 600)
  XGB_PCA_NCOMP   <- c(20, 50, 100)
  RF_PCA_NCOMP    <- c(20, 50, 100)
  XGB_PARAM_GRID  <- list(list(eta=0.1, max_depth=3, subsample=0.8, colsample_bytree=0.8, nrounds=200), list(eta=0.05, max_depth=4, subsample=0.8, colsample_bytree=0.8, nrounds=300))
  run_nested_cv <- function(X, y, algo, inner_k, outer_k, outer_rep, seed, verbose, dep_samp, dep_bins, indep_samp, stratify_obj) {
    n_obs <- length(y); if (n_obs < 3) stop(sprintf("Too few samples: %d", n_obs)); outer_k_eff <- min(outer_k, max(2, n_obs))
    outer_strata <- build_strata(y, dep_samp, dep_bins, indep_samp, stratify_obj)
    set.seed(seed); outer_splits <- vector("list", outer_rep * outer_k_eff); idx <- 1
    for (r in seq_len(outer_rep)) { folds_tr <- create_folds_balanced(n_obs, k = outer_k_eff, strata = outer_strata, seed = seed + r, returnTrain = TRUE); for (k_i in seq_len(outer_k_eff)) { tr_idx <- folds_tr[[k_i]]; te_idx <- setdiff(seq_len(n_obs), tr_idx); outer_splits[[idx]] <- list(train = tr_idx, test = te_idx); idx <- idx + 1 } }
    oof_pred <- rep(NA_real_, n_obs); cal_pred <- rep(NA_real_, n_obs); oof_fold <- rep(NA_integer_, n_obs); oof_inner_r2 <- rep(NA_real_, n_obs)
    pb <- NULL; if (isTRUE(verbose)) pb <- utils::txtProgressBar(min = 0, max = length(outer_splits), style = 3)
    for (i in seq_along(outer_splits)) {
      tr <- outer_splits[[i]]$train; te <- outer_splits[[i]]$test; Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]; ytr <- y[tr]
      inner_k_eff <- min(inner_k, max(2, length(ytr)))
      if (!is.null(stratify_obj)) { if (is.list(stratify_obj) || is.data.frame(stratify_obj)) strat_sub <- lapply(stratify_obj, function(v) v[tr]) else strat_sub <- stratify_obj[tr] } else strat_sub <- NULL
      inner_strata <- build_strata(ytr, dep_samp, dep_bins, indep_samp, strat_sub); inner_folds_tr <- create_folds_balanced(length(ytr), k = inner_k_eff, strata = inner_strata, seed = seed + i, returnTrain = TRUE)
      if (algo == "PLSR") {
        max_ncomp <- max(1, min(ncol(Xtr), length(ytr) - 1))
        grid <- tibble::tibble(K = intersect(PLS_K_GRID, seq_len(max_ncomp)))
        val_scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- rep(NA_real_, length(ytr))
          for (f in seq_along(inner_folds_tr)) {
            tr_in <- inner_folds_tr[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_plsr(Xtr[tr_in, , drop = FALSE], ytr[tr_in], ncomp = grid$K[g])
            preds[va_in] <- predict_plsr(m, Xtr[va_in, , drop = FALSE])
          }
          val_scores[g] <- rsq(ytr, preds)
        }
        best_idx <- which.max(val_scores)
        best_inner_R2 <- val_scores[best_idx]
        final_model <- fit_plsr(Xtr, ytr, ncomp = grid$K[best_idx])
        oof_pred[te] <- predict_plsr(final_model, Xte)
        cal_pred[tr] <- predict_plsr(final_model, Xtr)
      } else if (algo == "SVR") {
        grid <- tidyr::expand_grid(C = SVR_C_GRID, sigma = SVR_SIGMA_GRID)
        val_scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- rep(NA_real_, length(ytr))
          for (f in seq_along(inner_folds_tr)) {
            tr_in <- inner_folds_tr[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_svr(Xtr[tr_in, , drop = FALSE], ytr[tr_in], C = grid$C[g], sigma = grid$sigma[g])
            preds[va_in] <- predict_svr(m, Xtr[va_in, , drop = FALSE])
          }
          val_scores[g] <- rsq(ytr, preds)
        }
        best_idx <- which.max(val_scores)
        best_inner_R2 <- val_scores[best_idx]
        final_model <- fit_svr(Xtr, ytr, C = grid$C[best_idx], sigma = grid$sigma[best_idx])
        oof_pred[te] <- predict_svr(final_model, Xte)
        cal_pred[tr] <- predict_svr(final_model, Xtr)
      } else if (algo == "GLMNET" || algo == "RIDGE") {
        alphas <- if (algo == "RIDGE") c(0) else GLMNET_ALPHA_GRID
        grid <- tibble::tibble(alpha = alphas)
        val_scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- rep(NA_real_, length(ytr))
          for (f in seq_along(inner_folds_tr)) {
            tr_in <- inner_folds_tr[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            cv <- glmnet::cv.glmnet(as.matrix(Xtr[tr_in, , drop = FALSE]), ytr[tr_in], alpha = grid$alpha[g], standardize = TRUE)
            preds[va_in] <- as.numeric(predict(cv, newx = as.matrix(Xtr[va_in, , drop = FALSE]), s = cv$lambda.min))
          }
          val_scores[g] <- rsq(ytr, preds)
        }
        best_idx <- which.max(val_scores)
        best_inner_R2 <- val_scores[best_idx]
        cv <- glmnet::cv.glmnet(as.matrix(Xtr), ytr, alpha = grid$alpha[best_idx], standardize = TRUE)
        final_model <- cv
        oof_pred[te] <- as.numeric(predict(cv, newx = as.matrix(Xte), s = cv$lambda.min))
        cal_pred[tr] <- as.numeric(predict(cv, newx = as.matrix(Xtr), s = cv$lambda.min))
      } else if (algo == "KRR_RBF") {
        grid <- tidyr::expand_grid(lambda = KRR_LAMBDA_GRID, sigma = KRR_SIGMA_GRID)
        val_scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- rep(NA_real_, length(ytr))
          for (f in seq_along(inner_folds_tr)) {
            tr_in <- inner_folds_tr[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_krr_rbf(Xtr[tr_in, , drop = FALSE], ytr[tr_in], lambda = grid$lambda[g], sigma = grid$sigma[g])
            preds[va_in] <- predict_krr_rbf(m, Xtr[va_in, , drop = FALSE])
          }
          val_scores[g] <- rsq(ytr, preds)
        }
        best_idx <- which.max(val_scores)
        best_inner_R2 <- val_scores[best_idx]
        final_model <- fit_krr_rbf(Xtr, ytr, lambda = grid$lambda[best_idx], sigma = grid$sigma[best_idx])
        oof_pred[te] <- predict_krr_rbf(final_model, Xte)
        cal_pred[tr] <- predict_krr_rbf(final_model, Xtr)
      } else if (algo == "RF") {
        grid <- tidyr::expand_grid(mtry = RF_MTRY_GRID, trees = RF_TREES_GRID)
        val_scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- rep(NA_real_, length(ytr))
          for (f in seq_along(inner_folds_tr)) {
            tr_in <- inner_folds_tr[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_rf(Xtr[tr_in, , drop = FALSE], ytr[tr_in], mtry = grid$mtry[g], num.trees = grid$trees[g])
            preds[va_in] <- predict_rf(m, Xtr[va_in, , drop = FALSE])
          }
          val_scores[g] <- rsq(ytr, preds)
        }
        best_idx <- which.max(val_scores)
        best_inner_R2 <- val_scores[best_idx]
        final_model <- fit_rf(Xtr, ytr, mtry = grid$mtry[best_idx], num.trees = grid$trees[best_idx])
        oof_pred[te] <- predict_rf(final_model, Xte)
        cal_pred[tr] <- predict_rf(final_model, Xtr)
      } else if (algo == "XGB_PCA") {
        grid <- tidyr::expand_grid(ncomp = XGB_PCA_NCOMP, param_idx = seq_along(XGB_PARAM_GRID))
        val_scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- rep(NA_real_, length(ytr))
          for (f in seq_along(inner_folds_tr)) {
            tr_in <- inner_folds_tr[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_xgb_pca(Xtr[tr_in, , drop = FALSE], ytr[tr_in], ncomp = grid$ncomp[g], params = XGB_PARAM_GRID[[grid$param_idx[g]]])
            preds[va_in] <- predict_xgb_pca(m, Xtr[va_in, , drop = FALSE])
          }
          val_scores[g] <- rsq(ytr, preds)
        }
        best_idx <- which.max(val_scores)
        best_inner_R2 <- val_scores[best_idx]
        final_model <- fit_xgb_pca(Xtr, ytr, ncomp = grid$ncomp[best_idx], params = XGB_PARAM_GRID[[grid$param_idx[best_idx]]])
        oof_pred[te] <- predict_xgb_pca(final_model, Xte)
        cal_pred[tr] <- predict_xgb_pca(final_model, Xtr)
      } else if (algo == "RF_PCA") {
        grid <- tidyr::expand_grid(ncomp = RF_PCA_NCOMP, mtry = RF_MTRY_GRID, trees = RF_TREES_GRID)
        val_scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- rep(NA_real_, length(ytr))
          for (f in seq_along(inner_folds_tr)) {
            tr_in <- inner_folds_tr[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_rf_pca(Xtr[tr_in, , drop = FALSE], ytr[tr_in], ncomp = grid$ncomp[g], mtry = grid$mtry[g], num.trees = grid$trees[g])
            preds[va_in] <- predict_rf_pca(m, Xtr[va_in, , drop = FALSE])
          }
          val_scores[g] <- rsq(ytr, preds)
        }
        best_idx <- which.max(val_scores)
        best_inner_R2 <- val_scores[best_idx]
        final_model <- fit_rf_pca(Xtr, ytr, ncomp = grid$ncomp[best_idx], mtry = grid$mtry[best_idx], num.trees = grid$trees[best_idx])
        oof_pred[te] <- predict_rf_pca(final_model, Xte)
        cal_pred[tr] <- predict_rf_pca(final_model, Xtr)
      }
      oof_fold[te] <- i; oof_inner_r2[te] <- best_inner_R2; if (isTRUE(verbose) && !is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    if (!is.null(pb)) close(pb)
    tibble::tibble(idx = seq_along(y), pred = oof_pred, cal_pred = cal_pred, fold = oof_fold, inner_R2 = oof_inner_r2)
  }
  results <- list(); for (prep_name in names(X_list)) { for (algo in algos) { if (isTRUE(verbose)) cat(sprintf("\nRunning: %s + %s\n", prep_name, algo)); oof <- run_nested_cv(X_list[[prep_name]], y, algo, inner_k, outer_k, outer_rep, seed, verbose, dependent_sampling, dependent_bins, independent_sampling, stratify); results[[length(results) + 1]] <- dplyr::mutate(oof, prep = prep_name, algo = algo) } }
  res_df <- dplyr::bind_rows(results)
  res_df <- dplyr::left_join(res_df, tibble::tibble(idx = seq_along(y), truth = y), by = "idx")
  rsq <- function(truth, estimate) { ss_res <- sum((truth - estimate)^2); ss_tot <- sum((truth - mean(truth))^2); 1 - ss_res/ss_tot }
  summ <- res_df
  summ <- dplyr::group_by(summ, prep, algo)
  summ <- dplyr::summarize(summ, SEC = Metrics::rmse(truth[!is.na(cal_pred)], cal_pred[!is.na(cal_pred)]), RSQcal = rsq(truth[!is.na(cal_pred)], cal_pred[!is.na(cal_pred)]), SEP = Metrics::rmse(truth, pred), RSQval = rsq(truth, pred), Bias = mean(pred - truth, na.rm = TRUE), .groups = "drop")
  summ <- dplyr::arrange(summ, dplyr::desc(RSQval), SEP)
  list(results = res_df, summary = summ, best = summ[1, ])
}

# Internal: classification path
.st_model_classification <- function(X_list, y, algos, outer_k, outer_rep, inner_k, seed, verbose,
                                     independent_sampling, stratify) {
  algos <- match.arg(algos, choices = c("SVM","GLMNET","RF","XGB_PCA","PLSDA","RF_PCA"), several.ok = TRUE)
  independent_sampling <- match.arg(independent_sampling, choices = c("random", "stratified"))
  n_class <- nlevels(y)
  # helpers
  fit_svm <- function(X, y, C, sigma) {
    k <- kernlab::rbfdot(sigma = sigma)
    kernlab::ksvm(as.matrix(X), y, type = "C-svc", kernel = k, C = C, scaled = TRUE)
  }
  predict_svm <- function(model, Xnew) kernlab::predict(model, as.matrix(Xnew))
  fit_glmnet <- function(X, y, alpha) {
    fam <- if (nlevels(y) == 2) "binomial" else "multinomial"
    glmnet::cv.glmnet(as.matrix(X), y, family = fam, alpha = alpha, type.measure = "class", standardize = TRUE)
  }
  predict_glmnet <- function(model, Xnew) {
    if (model$glmnet.fit$call$family == "binomial") {
      as.factor(ifelse(as.numeric(predict(model, newx = as.matrix(Xnew), s = model$lambda.min, type = "response")) > 0.5, levels(model$glmnet.fit$y)[2], levels(model$glmnet.fit$y)[1]))
    } else {
      cls <- predict(model, newx = as.matrix(Xnew), s = model$lambda.min, type = "class")
      as.factor(cls)
    }
  }
  fit_rf <- function(X, y, mtry, num.trees) {
    df <- as.data.frame(X)
    colnames(df) <- paste0("X", seq_len(ncol(df)))
    df$y <- y
    ranger::ranger(y ~ ., data = df, mtry = mtry, num.trees = num.trees, probability = FALSE)
  }
  predict_rf <- function(model, Xnew) {
    df_new <- as.data.frame(Xnew)
    colnames(df_new) <- paste0("X", seq_len(ncol(df_new)))
    predict(model, data = df_new)$predictions
  }
  # RF on PCA scores (classification)
  fit_rf_pca <- function(X, y, ncomp, mtry, num.trees) {
    pr <- stats::prcomp(scale(X, center = TRUE, scale = TRUE), center = FALSE, scale. = FALSE)
    q <- min(ncomp, ncol(pr$x))
    df_pca <- as.data.frame(pr$x[, 1:q, drop = FALSE])
    colnames(df_pca) <- paste0("PC", seq_len(ncol(df_pca)))
    df_pca$y <- y
    # Ensure mtry doesn't exceed number of features
    mtry_adj <- min(mtry, ncol(df_pca) - 1)  # -1 because y is also in the data
    mdl <- ranger::ranger(y ~ ., data = df_pca, mtry = mtry_adj, num.trees = num.trees, probability = FALSE)
    list(pca = pr, ncomp = q, rf = mdl, levels = levels(y))
  }
  predict_rf_pca <- function(model, Xnew) {
    scores <- predict(model$pca, newdata = scale(Xnew, center = TRUE, scale = TRUE))[, 1:model$ncomp, drop = FALSE]
    df_scores <- as.data.frame(scores)
    colnames(df_scores) <- paste0("PC", seq_len(ncol(df_scores)))
    predict(model$rf, data = df_scores)$predictions
  }
  # PLS-DA via mixOmics
  fit_plsda <- function(X, y, ncomp) mixOmics::plsda(X, y, ncomp = ncomp)
  predict_plsda <- function(model, Xnew) as.factor(predict(model, Xnew)$class$max.dist)
  # folds (classification: stratify by y, plus optional independent stratify)
  to_factor <- function(x) { if (is.factor(x)) x else as.factor(x) }
  combine_stratify <- function(stratify_obj, n) {
    if (is.null(stratify_obj)) return(NULL)
    if (is.data.frame(stratify_obj) || is.list(stratify_obj)) { stopifnot(all(vapply(stratify_obj, length, integer(1)) == n)); f_list <- lapply(stratify_obj, to_factor); return(interaction(f_list, drop = TRUE)) } else { stopifnot(length(stratify_obj) == n); return(to_factor(stratify_obj)) }
  }
  build_strata <- function(y, indep_samp, stratify_obj) {
    s_indep <- if (indep_samp == "stratified") combine_stratify(stratify_obj, length(y)) else NULL
    if (!is.null(s_indep)) interaction(y, s_indep, drop = TRUE) else y
  }
  create_folds_balanced <- function(n, k, strata = NULL, seed = NULL, returnTrain = TRUE) {
    if (!is.null(seed)) set.seed(seed)
    if (is.null(strata)) return(caret::createFolds(seq_len(n), k = k, returnTrain = returnTrain))
    strata <- as.factor(strata); by_stratum <- split(seq_len(n), strata); folds_te <- vector("list", k); for (i in seq_len(k)) folds_te[[i]] <- integer(0)
    for (s in by_stratum) { s <- sample(s); parts <- split(s, rep(seq_len(k), length.out = length(s))); for (i in seq_len(k)) folds_te[[i]] <- c(folds_te[[i]], parts[[i]]) }
    if (isTRUE(returnTrain)) lapply(folds_te, function(te) setdiff(seq_len(n), te)) else folds_te
  }
  # grids
  C_GRID  <- c(0.5, 1, 2, 4); SIGMA_GRID <- c(0.01, 0.1, 1)
  GLMNET_ALPHA_GRID <- c(0, 0.25, 0.5, 0.75, 1)
  RF_MTRY_GRID    <- c(2, 4, 8); RF_TREES_GRID   <- c(300, 600)
  XGB_PCA_NCOMP   <- c(20, 50, 100)
  RF_PCA_NCOMP    <- c(20, 50, 100)
  XGB_PARAM_GRID  <- list(list(eta=0.1, max_depth=3, subsample=0.8, colsample_bytree=0.8, nrounds=200), list(eta=0.05, max_depth=4, subsample=0.8, colsample_bytree=0.8, nrounds=300))
  PLSDA_NCOMP_GRID <- 1:10
  run_nested <- function(X, y, algo, inner_k, outer_k, outer_rep, seed, verbose, indep_samp, stratify_obj) {
    n_obs <- length(y); outer_k_eff <- min(outer_k, max(2, n_obs)); strata_outer <- build_strata(y, indep_samp, stratify_obj)
    set.seed(seed); outer_splits <- vector("list", outer_rep * outer_k_eff); idx <- 1
    for (r in seq_len(outer_rep)) { folds_tr <- create_folds_balanced(n_obs, k = outer_k_eff, strata = strata_outer, seed = seed + r, returnTrain = TRUE); for (k_i in seq_len(outer_k_eff)) { tr_idx <- folds_tr[[k_i]]; te_idx <- setdiff(seq_len(n_obs), tr_idx); outer_splits[[idx]] <- list(train = tr_idx, test = te_idx); idx <- idx + 1 } }
    oof_pred <- factor(rep(NA_character_, n_obs), levels = levels(y)); oof_fold <- rep(NA_integer_, n_obs)
    pb <- NULL; if (isTRUE(verbose)) pb <- utils::txtProgressBar(min = 0, max = length(outer_splits), style = 3)
    for (i in seq_along(outer_splits)) {
      tr <- outer_splits[[i]]$train; te <- outer_splits[[i]]$test; Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]; ytr <- y[tr]
      inner_k_eff <- min(inner_k, max(2, length(ytr)))
      if (!is.null(stratify_obj)) { if (is.list(stratify_obj) || is.data.frame(stratify_obj)) strat_sub <- lapply(stratify_obj, function(v) v[tr]) else strat_sub <- stratify_obj[tr] } else strat_sub <- NULL
      strata_inner <- build_strata(ytr, indep_samp, strat_sub); folds_tr_in <- create_folds_balanced(length(ytr), k = inner_k_eff, strata = strata_inner, seed = seed + i, returnTrain = TRUE)
      if (algo == "SVM") {
        grid <- tidyr::expand_grid(C = C_GRID, sigma = SIGMA_GRID)
        scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- factor(rep(NA_character_, length(ytr)), levels = levels(ytr))
          for (f in seq_along(folds_tr_in)) {
            tr_in <- folds_tr_in[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_svm(Xtr[tr_in, , drop = FALSE], ytr[tr_in], C = grid$C[g], sigma = grid$sigma[g])
            preds[va_in] <- predict_svm(m, Xtr[va_in, , drop = FALSE])
          }
          scores[g] <- mean(preds == ytr)
        }
        best <- which.max(scores)
        final <- fit_svm(Xtr, ytr, C = grid$C[best], sigma = grid$sigma[best])
        oof_pred[te] <- predict_svm(final, Xte)
      } else if (algo == "GLMNET") {
        grid <- tibble::tibble(alpha = GLMNET_ALPHA_GRID)
        scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- factor(rep(NA_character_, length(ytr)), levels = levels(ytr))
          for (f in seq_along(folds_tr_in)) {
            tr_in <- folds_tr_in[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            cv <- fit_glmnet(Xtr[tr_in, , drop = FALSE], ytr[tr_in], alpha = grid$alpha[g])
            preds[va_in] <- predict_glmnet(cv, Xtr[va_in, , drop = FALSE])
          }
          scores[g] <- mean(preds == ytr)
        }
        best <- which.max(scores)
        cv <- fit_glmnet(Xtr, ytr, alpha = grid$alpha[best])
        oof_pred[te] <- predict_glmnet(cv, Xte)
      } else if (algo == "RF") {
        grid <- tidyr::expand_grid(mtry = c(8, 16, 32), trees = c(300, 600))
        scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- factor(rep(NA_character_, length(ytr)), levels = levels(ytr))
          for (f in seq_along(folds_tr_in)) {
            tr_in <- folds_tr_in[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_rf(Xtr[tr_in, , drop = FALSE], ytr[tr_in], mtry = grid$mtry[g], num.trees = grid$trees[g])
            preds[va_in] <- predict_rf(m, Xtr[va_in, , drop = FALSE])
          }
          scores[g] <- mean(preds == ytr)
        }
        best <- which.max(scores)
        final <- fit_rf(Xtr, ytr, mtry = grid$mtry[best], num.trees = grid$trees[best])
        oof_pred[te] <- predict_rf(final, Xte)
      } else if (algo == "XGB_PCA") {
        grid <- tidyr::expand_grid(ncomp = c(20, 50, 100), param_idx = 1:2)
        scores <- rep(NA_real_, nrow(grid))
        PARAMS <- list(
          list(eta = 0.1, max_depth = 3, subsample = 0.8, colsample_bytree = 0.8, nrounds = 200),
          list(eta = 0.05, max_depth = 4, subsample = 0.8, colsample_bytree = 0.8, nrounds = 300)
        )
        for (g in seq_len(nrow(grid))) {
          preds <- factor(rep(NA_character_, length(ytr)), levels = levels(ytr))
          for (f in seq_along(folds_tr_in)) {
            tr_in <- folds_tr_in[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_xgb_pca(Xtr[tr_in, , drop = FALSE], ytr[tr_in], ncomp = grid$ncomp[g], params = PARAMS[[grid$param_idx[g]]])
            preds[va_in] <- predict_xgb_pca(m, Xtr[va_in, , drop = FALSE])
          }
          scores[g] <- mean(preds == ytr)
        }
        best <- which.max(scores)
        final <- fit_xgb_pca(Xtr, ytr, ncomp = grid$ncomp[best], params = PARAMS[[grid$param_idx[best]]])
        oof_pred[te] <- predict_xgb_pca(final, Xte)
      } else if (algo == "RF_PCA") {
        grid <- tidyr::expand_grid(ncomp = RF_PCA_NCOMP, mtry = c(2,4,8), trees = c(300,600))
        scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- factor(rep(NA_character_, length(ytr)), levels = levels(ytr))
          for (f in seq_along(folds_tr_in)) {
            tr_in <- folds_tr_in[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_rf_pca(Xtr[tr_in, , drop = FALSE], ytr[tr_in], ncomp = grid$ncomp[g], mtry = grid$mtry[g], num.trees = grid$trees[g])
            preds[va_in] <- predict_rf_pca(m, Xtr[va_in, , drop = FALSE])
          }
          scores[g] <- mean(preds == ytr)
        }
        best <- which.max(scores)
        final <- fit_rf_pca(Xtr, ytr, ncomp = grid$ncomp[best], mtry = grid$mtry[best], num.trees = grid$trees[best])
        oof_pred[te] <- predict_rf_pca(final, Xte)
      } else if (algo == "PLSDA") {
        grid <- tibble::tibble(ncomp = PLSDA_NCOMP_GRID)
        scores <- rep(NA_real_, nrow(grid))
        for (g in seq_len(nrow(grid))) {
          preds <- factor(rep(NA_character_, length(ytr)), levels = levels(ytr))
          for (f in seq_along(folds_tr_in)) {
            tr_in <- folds_tr_in[[f]]
            va_in <- setdiff(seq_along(ytr), tr_in)
            m <- fit_plsda(Xtr[tr_in, , drop = FALSE], ytr[tr_in], ncomp = grid$ncomp[g])
            preds[va_in] <- predict_plsda(m, Xtr[va_in, , drop = FALSE])
          }
          scores[g] <- mean(preds == ytr)
        }
        best <- which.max(scores)
        final <- fit_plsda(Xtr, ytr, ncomp = grid$ncomp[best])
        oof_pred[te] <- predict_plsda(final, Xte)
      }
      oof_fold[te] <- i; if (isTRUE(verbose) && !is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    if (!is.null(pb)) close(pb)
    tibble::tibble(idx = seq_along(y), pred = oof_pred, fold = oof_fold)
  }
  results <- list(); for (prep_name in names(X_list)) { for (algo in algos) { if (isTRUE(verbose)) cat(sprintf("\nRunning: %s + %s\n", prep_name, algo)); oof <- run_nested(X_list[[prep_name]], y, algo, inner_k, outer_k, outer_rep, seed, verbose, independent_sampling, stratify); results[[length(results) + 1]] <- dplyr::mutate(oof, prep = prep_name, algo = algo) } }
  res_df <- dplyr::bind_rows(results)
  res_df <- dplyr::left_join(res_df, tibble::tibble(idx = seq_along(y), truth = y), by = "idx")
  # metrics: Accuracy, Balanced Accuracy, Macro F1
  bal_acc <- function(truth, pred) {
    tab <- table(truth, pred); rec <- diag(prop.table(tab, 1)); mean(rec, na.rm = TRUE)
  }
  f1_macro <- function(truth, pred) {
    lv <- levels(truth); f1s <- sapply(lv, function(l){
      tp <- sum(truth == l & pred == l); fp <- sum(truth != l & pred == l); fn <- sum(truth == l & pred != l)
      prec <- if ((tp+fp)==0) 0 else tp/(tp+fp); rec <- if ((tp+fn)==0) 0 else tp/(tp+fn); if ((prec+rec)==0) 0 else 2*prec*rec/(prec+rec)
    }); mean(f1s)
  }
  summ <- res_df
  summ <- dplyr::group_by(summ, prep, algo)
  summ <- dplyr::summarize(summ, Accuracy = mean(pred == truth, na.rm = TRUE), Balanced_Accuracy = bal_acc(truth, pred), F1_macro = f1_macro(truth, pred), .groups = "drop")
  summ <- dplyr::arrange(summ, dplyr::desc(Accuracy), dplyr::desc(Balanced_Accuracy))
  list(results = res_df, summary = summ, best = summ[1, ])
}
