#' One-call convenience: preprocess + model
#'
#' Runs spectral pretreatment(s) and nested cross‑validation with selected algorithms.
#'
#' @param data data.frame with spectral columns named by wavelength (e.g., "400","405",...)
#' @param target character, name of the target column in `data` (numeric → regression; factor/character → classification)
#' @param preprocess character vector of pretreatments to evaluate. Supported: "RAW","ABS","SNV","SG0","SG1","SG2","MSC".
#'   You can also pass chained pipelines like "ABS+SG1+SNV" or use [combos()] to generate combinations.
#' @param preprocess_params named list of parameter lists per pretreatment, e.g. `list(SG1 = list(w=17,p=2,m=1))`.
#'   For chained pipelines, step-specific parameters are looked up by step name (e.g., `SG1`).
#' @param noise_range optional numeric length-2 range to exclude (e.g., c(700,705))
#' @param wl_keep optional length-2 numeric or list of ranges to keep
#' @param algos character vector of algorithms. Regression: "PLSR","SVR","GLMNET","RIDGE","KRR_RBF","RF","XGB_PCA","RF_PCA". Classification: "SVM","GLMNET","RF","XGB_PCA","PLSDA","RF_PCA".
#' @param outer_k,outer_rep,inner_k integers controlling outer/inner folds and repeats
#' @param seed integer seed for reproducibility
#' @param verbose logical, print progress bars/messages
#'
#' @return list with: `results` (OOF predictions), `summary` (metrics per pipeline), `best` (top row)
#'
#' @examples
#' # Regression
#' # fit <- st_run(df, target = "Age", preprocess = c("ABS","SG1"),
#' #               preprocess_params = list(SG1 = list(w=17,p=2,m=1)),
#' #               algos = c("PLSR","SVR"), outer_k=5, inner_k=5, outer_rep=2)
#'
#' # Classification
#' # df$Class <- as.factor(df$Class)
#' # fit <- st_run(df, target = "Class", preprocess = c("ABS","SG1"),
#' #               algos = c("SVM","GLMNET"), outer_k=5, inner_k=5, outer_rep=2)
#'
#' @export
st_run <- function(data,
                   target,
                   preprocess = c("RAW","ABS","SNV","SG0","SG1","SG2","MSC"),
                   preprocess_params = list(),
                   noise_range = NULL,
                   wl_keep = NULL,
                   algos = c("PLSR","SVR"),
                   outer_k = 5,
                   outer_rep = 2,
                   inner_k = 5,
                   seed = 42,
                   verbose = FALSE,
                   dependent_sampling = c("random","stratified_bin"),
                   dependent_bins = 5,
                   independent_sampling = c("random","stratified"),
                   stratify = NULL) {
  dependent_sampling <- match.arg(dependent_sampling)
  independent_sampling <- match.arg(independent_sampling)
  # detect chained pipelines (contain '+')
  has_chain <- any(grepl("\\+", preprocess))
  if (!has_chain) {
    pp <- st_preprocess(
      data = data,
      target = target,
      preprocess = preprocess,
      preprocess_params = preprocess_params,
      noise_range = noise_range,
      wl_keep = wl_keep
    )
    return(st_model(
      X_list = pp$X,
      y = pp$y,
      algos = algos,
      outer_k = outer_k,
      outer_rep = outer_rep,
      inner_k = inner_k,
      seed = seed,
      verbose = verbose,
      dependent_sampling = dependent_sampling,
      dependent_bins = dependent_bins,
      independent_sampling = independent_sampling,
      stratify = stratify
    ))
  }

  # For chained pipelines: prepare RAW once (with range handling), then compose
  pp0 <- st_preprocess(
    data = data,
    target = target,
    preprocess = c("RAW"),
    preprocess_params = list(),
    noise_range = noise_range,
    wl_keep = wl_keep
  )
  X0 <- pp0$X[["RAW"]]
  y  <- pp0$y

  # Map of step functions defined in the package
  step_map <- list(
    RAW = prep_raw,
    ABS = prep_abs,
    SNV = prep_snv_abs,
    SG0 = prep_sg0_abs,
    SG1 = prep_sg1_abs,
    SG2 = prep_sg2_abs,
    MSC = prep_msc_abs
  )
  build_chain_fn <- function(steps) {
    function(X) {
      Xc <- X
      for (s in steps) {
        fn <- step_map[[s]]
        if (is.null(fn)) stop(sprintf("Unknown preprocess step: %s", s))
        params <- preprocess_params[[s]]
        if (is.null(params)) params <- list()
        Xc <- do.call(fn, c(list(X = Xc), params))
      }
      Xc
    }
  }
  # Normalize preprocess vector (e.g., from combos()) into character chains
  chains <- as.character(preprocess)
  X_list <- list()
  for (ch in chains) {
    steps <- unlist(strsplit(ch, "\\+"))
    steps <- steps[steps != ""]
    fn <- build_chain_fn(steps)
    X_list[[ch]] <- fn(X0)
  }
  st_model(
    X_list = X_list,
    y = y,
    algos = algos,
    outer_k = outer_k,
    outer_rep = outer_rep,
    inner_k = inner_k,
    seed = seed,
    verbose = verbose,
    dependent_sampling = dependent_sampling,
    dependent_bins = dependent_bins,
    independent_sampling = independent_sampling,
    stratify = stratify
  )
}

#' Generate chained pretreatment combinations
#'
#' Convenience helper to build vectors of chained pipelines (e.g., "ABS+SG1", "ABS+SG2+SNV").
#'
#' @param parts a list of character vectors; each element represents a choice set for a step.
#'   Use \code{NULL} or \code{NA} inside a choice set to indicate an optional step.
#' @return a character vector of pipeline strings joined by '+'
#' @examples
#' # ABS; then choose SG1 or SG2; optionally add SNV
#' # combos(list(c("ABS"), c("SG1","SG2"), c("SNV", NA)))
#' @export
combos <- function(parts) {
  if (!is.list(parts)) stop("parts must be a list of character vectors")
  # replace NULL with NA for expand.grid compatibility
  parts2 <- lapply(parts, function(v) if (is.null(v)) NA_character_ else v)
  eg <- do.call(expand.grid, c(parts2, stringsAsFactors = FALSE))
  out <- apply(eg, 1, function(row) {
    steps <- as.character(row)
    steps <- steps[!is.na(steps) & steps != ""]
    paste(steps, collapse = "+")
  })
  unique(out[nchar(out) > 0])
}

#' Plot best pipeline results (regression or classification)
#'
#' Selects the best pipeline (optionally within filters) and renders a task‑specific plot:
#' regression → truth vs OOF prediction; classification → confusion matrix heatmap.
#'
#' @param fit object returned by [st_run()] or [st_model()]
#' @param task one of "auto","regression","classification"; auto infers from `fit`
#' @param select_by selection metric. Regression: "RSQval" (default), "SEP", "Bias". Classification: "Accuracy" (default), "Balanced_Accuracy", "F1_macro".
#' @param preprocess_filter optional character vector limiting pretreatments considered
#' @param algo_filter optional character vector limiting algorithms considered
#' @param return_data logical; if TRUE, returns a list with `plot`, filtered `data`, and chosen `selection`
#'
#' @return a ggplot object (or list if `return_data=TRUE`)
#'
#' @examples
#' # p <- st_plot_best(fit)
#'
#' @export
st_plot_best <- function(fit,
                         task = c("auto","regression","classification"),
                         select_by = NULL,
                         preprocess_filter = NULL,
                         algo_filter = NULL,
                         return_data = FALSE,
                         normalize = c("row","none","col","all")) {
  stopifnot(is.list(fit), !is.null(fit$results), !is.null(fit$summary))
  task <- match.arg(task)
  normalize <- match.arg(normalize)
  res_df <- fit$results
  summ   <- fit$summary
  # detect task
  if (task == "auto") {
    if ("RSQval" %in% names(summ) || "SEP" %in% names(summ)) task <- "regression" else task <- "classification"
  }
  # filter
  if (!is.null(preprocess_filter)) { res_df <- dplyr::filter(res_df, prep %in% preprocess_filter); summ <- dplyr::filter(summ, prep %in% preprocess_filter) }
  if (!is.null(algo_filter)) { res_df <- dplyr::filter(res_df, algo %in% algo_filter); summ <- dplyr::filter(summ, algo %in% algo_filter) }
  stopifnot(nrow(summ) > 0)
  # choose metric
  if (is.null(select_by)) select_by <- if (task == "regression") "RSQval" else "Accuracy"
  # rank and select
  if (task == "regression") {
    ord <- dplyr::arrange(summ, dplyr::desc(.data[[select_by]]), if ("SEP" %in% names(summ)) SEP else .data[[select_by]])
  } else {
    ord <- dplyr::arrange(summ, dplyr::desc(.data[[select_by]]), if ("Balanced_Accuracy" %in% names(summ)) Balanced_Accuracy else .data[[select_by]])
  }
  best <- ord[1, , drop = FALSE]
  key <- paste(best$prep, best$algo, sep = "__")
  res_df_key <- dplyr::mutate(res_df, key = paste(prep, algo, sep = "__"))
  best_df <- dplyr::filter(res_df_key, key == !!key)
  if (task == "regression") {
    title_txt <- sprintf("%s + %s  (RSQval=%.3f, SEP=%.3f)", best$prep, best$algo, ifelse("RSQval" %in% names(best), best$RSQval, NA_real_), ifelse("SEP" %in% names(best), best$SEP, NA_real_))
    p <- ggplot2::ggplot(best_df, ggplot2::aes(truth, pred)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray50") +
      ggplot2::geom_smooth(method = "lm", se = FALSE) +
      ggplot2::labs(title = title_txt, x = "Truth", y = "OOF prediction") +
      ggplot2::theme_minimal()
  } else {
    # Ensure both truth and pred share the same, explicit levels
    lvl <- sort(unique(c(as.character(best_df$truth), as.character(best_df$pred))))
    best_df$truth <- factor(best_df$truth, levels = lvl)
    best_df$pred  <- factor(best_df$pred,  levels = lvl)

    # Build confusion matrix including zero-count cells
    cm <- dplyr::count(best_df, truth, pred, name = "n")
    cm <- tidyr::complete(cm,
                          truth = levels(best_df$truth),
                          pred  = levels(best_df$pred),
                          fill  = list(n = 0L))

    # Compute normalization
    if (normalize == "none") {
      cm$fill <- cm$n
      cm$label <- ifelse(cm$n == 0, "", as.character(cm$n))
      fill_name <- "count"
      subtitle  <- "Confusion matrix (counts)"
    } else if (normalize == "row") {
      cm <- dplyr::group_by(cm, truth)
      cm <- dplyr::mutate(cm, denom = sum(n), pct = ifelse(denom > 0, n/denom, 0))
      cm <- dplyr::ungroup(cm)
      cm$fill  <- cm$pct
      cm$label <- scales::percent(cm$pct, accuracy = 1)
      fill_name <- "row %"
      subtitle  <- "Row-normalized (per observed class)"
    } else if (normalize == "col") {
      cm <- dplyr::group_by(cm, pred)
      cm <- dplyr::mutate(cm, denom = sum(n), pct = ifelse(denom > 0, n/denom, 0))
      cm <- dplyr::ungroup(cm)
      cm$fill  <- cm$pct
      cm$label <- scales::percent(cm$pct, accuracy = 1)
      fill_name <- "col %"
      subtitle  <- "Column-normalized (per predicted class)"
    } else {
      total_n <- sum(cm$n)
      cm$pct <- ifelse(total_n > 0, cm$n/total_n, 0)
      cm$fill  <- cm$pct
      cm$label <- scales::percent(cm$pct, accuracy = 1)
      fill_name <- "overall %"
      subtitle  <- "Overall-normalized"
    }

    p <- ggplot2::ggplot(cm, ggplot2::aes(pred, truth, fill = fill)) +
      ggplot2::geom_tile(color = "grey80") +
      ggplot2::geom_text(ggplot2::aes(label = label)) +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue", name = fill_name, limits = if (normalize == "none") NULL else c(0,1)) +
      ggplot2::labs(title = sprintf("%s + %s  (Acc=%.3f)", best$prep, best$algo, ifelse("Accuracy" %in% names(best), best$Accuracy, NA_real_)), subtitle = subtitle, x = "Predicted", y = "Observed") +
      ggplot2::theme_minimal()
  }
  if (isTRUE(return_data)) return(list(plot = p, data = best_df, selection = best)) else return(p)
}

#' Fit final model on full data for a chosen or best pipeline
#'
#' Refits a single pipeline (pretreatment + algorithm) on the full dataset,
#' selecting hyperparameters by inner CV, and returns a serializable model object.
#'
#' @param data data.frame with spectral columns
#' @param target character, target column name (numeric → regression; factor/character → classification)
#' @param preprocess single pretreatment method (e.g., "SG1")
#' @param algo single algorithm (see [st_run()] for supported values)
#' @param preprocess_params list of parameters for the pretreatment (e.g., list(w=17,p=2,m=1))
#' @param inner_k inner folds for hyperparameter selection
#' @param seed integer
#'
#' @return an object of class `spectrotune_model` with fields: `preprocess`, `preprocess_params`, `algo`, `model`, `meta`
#'
#' @examples
#' # final <- st_fit_final(df, target = "Age", preprocess = "SG1", algo = "SVR",
#' #                       preprocess_params = list(w=17,p=2,m=1), inner_k=5)
#' # st_save_final(final, "spectrotune_final_model.rds")
#'
#' @export
st_fit_final <- function(data, target,
                         preprocess,
                         algo,
                         preprocess_params = list(),
                         inner_k = 5,
                         seed = 42) {
  stopifnot(is.data.frame(data), length(preprocess) == 1, length(algo) == 1)
  # preprocess full data
  pp <- st_preprocess(data, target = target, preprocess = preprocess, preprocess_params = setNames(list(preprocess_params), preprocess))
  X <- pp$X[[preprocess]]; y <- pp$y
  n <- nrow(X)
  set.seed(seed)
  # detect task
  is_cls <- is.factor(y) || is.character(y)
  # simple inner CV on full data for hyperparameter selection
  folds <- caret::createFolds(if (is_cls) as.factor(y) else seq_len(n), k = min(inner_k, max(2, n)), returnTrain = TRUE)
  # helpers: reuse minimal versions
  fit_plsr <- function(X, y, ncomp) { df <- as.data.frame(X); df$y <- y; pls::plsr(y ~ ., data = df, ncomp = ncomp, validation = "none", scale = TRUE) }
  predict_plsr <- function(model, Xnew) as.numeric(stats::predict(model, newdata = as.data.frame(Xnew), ncomp = model$ncomp))
  fit_svr_r <- function(X, y, C, sigma) { k <- kernlab::rbfdot(sigma = sigma); kernlab::ksvm(as.matrix(X), y, type = "eps-svr", kernel = k, C = C, scaled = TRUE) }
  pred_svr_r <- function(m, Xn) as.numeric(kernlab::predict(m, as.matrix(Xn)))
  fit_svm_c <- function(X, y, C, sigma) { k <- kernlab::rbfdot(sigma = sigma); kernlab::ksvm(as.matrix(X), as.factor(y), type = "C-svc", kernel = k, C = C, scaled = TRUE) }
  pred_svm_c <- function(m, Xn) predict(m, as.matrix(Xn))
  fit_glmnet_r <- function(X, y, alpha) glmnet::cv.glmnet(as.matrix(X), y, alpha = alpha, standardize = TRUE)
  pred_glmnet_r <- function(m, Xn) as.numeric(predict(m, newx = as.matrix(Xn), s = m$lambda.min))
  fit_glmnet_c <- function(X, y, alpha) glmnet::cv.glmnet(as.matrix(X), as.factor(y), family = if (nlevels(as.factor(y))==2) "binomial" else "multinomial", alpha = alpha, type.measure = "class", standardize = TRUE)
  pred_glmnet_c <- function(m, Xn) { if (m$glmnet.fit$call$family == "binomial") as.factor(ifelse(as.numeric(predict(m, newx = as.matrix(Xn), s = m$lambda.min, type = "response"))>0.5, levels(m$glmnet.fit$y)[2], levels(m$glmnet.fit$y)[1])) else as.factor(predict(m, newx = as.matrix(Xn), s = m$lambda.min, type = "class")) }
  fit_rf <- function(X, y, mtry, num.trees, cls = FALSE) ranger::ranger(y ~ ., data = cbind.data.frame(y = if (cls) as.factor(y) else y, as.data.frame(X)), mtry = mtry, num.trees = num.trees, probability = cls)
  pred_rf_r <- function(m, Xn) as.numeric(predict(m, data = as.data.frame(Xn))$predictions)
  pred_rf_c <- function(m, Xn) predict(m, data = as.data.frame(Xn))$predictions
  # grids
  PLS_K_GRID  <- 4:18; SVR_C_GRID  <- c(0.5,1,2,4); SVR_SIGMA_GRID <- c(0.01,0.1,1,10)
  GLMNET_ALPHA_GRID <- c(0,0.25,0.5,0.75,1); RF_MTRY_GRID <- c(8,16,32); RF_TREES_GRID <- c(300,600)
  select_best <- function(score_vec, which_max = TRUE) if (which_max) which.max(score_vec) else which.min(score_vec)
  if (!is_cls) {
    # regression
    if (algo == "PLSR") {
      Ks <- intersect(PLS_K_GRID, seq_len(min(ncol(X), nrow(X)-1)))
      scores <- sapply(Ks, function(K) {
        preds <- rep(NA_real_, nrow(X))
        for (f in seq_along(folds)) {
          tr <- folds[[f]]; va <- setdiff(seq_len(nrow(X)), tr)
          m <- fit_plsr(X[tr,,drop=FALSE], y[tr], ncomp = K)
          preds[va] <- predict_plsr(m, X[va,,drop=FALSE])
        }
        1 - sum((y - preds)^2) / sum((y - mean(y))^2)
      })
      bestK <- Ks[ select_best(scores, TRUE) ]; final <- fit_plsr(X, y, ncomp = bestK); hp <- list(K = bestK)
    } else if (algo == "SVR") {
      grid <- expand.grid(C = SVR_C_GRID, sigma = SVR_SIGMA_GRID)
      scores <- apply(grid, 1, function(g){
        preds <- rep(NA_real_, nrow(X))
        for (f in seq_along(folds)) {
          tr <- folds[[f]]
          va <- setdiff(seq_len(nrow(X)), tr)
          m <- fit_svr_r(X[tr,,drop=FALSE], y[tr], C = g[["C"]], sigma = g[["sigma"]])
          preds[va] <- pred_svr_r(m, X[va,,drop=FALSE])
        }
        1 - sum((y - preds)^2) / sum((y - mean(y))^2)
      })
      j <- select_best(scores, TRUE); final <- fit_svr_r(X, y, C = grid$C[j], sigma = grid$sigma[j]); hp <- list(C = grid$C[j], sigma = grid$sigma[j])
    } else if (algo == "GLMNET" || algo == "RIDGE") {
      alphas <- if (algo == "RIDGE") c(0) else GLMNET_ALPHA_GRID
      scores <- sapply(alphas, function(a) {
        preds <- rep(NA_real_, nrow(X))
        for (f in seq_along(folds)) {
          tr <- folds[[f]]
          va <- setdiff(seq_len(nrow(X)), tr)
          cv <- fit_glmnet_r(X[tr,,drop=FALSE], y[tr], alpha = a)
          preds[va] <- pred_glmnet_r(cv, X[va,,drop=FALSE])
        }
        1 - sum((y - preds)^2) / sum((y - mean(y))^2)
      })
      a <- alphas[ select_best(scores, TRUE) ]; final <- fit_glmnet_r(X, y, alpha = a); hp <- list(alpha = a, lambda = final$lambda.min)
    } else if (algo == "RF") {
      grid <- expand.grid(mtry = RF_MTRY_GRID, trees = RF_TREES_GRID)
      scores <- apply(grid, 1, function(g) {
        preds <- rep(NA_real_, nrow(X))
        for (f in seq_along(folds)) {
          tr <- folds[[f]]
          va <- setdiff(seq_len(nrow(X)), tr)
          m <- fit_rf(X[tr,,drop=FALSE], y[tr], mtry = g[["mtry"]], num.trees = g[["trees"]], cls = FALSE)
          preds[va] <- pred_rf_r(m, X[va,,drop=FALSE])
        }
        1 - sum((y - preds)^2) / sum((y - mean(y))^2)
      })
      j <- select_best(scores, TRUE); final <- fit_rf(X, y, mtry = grid$mtry[j], num.trees = grid$trees[j], cls = FALSE); hp <- list(mtry = grid$mtry[j], trees = grid$trees[j])
    } else {
      stop("Final fit currently supports PLSR, SVR, GLMNET/RIDGE, RF for regression.")
    }
  } else {
    # classification (Accuracy)
    y <- as.factor(y)
    if (algo == "SVM") {
      grid <- expand.grid(C = SVR_C_GRID, sigma = SVR_SIGMA_GRID)
      acc <- apply(grid, 1, function(g) {
        preds <- factor(rep(NA_character_, nrow(X)), levels = levels(y))
        for (f in seq_along(folds)) {
          tr <- folds[[f]]
          va <- setdiff(seq_len(nrow(X)), tr)
          m <- fit_svm_c(X[tr,,drop=FALSE], y[tr], C = g[["C"]], sigma = g[["sigma"]])
          preds[va] <- pred_svm_c(m, X[va,,drop=FALSE])
        }
        mean(preds == y)
      })
      j <- select_best(acc, TRUE); final <- fit_svm_c(X, y, C = grid$C[j], sigma = grid$sigma[j]); hp <- list(C = grid$C[j], sigma = grid$sigma[j])
    } else if (algo == "GLMNET") {
      acc <- sapply(GLMNET_ALPHA_GRID, function(a) {
        preds <- factor(rep(NA_character_, nrow(X)), levels = levels(y))
        for (f in seq_along(folds)) {
          tr <- folds[[f]]
          va <- setdiff(seq_len(nrow(X)), tr)
          cv <- fit_glmnet_c(X[tr,,drop=FALSE], y[tr], alpha = a)
          preds[va] <- pred_glmnet_c(cv, X[va,,drop=FALSE])
        }
        mean(preds == y)
      })
      a <- GLMNET_ALPHA_GRID[ select_best(acc, TRUE) ]; final <- fit_glmnet_c(X, y, alpha = a); hp <- list(alpha = a, lambda = final$lambda.min)
    } else if (algo == "RF") {
      grid <- expand.grid(mtry = RF_MTRY_GRID, trees = RF_TREES_GRID)
      acc <- apply(grid, 1, function(g) {
        preds <- factor(rep(NA_character_, nrow(X)), levels = levels(y))
        for (f in seq_along(folds)) {
          tr <- folds[[f]]
          va <- setdiff(seq_len(nrow(X)), tr)
          m <- fit_rf(X[tr,,drop=FALSE], y[tr], mtry = g[["mtry"]], num.trees = g[["trees"]], cls = TRUE)
          preds[va] <- pred_rf_c(m, X[va,,drop=FALSE])
        }
        mean(preds == y)
      })
      j <- select_best(acc, TRUE); final <- fit_rf(X, y, mtry = grid$mtry[j], num.trees = grid$trees[j], cls = TRUE); hp <- list(mtry = grid$mtry[j], trees = grid$trees[j])
    } else {
      stop("Final fit currently supports SVM, GLMNET, RF for classification.")
    }
  }
  structure(list(
    preprocess = preprocess,
    preprocess_params = preprocess_params,
    algo = algo,
    model = final,
    meta = list(target = target, n = nrow(X), seed = seed, hp = hp)
  ), class = "spectrotune_model")
}

#' Save final model to an RDS file
#'
#' @param final_model object returned by [st_fit_final()]
#' @param file path to `.rds` file
#' @return invisibly returns `file`
#' @examples
#' # st_save_final(final, "spectrotune_final_model.rds")
#' @export
st_save_final <- function(final_model, file) {
  stopifnot(inherits(final_model, "spectrotune_model"))
  saveRDS(final_model, file = file)
  invisible(file)
}
