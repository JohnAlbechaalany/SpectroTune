#' SpectroTune: Spectral pretreatment, nested CV benchmarking, plotting, and export
#'
#' SpectroTune streamlines spectral machine learning:
#' pretreatments, wavelength range controls, nested cross‑validation for
#' regression and classification, plotting best pipelines, and exporting final models.
#'
#' @section Pretreatments:
#' \itemize{
#'   \item \strong{RAW}, \strong{ABS} (-log10), \strong{SNV}, \strong{MSC}
#'   \item Savitzky–Golay: \strong{SG0} (smoothing), \strong{SG1}, \strong{SG2} (derivatives)
#'   \item Tunable via \code{preprocess_params} (e.g., \code{w}, \code{p}, \code{m})
#' }
#'
#' @section Algorithms:
#' \itemize{
#'   \item \emph{Regression}: PLSR, SVR (RBF), GLMNET (Elastic Net), RIDGE, Kernel Ridge (RBF),
#'         Random Forest, XGBoost on PCA scores (XGB_PCA), Random Forest on PCA scores (RF_PCA)
#'   \item \emph{Classification}: SVM (RBF), GLMNET (logistic/multinomial), Random Forest,
#'         XGBoost on PCA scores (XGB_PCA), PLS-DA, Random Forest on PCA scores (RF_PCA)
#' }
#'
#' @section Sampling:
#' \itemize{
#'   \item Dependent (numeric y): \code{"random"} or \code{"stratified_bin"} (quantile bins via \code{dependent_bins})
#'   \item Independent covariates: \code{"random"} or \code{"stratified"} (balance \code{stratify} factors)
#' }
#'
#' @section User-facing functions:
#' \itemize{
#'   \item \code{st_preprocess()}: apply pretreatments and wavelength range controls
#'   \item \code{st_model()}: run nested CV on provided pretreatments (advanced)
#'   \item \code{st_run()}: one-call wrapper (preprocess + model)
#'   \item \code{st_plot_best()}: plot best pipeline (scatter or confusion matrix)
#'   \item \code{st_fit_final()}, \code{st_save_final()}: refit and save a final model
#' }
#'
#' @section Getting started:
#' See \code{?st_run} examples for regression and classification end-to-end usage.
#'
#' @keywords internal
"_PACKAGE"

