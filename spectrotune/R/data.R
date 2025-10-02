#' Sample spectral dataset for demonstrations
#'
#' A synthetic dataset containing NIR spectra and associated target variables
#' for regression and classification examples.
#'
#' @format A data frame with 150 rows and 1002 variables:
#' \describe{
#'   \item{ID}{Sample identifier}
#'   \item{Age}{Numeric target for regression (20-80 years)}
#'   \item{Class}{Factor target for classification (A, B, C)}
#'   \item{breed}{Factor covariate for stratification (Breed1, Breed2, Breed3)}
#'   \item{sex}{Factor covariate for stratification (M, F)}
#'   \item{400-2499}{Spectral reflectance values at 5nm intervals (420 wavelengths)}
#' }
#'
#' @details
#' This dataset is synthetically generated for demonstration purposes.
#' The spectra represent NIR reflectance measurements from 400-2499 nm
#' with realistic noise and baseline variations.
#'
#' @examples
#' data(sample_spectra)
#' head(sample_spectra[, 1:10])
#' summary(sample_spectra$Age)
#' table(sample_spectra$Class)
#'
#' # Use in st_run
#' \dontrun{
#' fit <- st_run(
#'   data = sample_spectra,
#'   target = "Age",
#'   preprocess = c("ABS", "SG1"),
#'   algos = c("PLSR", "SVR"),
#'   outer_k = 5,
#'   inner_k = 5,
#'   seed = 42
#' )
#' }
#'
#' @source Synthetic data generated for SpectroTune package demonstrations
"sample_spectra"
