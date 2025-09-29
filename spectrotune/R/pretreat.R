#' Spectral pretreatments with optional wavelength controls
#' @param data data.frame with spectral columns named by wavelength (e.g., "400","405",...)
#' @param preprocess character vector of methods: c("RAW","ABS","SNV","SG0","SG1","SG2","MSC")
#' @param preprocess_params named list of per-method parameter lists, e.g. list(SG1 = list(w=17,p=2,m=1))
#' @param noise_range numeric length-2 to exclude (e.g., c(700,705)) or NULL
#' @param wl_keep single length-2 or list of length-2 numeric ranges to keep, or NULL
#' @param target optional target column name to attach
#' @return list(X = named list of matrices by method, y = optional vector, wl = numeric)
#' @export
st_preprocess <- function(data,
                          preprocess = c("RAW","ABS","SNV","SG0","SG1","SG2","MSC"),
                          preprocess_params = list(),
                          noise_range = NULL,
                          wl_keep = NULL,
                          target = NULL) {
  stopifnot(is.data.frame(data))
  get_spectral_cols <- function(df) {
    cols <- colnames(df)
    cols[ suppressWarnings(!is.na(as.numeric(cols))) ]
  }
  exclude_noise_range <- function(X, wavelengths, rng) {
    keep <- !(wavelengths >= rng[1] & wavelengths <= rng[2])
    list(X = X[, keep, drop = FALSE], wl = wavelengths[keep])
  }
  keep_wavelength_ranges <- function(X, wavelengths, wl_keep = NULL) {
    if (is.null(wl_keep)) return(list(X = X, wl = wavelengths))
    if (is.list(wl_keep)) {
      mask <- Reduce(`|`, lapply(wl_keep, function(rng) wavelengths >= rng[1] & wavelengths <= rng[2]))
    } else {
      mask <- (wavelengths >= wl_keep[1] & wavelengths <= wl_keep[2])
    }
    list(X = X[, mask, drop = FALSE], wl = wavelengths[mask])
  }
  wl_cols <- get_spectral_cols(data)
  stopifnot(length(wl_cols) > 0)
  wl <- as.numeric(wl_cols)
  X <- as.matrix(data[, wl_cols, drop = FALSE])
  if (!is.null(noise_range)) {
    exc <- exclude_noise_range(X, wl, noise_range)
    X <- exc$X; wl <- exc$wl
  }
  if (!is.null(wl_keep)) {
    kept <- keep_wavelength_ranges(X, wl, wl_keep)
    X <- kept$X; wl <- kept$wl
  }
  prep_abs <- function(X) -log10(pmax(X, .Machine$double.eps))
  prep_raw <- function(X) X
  prep_snv_abs <- function(X) prospectr::standardNormalVariate(prep_abs(X))
  prep_sg_abs <- function(X, w, p, m) prospectr::savitzkyGolay(prep_abs(X), m = m, p = p, w = w)
  prep_msc_abs <- function(X) prospectr::msc(prep_abs(X))
  prep_map <- list(
    RAW = function(X, ...) prep_raw(X),
    ABS = function(X, ...) prep_abs(X),
    SNV = function(X, ...) prep_snv_abs(X),
    SG0 = function(X, w = 11, p = 2, m = 0, ...) prep_sg_abs(X, w, p, m),
    SG1 = function(X, w = 11, p = 2, m = 1, ...) prep_sg_abs(X, w, p, m),
    SG2 = function(X, w = 15, p = 2, m = 2, ...) prep_sg_abs(X, w, p, m),
    MSC = function(X, ...) prep_msc_abs(X)
  )
  preprocess <- match.arg(preprocess, several.ok = TRUE, choices = names(prep_map))
  out_list <- lapply(preprocess, function(name) {
    fn <- prep_map[[name]]
    params <- preprocess_params[[name]]
    if (is.null(params)) params <- list()
    do.call(fn, c(list(X = X), params))
  })
  names(out_list) <- preprocess
  y <- if (!is.null(target) && target %in% names(data)) data[[target]] else NULL
  list(X = out_list, y = y, wl = wl)
}
