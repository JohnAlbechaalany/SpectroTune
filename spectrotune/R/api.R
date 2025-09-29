#' One-call convenience: preprocess + model
#' @param data data.frame with spectral columns named by wavelength and a target column
#' @param target character, name of the target column in data
#' @param preprocess character vector of pretreatments
#' @param preprocess_params named list of parameter lists for pretreatments
#' @param noise_range numeric length-2 to exclude (or NULL)
#' @param wl_keep single length-2 or list of ranges to keep (or NULL)
#' @param algos character vector among c("PLSR","SVR")
#' @param outer_k,outer_rep,inner_k integers for nested CV
#' @param seed integer
#' @param verbose logical, print progress
#' @return list(results, summary, best) from st_model
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
                   verbose = FALSE) {
  pp <- st_preprocess(
    data = data,
    target = target,
    preprocess = preprocess,
    preprocess_params = preprocess_params,
    noise_range = noise_range,
    wl_keep = wl_keep
  )
  st_model(
    X_list = pp$X,
    y = pp$y,
    algos = algos,
    outer_k = outer_k,
    outer_rep = outer_rep,
    inner_k = inner_k,
    seed = seed,
    verbose = verbose
  )
}
