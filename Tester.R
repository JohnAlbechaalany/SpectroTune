# 1) Clean load state
if ("spectrotune" %in% loadedNamespaces()) detach("package:spectrotune", unload=TRUE)

# 2) Regenerate docs and ensure exports
roxygen2::roxygenise("spectrotune")     # or devtools::document("spectrotune")

# 3) Reinstall and reload
devtools::install("spectrotune", force=TRUE, dependencies=TRUE)
library(spectrotune)

# 4) Verify functions are exported
ls("package:spectrotune")
getNamespaceExports("spectrotune")

install.packages("roxygen2") 
devtools::document("spectrotune")
devtools::install("spectrotune")
devtools::document("spectrotune")
help(package = "spectrotune")
?st_run
library(spectrotune)
?st_run   # or press F1 on st_run
library(spectrotune)
library(readxl)

fit <- st_run(
  data = df,
  target = "",
  preprocess = c("ABS"),
  preprocess_params = list(SG1 = list(w = 17, p = 2, m = 1)),
  algos = c("PLSR","SVR"),
  outer_k = 5, outer_rep = 2, inner_k = 5, seed = 42,
  verbose = TRUE
)

fit$summary
