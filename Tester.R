devtools::document("spectrotune")
devtools::install("spectrotune")
library(spectrotune)
library(readxl)
library(spectrotune)

df <- readxl::read_excel("C:/Users/user/Desktop/df25.xlsx", sheet = "Overall")
str(df)

fit <- st_run(
  data = df,
  target = "Age",
  preprocess = c("ABS"),
  preprocess_params = list(SG1 = list(w = 17, p = 2, m = 1)),
  algos = c("PLSR","SVR"),
  outer_k = 5, outer_rep = 2, inner_k = 5, seed = 42,
  verbose = TRUE
)

fit$summary
