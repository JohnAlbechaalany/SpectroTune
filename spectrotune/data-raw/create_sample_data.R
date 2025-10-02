# Generate sample spectral data for SpectroTune package
# This script creates the sample_spectra dataset

set.seed(42)

# Create wavelength vector (400-2499 nm, 5nm intervals)
wavelengths <- seq(400, 2499, by = 5)
n_wl <- length(wavelengths)

# Generate 150 samples
n_samples <- 150

# Create sample IDs
sample_ids <- paste0("Sample_", sprintf("%03d", 1:n_samples))

# Generate synthetic spectra with realistic characteristics
spectra <- matrix(0, nrow = n_samples, ncol = n_wl)

for (i in 1:n_samples) {
  # Base spectrum with water absorption bands
  base_spectrum <- 0.3 + 0.4 * exp(-((wavelengths - 1000) / 300)^2) +
                   0.2 * exp(-((wavelengths - 1500) / 200)^2) +
                   0.1 * exp(-((wavelengths - 2000) / 150)^2)
  
  # Add noise
  noise <- rnorm(n_wl, 0, 0.02)
  
  # Add baseline drift
  baseline <- 0.05 * (wavelengths - 400) / (2499 - 400)
  
  # Combine
  spectra[i, ] <- pmax(0, pmin(1, base_spectrum + noise + baseline))
}

##
## Generate target variables with real signal in the spectra
## --------------------------------------------------------
# Latent factors that drive both spectra and targets
z1 <- rnorm(n_samples)
z2 <- rnorm(n_samples)

# Build informative spectral bands (two Gaussian peaks)
peak <- function(w0, bw) exp(-((wavelengths - w0) / bw)^2)
band1 <- peak(1100, 120)    # CH/OH band
band2 <- peak(1650, 140)    # protein band

# Inject signal from latent factors into spectra
for (i in 1:n_samples) {
  spectra[i, ] <- spectra[i, ] +
    0.18 * z1[i] * band1 +
    0.14 * z2[i] * band2
}
# keep matrix shape when clipping values
spectra <- pmax(spectra, 0)
spectra <- pmin(spectra, 1)

# Age (regression target) linked to latent factors
age <- 45 + 12 * z1 + 8 * z2 + rnorm(n_samples, sd = 5)
age <- round(pmax(20, pmin(80, age)))

# Class (classification target) derived from a linear score of factors
score <- 0.9 * z1 - 0.7 * z2 + rnorm(n_samples, sd = 0.3)
cuts <- quantile(score, probs = c(0.33, 0.66))
class <- ifelse(score <= cuts[1], "C", ifelse(score <= cuts[2], "B", "A"))
class_labels <- c("A","B","C")
class <- factor(class, levels = class_labels)

# Covariates for stratification
breed <- sample(c("Breed1", "Breed2", "Breed3"), n_samples, replace = TRUE, 
                prob = c(0.5, 0.3, 0.2))
breed <- factor(breed)

sex <- sample(c("M", "F"), n_samples, replace = TRUE, prob = c(0.6, 0.4))
sex <- factor(sex)

# Create data frame
colnames(spectra) <- as.character(wavelengths)

sample_spectra <- data.frame(
  ID = sample_ids,
  Age = age,
  Class = class,
  breed = breed,
  sex = sex,
  spectra,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Save as package data (no usethis dependency)
pkg_data_dir <- file.path(getwd(), "data")
if (!dir.exists(pkg_data_dir)) dir.create(pkg_data_dir, recursive = TRUE)
save(sample_spectra, file = file.path(pkg_data_dir, "sample_spectra.rda"))

cat("Sample spectral data generated and saved to", file.path(pkg_data_dir, "sample_spectra.rda"), "\n")
cat("Dataset contains", n_samples, "samples and", n_wl, "wavelengths\n")
