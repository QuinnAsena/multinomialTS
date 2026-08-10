# multinomialTS

`multinomialTS` fits state-space models to multinomially distributed
data, such as community count data or palaeoecological time series. It
estimates taxon-taxon interactions, driver-taxon relationships, and
cross-correlated error, and it handles time gaps in the response where
covariates are observed.

See [the vignettes on the package
website](https://quinnasena.github.io/multinomialTS/index.html) on how
to use the model.

## Quick start

``` r

library(multinomialTS)

data(story_pollen_matrix)
data(story_char_matrix)

# Rows of Y with counts; X must cover every time step the model iterates over.
sample_idx <- which(rowSums(story_pollen_matrix) != 0)
Y <- story_pollen_matrix[sample_idx, ]
X <- scale(story_char_matrix)
n <- ncol(Y)

# mnGLMM provides starting values for mnTS.
glmm_mod <- mnGLMM(Y = Y, X = X[sample_idx, , drop = FALSE], V.fixed = diag(n))

V.fixed <- matrix(NA, n, n); V.fixed[1] <- 1
B.fixed <- matrix(NA, ncol(X), n); B.fixed[, 1] <- 0
C.start <- 0.5 * diag(n); C.fixed <- C.start
C.fixed[C.fixed != 0] <- NA

ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
               B0.start = glmm_mod$B[1, , drop = FALSE],
               B.start  = glmm_mod$B[2, , drop = FALSE], B.fixed = B.fixed,
               C.start = C.start, C.fixed = C.fixed,
               V.start = glmm_mod$V, V.fixed = V.fixed,
               sigma.start = glmm_mod$sigma)
summary(ts_mod)

# Bootstrapping gives better standard errors than the Wald approximation.
# boot_out <- boot(ts_mod, reps = 1000)
```

## Installation

### Windows

The best way to install at the moment is to compile from source:

- If you have Rtools44 or (Rtools45 for R 4.5) and devtools installed,
  you can build the latest version of the package directly from github
  using:
  `devtools::install_github("https://github.com/QuinnAsena/multinomialTS")`

Alternatively, if you are having issues with the build tools, install
the binaries:

- `install.packages("https://github.com/QuinnAsena/multinomialTS/releases/download/v1.1.0/multinomialTS_1.1.0-windows.zip", repos = NULL, type = "win.binary")`

### macOS

The best way to install at the moment is to compile from source:

- If you have xcode-select and devtools installed, you can build the
  latest version of the package directly from github using:
  `devtools::install_github("https://github.com/QuinnAsena/multinomialTS")`

- On apple it is nice and easy to download xcode-select by opening a
  terminal and copying this code: `xcode-select --install`. Then try and
  run:
  `devtools::install_github("https://github.com/QuinnAsena/multinomialTS")`

Alternatively, if you are having issues with the build tools, install
the binaries:

For new M-series macs, use:

- `install.packages("https://github.com/QuinnAsena/multinomialTS/releases/download/v1.1.0/multinomialTS_1.1.0-macOS-arm64.tgz", repos = NULL)`

For Intel macs, use:

- `install.packages("https://github.com/QuinnAsena/multinomialTS/releases/download/v1.1.0/multinomialTS_1.1.0-macOS-x86_64.tgz", repos = NULL)`

### Linux

Compile from source. You need a C++ compiler (`g++` or `clang++`) plus
BLAS and LAPACK development headers, which most distributions package:

- Debian/Ubuntu:
  `sudo apt-get install build-essential libblas-dev liblapack-dev`
- Fedora/RHEL: `sudo dnf install gcc-c++ blas-devel lapack-devel`

Then:

- `devtools::install_github("https://github.com/QuinnAsena/multinomialTS")`

## Citation

`citation("multinomialTS")` gives the current reference.
