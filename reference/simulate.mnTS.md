# Simulate Data from a Multinomial State-Space Model

Generates simulated datasets from a fitted mnTS object.

## Usage

``` r
# S3 method for class 'mnTS'
simulate(object, nsim = 1, seed = NULL, size = NULL, ...)
```

## Arguments

- object:

  An object of class `"mnTS"`, as returned by the
  [`mnTS`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md)
  function.

- nsim:

  Number of simulated datasets to generate (default 1).

- seed:

  Optional seed for the random number generator. Passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html); the previous RNG
  state is restored on exit.

- size:

  Multinomial sample size (total counts) to draw at each sampled time
  step. Defaults to `NULL`, which reuses the observed totals from the
  data the model was fitted to, so simulated replicates carry the same
  sampling effort as the data. May also be a single number, or a vector
  with one entry per element of `object$Tsample`.

- ...:

  Additional arguments (currently unused).

## Value

When `nsim = 1`, a list with the simulated response `Y` and the model
components used to generate it (`X`, `B0`, `C`, `sigma`, `V`, `B`). When
`nsim > 1`, an unnamed list of `nsim` such lists. Rows of `Y` at time
steps not in `object$Tsample` are `NA`, matching the observation design
of the input data.

## Examples

``` r
# A small, fast example: the reference taxon plus two others over the first
# 60 time steps. See the package vignette for a full analysis.
data(story_pollen_matrix)
data(story_char_matrix)

Y_all <- story_pollen_matrix[, 1:3]
sample_idx <- which(rowSums(Y_all) != 0)
sample_idx <- sample_idx[sample_idx <= 60]
Y <- Y_all[sample_idx, ]
X <- scale(story_char_matrix)[1:60, , drop = FALSE]
n <- ncol(Y)

# NA = estimate, a number = hold fixed at that value.
B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
B.fixed <- matrix(NA, nrow = ncol(X), ncol = n)
B.fixed[, 1] <- 0                    # reference taxon
V.fixed <- matrix(NA, n, n)
V.fixed[1] <- 1                      # identifiability constraint
C.start <- 0.5 * diag(n)             # self-regulation only
C.fixed <- C.start
C.fixed[C.fixed != 0] <- NA

# \donttest{
ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
               B0.fixed = B0.fixed, B.fixed = B.fixed,
               C.start = C.start, C.fixed = C.fixed,
               V.fixed = V.fixed, dispersion.fixed = 1)
sim <- simulate(ts_mod, seed = 1)
# Simulated totals match the observed sampling effort.
rowSums(sim$Y[sample_idx, ])[1:5]
#> [1]  85  70 100  84 106
rowSums(ts_mod$Y)[1:5]
#> [1]  85  70 100  84 106
# }
```
