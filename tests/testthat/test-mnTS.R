# Baseline: mnTS runs with example matrices
test_that("mnTS works with example matrices", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)
  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed <- diag(n)
  B.fixed <- matrix(c(rep(0,p), rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p), rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(
    Y = Y[Tsample, ], X = X[Tsample, , drop = FALSE],
    B.start = B.start, B.fixed = B.fixed, V.fixed = V.fixed
  )

  B0.start <- glmm_mod$B[1, , drop = FALSE]
  B.start <- glmm_mod$B[2, , drop = FALSE]
  sigma.start <- glmm_mod$sigma

  V.fixed <- matrix(NA, n, n); V.fixed[1] <- 1
  V.start <- glmm_mod$V

  B.fixed <- matrix(NA, ncol(X), n); B.fixed[,1] <- 0
  B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  C.start <- .5 * diag(n)
  C.fixed <- C.start; C.fixed[C.fixed != 0] <- NA

  ts_mod <- multinomialTS::mnTS(
    Y = Y[Tsample, ], X = X, Tsample = Tsample,
    B0.start = B0.start, B.start = B.start,
    C.start = C.start, C.fixed = C.fixed,
    B0.fixed = B0.fixed, B.fixed = B.fixed,
    V.fixed = V.fixed, V.start = V.start,
    dispersion.fixed = 1, maxit.optim = 1e5,
    sigma.start = sigma.start)

  expect_s3_class(ts_mod, "mnTS")
  expect_s3_class(glmm_mod, "mnGLMM")
})


# Wrong optimiser method
test_that("mnTS prints error for optimiser method", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)
  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed <- diag(n)
  B.fixed <- matrix(c(rep(0,p), rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p), rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(
    Y = Y[Tsample, ], X = X[Tsample, , drop = FALSE],
    B.start = B.start, B.fixed = B.fixed, V.fixed = V.fixed)

  B0.start <- glmm_mod$B[1, , drop = FALSE]
  B.start <- glmm_mod$B[2, , drop = FALSE]
  sigma.start <- glmm_mod$sigma

  V.fixed <- matrix(NA, n, n); V.fixed[1] <- 1
  V.start <- glmm_mod$V

  B.fixed <- matrix(NA, ncol(X), n); B.fixed[,1] <- 0
  B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  C.start <- .5 * diag(n)
  C.fixed <- C.start; C.fixed[C.fixed != 0] <- NA

  expect_error(
    ts_mod <- multinomialTS::mnTS(
      Y = Y[Tsample, ], X = X, Tsample = Tsample,
      B0.start = B0.start, B.start = B.start,
      C.start = C.start, C.fixed = C.fixed,
      B0.fixed = B0.fixed, B.fixed = B.fixed,
      V.fixed = V.fixed, V.start = V.start,
      dispersion.fixed = 1, maxit.optim = 1e5,
      sigma.start = sigma.start, method = "foo-bar"),
    regexp = "Acceptable methods are Nelder-Mead {optim},  BFGS {optim}, and bobyqa {minqa}.",
    fixed = TRUE
  )
})


# Incorrect direction of Y
test_that("mnTS prints error for incorrect direction of Y", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)
  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed <- diag(n)
  B.fixed <- matrix(c(rep(0,p), rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p), rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(
    Y = Y[Tsample, ], X = X[Tsample, , drop = FALSE],
    B.start = B.start, B.fixed = B.fixed, V.fixed = V.fixed
  )

  B0.start <- glmm_mod$B[1, , drop = FALSE]
  B.start <- glmm_mod$B[2, , drop = FALSE]
  sigma.start <- glmm_mod$sigma

  V.fixed <- matrix(NA, n, n); V.fixed[1] <- 1
  V.start <- glmm_mod$V

  B.fixed <- matrix(NA, ncol(X), n); B.fixed[,1] <- 0
  B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  C.start <- .5 * diag(n)
  C.fixed <- C.start; C.fixed[C.fixed != 0] <- NA

  expect_error(
    ts_mod <- multinomialTS::mnTS(
      Y = t(Y[Tsample, ]), X = X, Tsample = Tsample,
      B0.start = B0.start, B.start = B.start,
      C.start = C.start, C.fixed = C.fixed,
      B0.fixed = B0.fixed, B.fixed = B.fixed,
      V.fixed = V.fixed, V.start = V.start,
      dispersion.fixed = 1, maxit.optim = 1e5,
      sigma.start = sigma.start),
    regexp = "Time should run in the vertical direction."
  )
})


# Incorrect dimensions of V
test_that("mnTS prints error for incorrect dimensions of V", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)
  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed <- diag(n)
  B.fixed <- matrix(c(rep(0,p), rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p), rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(
    Y = Y[Tsample, ], X = X[Tsample, , drop = FALSE],
    B.start = B.start, B.fixed = B.fixed, V.fixed = V.fixed
  )

  B0.start <- glmm_mod$B[1, , drop = FALSE]
  B.start <- glmm_mod$B[2, , drop = FALSE]
  sigma.start <- glmm_mod$sigma

  V.fixed <- matrix(NA, n, n); V.fixed[1] <- 1
  V.fixed <- V.fixed[, 1:(ncol(V.fixed)-1)]
  V.start <- glmm_mod$V

  B.fixed <- matrix(NA, ncol(X), n); B.fixed[,1] <- 0
  B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  C.start <- .5 * diag(n)
  C.fixed <- C.start; C.fixed[C.fixed != 0] <- NA

  expect_error(
    ts_mod <- multinomialTS::mnTS(
      Y = Y[Tsample, ], X = X, Tsample = Tsample,
      B0.start = B0.start, B.start = B.start,
      C.start = C.start, C.fixed = C.fixed,
      B0.fixed = B0.fixed, B.fixed = B.fixed,
      V.fixed = V.fixed, V.start = V.start,
      dispersion.fixed = 1, maxit.optim = 1e5,
      sigma.start = sigma.start),
    regexp = "V.fixed should be diagonal with dimensions equal to the number of categories in Y."
  )
})


# Incorrect dimensions of V
test_that("mnTS prints error for incorrect dimensions of V", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)
  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed <- diag(n)
  B.fixed <- matrix(c(rep(0,p), rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p), rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(
    Y = Y[Tsample, ], X = X[Tsample, , drop = FALSE],
    B.start = B.start, B.fixed = B.fixed, V.fixed = V.fixed
  )

  B0.start <- glmm_mod$B[1, , drop = FALSE]
  B.start <- glmm_mod$B[2, , drop = FALSE]
  sigma.start <- glmm_mod$sigma

  V.fixed <- matrix(NA, n, n); V.fixed[1] <- 1
  V.fixed[1, 2] <- 0.5
  V.start <- glmm_mod$V

  B.fixed <- matrix(NA, ncol(X), n); B.fixed[,1] <- 0
  B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  C.start <- .5 * diag(n)
  C.fixed <- C.start; C.fixed[C.fixed != 0] <- NA

  expect_error(
    ts_mod <- multinomialTS::mnTS(
      Y = Y[Tsample, ], X = X, Tsample = Tsample,
      B0.start = B0.start, B.start = B.start,
      C.start = C.start, C.fixed = C.fixed,
      B0.fixed = B0.fixed, B.fixed = B.fixed,
      V.fixed = V.fixed, V.start = V.start,
      dispersion.fixed = 1, maxit.optim = 1e5,
      sigma.start = sigma.start),
    regexp = "V.fixed should be symmetric."
  )
})


# Sigma < 0
# Wrong optimiser method
test_that("mnTS prints error for sigma < 0", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)
  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed <- diag(n)
  B.fixed <- matrix(c(rep(0,p), rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p), rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(
    Y = Y[Tsample, ], X = X[Tsample, , drop = FALSE],
    B.start = B.start, B.fixed = B.fixed, V.fixed = V.fixed)

  B0.start <- glmm_mod$B[1, , drop = FALSE]
  B.start <- glmm_mod$B[2, , drop = FALSE]
  sigma.start <- -0.1

  V.fixed <- matrix(NA, n, n); V.fixed[1] <- 1
  V.start <- glmm_mod$V

  B.fixed <- matrix(NA, ncol(X), n); B.fixed[,1] <- 0
  B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  C.start <- .5 * diag(n)
  C.fixed <- C.start; C.fixed[C.fixed != 0] <- NA

  expect_error(
    ts_mod <- multinomialTS::mnTS(
      Y = Y[Tsample, ], X = X, Tsample = Tsample,
      B0.start = B0.start, B.start = B.start,
      C.start = C.start, C.fixed = C.fixed,
      B0.fixed = B0.fixed, B.fixed = B.fixed,
      V.fixed = V.fixed, V.start = V.start,
      dispersion.fixed = 1, maxit.optim = 1e5,
      sigma.start = sigma.start),
    regexp = "Sigma is less than 0"
  )
})

