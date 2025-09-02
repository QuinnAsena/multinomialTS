# Testing mnGLMM function runs with correct inputs
test_that("mnGLMM works with example matrices", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)

  story_pollen_matrix[Tsample, ]

  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed = diag(n)
  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(
    Y = Y[Tsample, ], X = X[Tsample, ,drop = F],
    B.start = B.start, B.fixed = B.fixed,
    V.fixed = V.fixed)

  expect_s3_class(glmm_mod, "mnGLMM")
})

# Testing mnGLMM function runs with incorrect method for optimiser
test_that("mnGLMM prints error for optimiser method", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)

  story_pollen_matrix[Tsample, ]

  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed = diag(n)
  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  expect_error(
    glmm_mod <- multinomialTS::mnGLMM(
      Y = Y[Tsample, ], X = X[Tsample, ,drop = F],
      B.start = B.start, B.fixed = B.fixed,
      V.fixed = V.fixed, method = "foo-bar"),
  regexp = "Acceptable methods are Nelder-Mead {optim},  BFGS {optim}, and bobyqa {minqa}.",
  fixed = TRUE
  )
})


# Testing mnGLMM function runs with incorrect dimensions of Y
test_that("mnGLMM prints error for incorrect dimensions of X and Y", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)

  story_pollen_matrix[Tsample, ]

  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix[Tsample, ]
  Y <- Y[1:(nrow(Y)-10), ]
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed = diag(n)
  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  expect_error(
    glmm_mod <- multinomialTS::mnGLMM(
      Y = Y, X = X[Tsample, ,drop = F],
      B.start = B.start, B.fixed = B.fixed,
      V.fixed = V.fixed),
    regexp = "Variables X and Y have different lengths."
  )
})


# Testing mnGLMM function runs with incorrect dimensions of X
test_that("mnGLMM prints error for incorrect dimensions of X and Y", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)

  story_pollen_matrix[Tsample, ]

  X <- scale(story_char_matrix[1:(nrow(story_char_matrix)-10), ])
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed = diag(n)
  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  expect_error(
    glmm_mod <- multinomialTS::mnGLMM(
      Y = Y[Tsample, ], X = X,
      B.start = B.start, B.fixed = B.fixed,
      V.fixed = V.fixed),
    regexp = "Variables X and Y have different lengths."
  )
})


# Testing mnGLMM function runs with incorrect dimensions of V
test_that("mnGLMM prints error if V is the incorrect dimensions", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)

  story_pollen_matrix[Tsample, ]

  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed = diag(n-1)
  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  expect_error(
    glmm_mod <- multinomialTS::mnGLMM(
      Y = Y[Tsample, ], X = X[Tsample, ,drop = F],
      B.start = B.start, B.fixed = B.fixed,
      V.fixed = V.fixed),
    regexp = "V.fixed should have dimensions equal to the number of categories in Y."
  )
})

# Testing mnGLMM function runs with incorrect dimensions of V
test_that("mnGLMM prints error if V is not symetric", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)

  story_pollen_matrix[Tsample, ]

  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed = diag(n)
  V.fixed[1,2] <- 0.5
  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  expect_error(
    glmm_mod <- multinomialTS::mnGLMM(
      Y = Y[Tsample, ], X = X[Tsample, ,drop = F],
      B.start = B.start, B.fixed = B.fixed,
      V.fixed = V.fixed),
    regexp = "V.fixed should be symmetric."
  )
})

# Testing mnGLMM function runs with incorrect dimensions of B
test_that("mnGLMM prints error if B is the incorrect dimensions", {
  data("story_char_matrix", package = "multinomialTS")
  data("story_pollen_matrix", package = "multinomialTS")

  Tsample <- which(rowSums(story_pollen_matrix) != 0)

  story_pollen_matrix[Tsample, ]

  X <- scale(story_char_matrix)
  Y <- story_pollen_matrix
  p <- ncol(X) + 1
  n <- ncol(Y)

  V.fixed = diag(n)
  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.fixed <- B.fixed[ ,1:(ncol(B.fixed)-1)]
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  expect_error(
    glmm_mod <- multinomialTS::mnGLMM(
      Y = Y[Tsample, ], X = X[Tsample, ,drop = F],
      B.start = B.start, B.fixed = B.fixed,
      V.fixed = V.fixed),
    regexp = "dimensions of B.fixed should be c(ncol(X)+1, ncol(Y)).",
    fixed = TRUE
  )
})
