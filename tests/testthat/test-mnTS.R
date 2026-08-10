# Baseline fit and input checks for mnTS.
# The shared setup lives in helper-fixtures.R.

fit_mnTS <- function(..., s = story_setup()) {
  g <- story_glmm()
  args <- list(Y = s$Y[s$Tsample, ], X = s$X, Tsample = s$Tsample,
               B0.start = g$B[1, , drop = FALSE], B0.fixed = s$B0.fixed,
               B.start = g$B[2, , drop = FALSE], B.fixed = s$B.fixed,
               C.start = s$C.start, C.fixed = s$C.fixed,
               V.start = g$V, V.fixed = s$V.fixed,
               sigma.start = g$sigma,
               dispersion.fixed = 1, maxit.optim = 1e5)
  do.call(mnTS, utils::modifyList(args, list(...)))
}

test_that("mnTS works with example matrices", {
  s <- story_setup()
  m <- story_mnTS()
  expect_s3_class(m, "mnTS")
  expect_s3_class(story_glmm(), "mnGLMM")

  # Shapes the documentation promises.
  expect_equal(dim(m$B0), c(1L, s$n))
  expect_equal(dim(m$C), c(s$n, s$n))
  expect_equal(dim(m$V), c(s$n, s$n))
  expect_equal(nrow(m$mu), max(s$Tsample))
  expect_equal(nrow(m$y), max(s$Tsample))
  # Fitted probabilities must be a valid simplex at every time step.
  expect_equal(unname(rowSums(m$mu)), rep(1, nrow(m$mu)), tolerance = 1e-8)
  expect_equal(m$AIC, -2 * m$logLik + 2 * m$npar)
})

test_that("mnTS prints error for optimiser method", {
  expect_error(
    fit_mnTS(method = "foo-bar", maxit.optim = 200),
    regexp = "Acceptable methods are Nelder-Mead {optim},  BFGS {optim}, and bobyqa {minqa}.",
    fixed = TRUE)
})

test_that("mnTS prints error for incorrect direction of Y", {
  s <- story_setup()
  expect_error(
    fit_mnTS(Y = t(s$Y[s$Tsample, ]), maxit.optim = 200, s = s),
    regexp = "Time should run in the vertical direction.")
})

test_that("mnTS prints error for incorrect dimensions of V", {
  s <- story_setup()
  V.fixed <- s$V.fixed[, 1:(ncol(s$V.fixed) - 1)]
  expect_error(
    fit_mnTS(V.fixed = V.fixed, maxit.optim = 200, s = s),
    regexp = "V.fixed should be diagonal with dimensions equal to the number of categories in Y.")
})

test_that("mnTS prints error if V is not symmetric", {
  s <- story_setup()
  V.fixed <- s$V.fixed
  V.fixed[1, 2] <- 0.5
  expect_error(
    fit_mnTS(V.fixed = V.fixed, maxit.optim = 200, s = s),
    regexp = "V.fixed should be symmetric.")
})

test_that("mnTS prints error for sigma < 0", {
  expect_error(
    fit_mnTS(sigma.start = -0.1, maxit.optim = 200),
    regexp = "Sigma is less than 0")
})
