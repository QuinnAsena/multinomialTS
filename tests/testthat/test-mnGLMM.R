# Baseline fit and input checks for mnGLMM.
# The shared setup lives in helper-fixtures.R.

fit_glmm <- function(..., s = story_setup()) {
  args <- list(Y = s$Y[s$Tsample, ], X = s$X[s$Tsample, , drop = FALSE],
               B.start = s$B.start.glmm, B.fixed = s$B.fixed.glmm,
               V.fixed = s$V.fixed.glmm)
  do.call(mnGLMM, utils::modifyList(args, list(...)))
}

test_that("mnGLMM works with example matrices", {
  g <- story_glmm()
  expect_s3_class(g, "mnGLMM")
  # Regression guard on the fitted likelihood, not just the class.
  expect_equal(g$logLik, 464.0598885, tolerance = 1e-6)
  expect_equal(dim(g$B), c(story_setup()$p, story_setup()$n))
  expect_true(all(is.finite(g$se)))
})

test_that("mnGLMM prints error for optimiser method", {
  expect_error(
    fit_glmm(method = "foo-bar"),
    regexp = "Acceptable methods are Nelder-Mead {optim},  BFGS {optim}, and bobyqa {minqa}.",
    fixed = TRUE)
})

test_that("mnGLMM prints error when Y has fewer rows than X", {
  s <- story_setup()
  Y <- s$Y[s$Tsample, ]
  expect_error(
    fit_glmm(Y = Y[1:(nrow(Y) - 10), ], s = s),
    regexp = "Variables X and Y have different lengths.")
})

test_that("mnGLMM prints error when X has fewer rows than Y", {
  s <- story_setup()
  expect_error(
    fit_glmm(X = scale(get("story_char_matrix")[1:50, , drop = FALSE]), s = s),
    regexp = "Variables X and Y have different lengths.")
})

test_that("mnGLMM prints error if V is the incorrect dimensions", {
  s <- story_setup()
  expect_error(
    fit_glmm(V.fixed = diag(s$n - 1), s = s),
    regexp = "V.fixed should have dimensions equal to the number of categories in Y.")
})

test_that("mnGLMM prints error if V is not symmetric", {
  s <- story_setup()
  V.fixed <- diag(s$n)
  V.fixed[1, 2] <- 0.5
  expect_error(
    fit_glmm(V.fixed = V.fixed, s = s),
    regexp = "V.fixed should be symmetric.")
})

test_that("mnGLMM prints error if B is the incorrect dimensions", {
  s <- story_setup()
  expect_error(
    fit_glmm(B.fixed = s$B.fixed.glmm[, 1:(s$n - 1)], s = s),
    regexp = "dimensions of B.fixed should be c(ncol(X)+1, ncol(Y)).",
    fixed = TRUE)
})
