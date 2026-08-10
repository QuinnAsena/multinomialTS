# Input validation. Each of these used to either fail with an opaque error from
# Armadillo, or worse, be accepted silently.

test_that("mnTS rejects an unsorted or mis-sized Tsample", {
  s <- story_setup()
  fit <- function(Y, Tsample) {
    mnTS(Y = Y, X = s$X, Tsample = Tsample,
         B0.fixed = s$B0.fixed, B.fixed = s$B.fixed,
         C.start = s$C.start, C.fixed = s$C.fixed,
         V.fixed = s$V.fixed, maxit.optim = 200)
  }

  # Reversed Tsample used to return a model fitted over a single time step.
  rev_T <- rev(s$Tsample)
  expect_error(fit(s$Y[rev_T, ], rev_T), "strictly increasing")

  expect_error(fit(s$Y[s$Tsample, ], s$Tsample[-1]),
               "must equal nrow\\(Y\\)")
  expect_error(fit(s$Y[s$Tsample, ], s$Tsample + 0.5),
               "whole-number row indices")
  expect_error(fit(s$Y[s$Tsample, ], replace(s$Tsample, 1, -1)),
               "whole-number row indices")
  expect_error(fit(s$Y[s$Tsample, ], replace(s$Tsample, 1, NA)),
               "no missing values")
})

test_that("mnTS rejects an X that does not cover every time step", {
  s <- story_setup()
  expect_error(
    mnTS(Y = s$Y[s$Tsample, ], X = s$X[1:50, , drop = FALSE],
         Tsample = s$Tsample, B0.fixed = s$B0.fixed, B.fixed = s$B.fixed,
         C.start = s$C.start, C.fixed = s$C.fixed, V.fixed = s$V.fixed,
         maxit.optim = 200),
    "must cover every time step")
})

test_that("both fitters reject a non-positive dispersion", {
  s <- story_setup()

  # dispersion is log-transformed internally, so 0 produced -Inf in par.fixed
  # and surfaced as "Mat::elem(): size mismatch" from Armadillo.
  expect_error(
    mnTS(Y = s$Y[s$Tsample, ], X = s$X, Tsample = s$Tsample,
         B0.fixed = s$B0.fixed, B.fixed = s$B.fixed,
         C.start = s$C.start, C.fixed = s$C.fixed, V.fixed = s$V.fixed,
         dispersion.fixed = 0, maxit.optim = 200),
    "dispersion.fixed must be a positive, finite value")

  expect_error(
    mnGLMM(Y = s$Y[s$Tsample, ], X = s$X[s$Tsample, , drop = FALSE],
           B.start = s$B.start.glmm, B.fixed = s$B.fixed.glmm,
           V.fixed = s$V.fixed.glmm, dispersion.fixed = -1),
    "dispersion.fixed must be a positive, finite value")

  expect_error(
    mnGLMM(Y = s$Y[s$Tsample, ], X = s$X[s$Tsample, , drop = FALSE],
           B.start = s$B.start.glmm, B.fixed = s$B.fixed.glmm,
           V.fixed = s$V.fixed.glmm, dispersion.start = 0),
    "dispersion.start must be a positive, finite value")
})

test_that("REML has been removed from mnGLMM", {
  s <- story_setup()
  # REML was documented but never reached the optimiser, and its implementation
  # errored when called directly. Passing it must now fail loudly rather than
  # silently returning an ML fit.
  expect_error(
    mnGLMM(Y = s$Y[s$Tsample, ], X = s$X[s$Tsample, , drop = FALSE],
           B.start = s$B.start.glmm, B.fixed = s$B.fixed.glmm,
           V.fixed = s$V.fixed.glmm, REML = TRUE),
    "unused argument")
})
