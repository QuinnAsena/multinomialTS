# Shared fixtures for the test suite.
#
# Fitting mnGLMM/mnTS is slow, so the standard Story Lake setup is built once
# here and cached rather than repeated in every test_that() block.

# Matrices and the parameter templates used by nearly every test.
story_setup <- function() {
  Y <- get("story_pollen_matrix")
  X <- scale(get("story_char_matrix"))
  Tsample <- which(rowSums(Y) != 0)
  n <- ncol(Y)
  p <- ncol(X) + 1

  V.fixed.glmm <- diag(n)
  B.fixed.glmm <- matrix(c(rep(0, p), rep(NA, (n - 1) * p)), p, n)
  B.start.glmm <- matrix(c(rep(0, p), rep(0.01, (n - 1) * p)), p, n)

  V.fixed <- matrix(NA, n, n)
  V.fixed[1] <- 1
  B.fixed <- matrix(NA, ncol(X), n)
  B.fixed[, 1] <- 0
  B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
  C.start <- 0.5 * diag(n)
  C.fixed <- C.start
  C.fixed[C.fixed != 0] <- NA

  list(Y = Y, X = X, Tsample = Tsample, n = n, p = p,
       V.fixed.glmm = V.fixed.glmm, B.fixed.glmm = B.fixed.glmm,
       B.start.glmm = B.start.glmm,
       V.fixed = V.fixed, B.fixed = B.fixed, B0.fixed = B0.fixed,
       C.start = C.start, C.fixed = C.fixed)
}

# Cache the two fitted models so the whole suite pays for each fit only once.
.fit_cache <- new.env(parent = emptyenv())

story_glmm <- function() {
  if (is.null(.fit_cache$glmm)) {
    s <- story_setup()
    .fit_cache$glmm <- mnGLMM(
      Y = s$Y[s$Tsample, ], X = s$X[s$Tsample, , drop = FALSE],
      B.start = s$B.start.glmm, B.fixed = s$B.fixed.glmm,
      V.fixed = s$V.fixed.glmm)
  }
  .fit_cache$glmm
}

story_mnTS <- function() {
  if (is.null(.fit_cache$mnTS)) {
    s <- story_setup()
    g <- story_glmm()
    .fit_cache$mnTS <- mnTS(
      Y = s$Y[s$Tsample, ], X = s$X, Tsample = s$Tsample,
      B0.start = g$B[1, , drop = FALSE], B0.fixed = s$B0.fixed,
      B.start = g$B[2, , drop = FALSE], B.fixed = s$B.fixed,
      C.start = s$C.start, C.fixed = s$C.fixed,
      V.start = g$V, V.fixed = s$V.fixed,
      sigma.start = g$sigma,
      dispersion.fixed = 1, maxit.optim = 1e5)
  }
  .fit_cache$mnTS
}
