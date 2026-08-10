# The likelihoods are evaluated in C++ at fit time. These tests check the C++
# against the pure-R reference implementations in helper-reference-impl.R.
# Both are deterministic, so they must agree to machine precision.

tol <- 1e-10

# mnTS ------------------------------------------------------------------------

mnTS_pars <- function(n) {
  # Layout: B0(n) | vec(C)(n*n) | sigma | vec(L)(n*n) | log(dispersion)
  c(rep(-1, n),                  # B0
    as.numeric(diag(0.5, n)),    # C, spectral radius 0.5 (stationary)
    0.5,                         # sigma
    as.numeric(diag(1, n)),      # L
    log(1))                      # log dispersion
}

# mnTS.ml() takes Y as n x Tmax, so transpose the stored matrix.
mnTS_Y <- function() t(as.matrix(get("story_pollen_wide")[, 3:7]))

expect_cpp_matches_r <- function(par, par_fixed, Y, X, Tsample) {
  r_val <- mnTS.ml(par, par.fixed = par_fixed, Y = Y, X = X, Tsample = Tsample)
  cpp_val <- mnTS_ml_cpp_listout(par, par_fixed, Y, X, Tsample)$LL
  expect_equal(as.numeric(cpp_val), as.numeric(r_val), tolerance = tol)
}

test_that("mnTS C++ matches the R reference: stationary C, no covariates", {
  Y <- mnTS_Y(); n <- nrow(Y)
  par <- mnTS_pars(n)
  expect_cpp_matches_r(par, rep(NA_real_, length(par)), Y, NULL, seq_len(ncol(Y)))
})

test_that("mnTS C++ matches the R reference: non-stationary C", {
  Y <- mnTS_Y(); n <- nrow(Y)
  par_all <- mnTS_pars(n)
  par_fixed <- rep(NA_real_, length(par_all))
  # Fix C to a diagonal with spectral radius 1.5; the filter takes the
  # non-stationary initialisation branch (PP <- Se).
  par_fixed[n + seq_len(n * n)] <- as.numeric(diag(1.5, n))
  expect_cpp_matches_r(par_all[is.na(par_fixed)], par_fixed, Y, NULL, seq_len(ncol(Y)))
})

test_that("mnTS C++ matches the R reference: with a covariate", {
  Y <- mnTS_Y(); n <- nrow(Y); Tmax <- ncol(Y)
  X <- matrix(scale(seq_len(Tmax)), nrow = 1, ncol = Tmax)   # p x Tmax
  par <- c(mnTS_pars(n), rep(0, n * nrow(X)))
  expect_cpp_matches_r(par, rep(NA_real_, length(par)), Y, X, seq_len(Tmax))
})

test_that("mnTS C++ matches the R reference: sparse Tsample", {
  Y <- mnTS_Y(); n <- nrow(Y)
  Tsample <- seq(1, ncol(Y), by = 3)
  par <- mnTS_pars(n)
  expect_cpp_matches_r(par, rep(NA_real_, length(par)), Y[, Tsample], NULL, Tsample)
})

test_that("mnTS C++ matches the R reference: NA column in Y is skipped", {
  Y <- mnTS_Y(); n <- nrow(Y)
  Y[, 10] <- NA_real_
  par <- mnTS_pars(n)
  expect_cpp_matches_r(par, rep(NA_real_, length(par)), Y, NULL, seq_len(ncol(Y)))
})

# mnGLMM ----------------------------------------------------------------------

mnGLMM_Y <- function() as.matrix(get("story_pollen_wide")[, 3:7])

mnGLMM_pars <- function(n, p) {
  # Layout: B(n*p) | L(n*n) | sigma | log(dispersion)
  c(rep(0, n * p), as.numeric(diag(1, n)), 0.1, log(1))
}

expect_glmm_cpp_matches_r <- function(par, par_fixed, Y, X) {
  r_val <- mnGLMM.ml(par, par.fixed = par_fixed, Y = Y, X = X, fitted.values = FALSE)
  cpp_val <- mnGLMM_ml_cpp_listout(par, par_fixed, Y, X)$LL
  expect_equal(as.numeric(cpp_val), as.numeric(r_val), tolerance = tol)
}

test_that("mnGLMM C++ matches the R reference: intercept only", {
  Y <- mnGLMM_Y(); n <- ncol(Y)
  X <- matrix(1, nrow = nrow(Y), ncol = 1, dimnames = list(NULL, "(intercept)"))
  par <- mnGLMM_pars(n, ncol(X))
  expect_glmm_cpp_matches_r(par, rep(NA_real_, length(par)), Y, X)
})

test_that("mnGLMM C++ matches the R reference: intercept plus covariate", {
  Y <- mnGLMM_Y(); n <- ncol(Y)
  X <- cbind(1, scale(seq_len(nrow(Y))))
  colnames(X) <- c("(intercept)", "time")
  par <- mnGLMM_pars(n, ncol(X))
  expect_glmm_cpp_matches_r(par, rep(NA_real_, length(par)), Y, X)
})

test_that("mnGLMM C++ matches the R reference: sigma or dispersion fixed", {
  Y <- mnGLMM_Y(); n <- ncol(Y)
  X <- matrix(1, nrow = nrow(Y), ncol = 1, dimnames = list(NULL, "(intercept)"))
  par_all <- mnGLMM_pars(n, ncol(X))

  par_fixed <- rep(NA_real_, length(par_all))
  par_fixed[n * ncol(X) + n * n + 1] <- 0.1                 # sigma fixed
  expect_glmm_cpp_matches_r(par_all[is.na(par_fixed)], par_fixed, Y, X)

  par_fixed <- rep(NA_real_, length(par_all))
  par_fixed[n * ncol(X) + n * n + 2] <- log(1)              # dispersion fixed
  expect_glmm_cpp_matches_r(par_all[is.na(par_fixed)], par_fixed, Y, X)
})

test_that("mnGLMM C++ matches the R reference: random parameter values", {
  Y <- mnGLMM_Y(); n <- ncol(Y)
  X <- matrix(1, nrow = nrow(Y), ncol = 1, dimnames = list(NULL, "(intercept)"))
  set.seed(99)
  par <- c(rnorm(n * ncol(X), 0, 0.5),
           as.numeric(diag(runif(n, 0.5, 1.5))),
           abs(rnorm(1, 0.2, 0.1)),
           log(runif(1, 0.5, 2)))
  expect_glmm_cpp_matches_r(par, rep(NA_real_, length(par)), Y, X)
})

test_that("mnGLMM C++ matches the R reference on the fitted-value outputs", {
  Y <- mnGLMM_Y(); n <- ncol(Y)
  X <- matrix(1, nrow = nrow(Y), ncol = 1, dimnames = list(NULL, "(intercept)"))
  par <- mnGLMM_pars(n, ncol(X))
  par_fixed <- rep(NA_real_, length(par))

  r_fit <- mnGLMM.ml(par, par.fixed = par_fixed, Y = Y, X = X, fitted.values = TRUE)
  cpp_fit <- mnGLMM_ml_cpp_listout(par, par_fixed, Y, X)

  expect_equal(unname(cpp_fit$B), unname(r_fit$B), tolerance = tol)
  expect_equal(unname(cpp_fit$mu), unname(r_fit$mu), tolerance = tol)
  expect_equal(as.numeric(cpp_fit$logLik), as.numeric(r_fit$logLik), tolerance = tol)
})
