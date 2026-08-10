# Coverage for the S3 methods and the argument combinations that previously
# errored or returned silently wrong results.
#
# Tests that only check structure use a low maxit.optim: they do not depend on
# the optimiser converging, and a full fit costs roughly a minute.

QUICK <- 1000

test_that("summary and print work and return their input invisibly", {
  m <- story_mnTS()
  g <- story_glmm()

  expect_output(summary(m), "Call: mnTS")
  expect_output(summary(g), "Call: mnGLMM")
  expect_output(print(m), "Call: mnTS")
  expect_output(print(g), "Call: mnGLMM")

  # capture.output keeps the summaries out of the test log.
  invisible(capture.output({
    vm <- withVisible(print(m))
    vg <- withVisible(print(g))
  }))
  expect_false(vm$visible)
  expect_false(vg$visible)
  expect_identical(vm$value, m)
  expect_identical(vg$value, g)
})

test_that("print methods are registered for dispatch", {
  # print.mnTS used to exist without a NAMESPACE entry, so an installed package
  # fell through to print.default.
  expect_true("print.mnTS" %in% as.character(utils::methods(class = "mnTS")))
  expect_true("print.mnGLMM" %in% as.character(utils::methods(class = "mnGLMM")))
})

test_that("coef returns a coefficient table", {
  ct <- coef(story_mnTS())
  expect_true(is.matrix(ct))
  expect_equal(colnames(ct), c("Coef.", "se", "t", "P"))
  expect_gt(nrow(ct), 0)
})

# compute.information.matrix = FALSE -----------------------------------------

test_that("mnTS runs with compute.information.matrix = FALSE", {
  s <- story_setup()
  invisible(capture.output(
    m <- mnTS(Y = s$Y[s$Tsample, ], X = s$X, Tsample = s$Tsample,
              B0.fixed = s$B0.fixed, B.fixed = s$B.fixed,
              C.start = s$C.start, C.fixed = s$C.fixed,
              V.fixed = s$V.fixed, compute.information.matrix = FALSE,
              maxit.optim = QUICK)
  ))
  expect_s3_class(m, "mnTS")
  expect_null(m$se)
  expect_null(m$inv.information.matrix)
  # coef() and summary() must tolerate the missing information matrix.
  expect_equal(colnames(coef(m)), "Coef.")
  expect_output(summary(m), "Call: mnTS")
})

test_that("mnGLMM runs with compute.information.matrix = FALSE", {
  s <- story_setup()
  g <- mnGLMM(Y = s$Y[s$Tsample, ], X = s$X[s$Tsample, , drop = FALSE],
              B.start = s$B.start.glmm, B.fixed = s$B.fixed.glmm,
              V.fixed = s$V.fixed.glmm, compute.information.matrix = FALSE)
  expect_s3_class(g, "mnGLMM")
  expect_null(g$se)
  expect_output(summary(g), "Call: mnGLMM")
})

# no covariates ---------------------------------------------------------------

test_that("mnTS fits with X = NULL", {
  s <- story_setup()
  invisible(capture.output(
    m <- mnTS(Y = s$Y[s$Tsample, ], X = NULL, Tsample = s$Tsample,
              B0.fixed = s$B0.fixed, C.start = s$C.start, C.fixed = s$C.fixed,
              V.fixed = s$V.fixed, maxit.optim = QUICK)
  ))
  expect_s3_class(m, "mnTS")
  expect_null(m$X)
  expect_null(m$B)
  # B must not fall through to partial matching and pick up B0.
  expect_equal(dim(m$B0), c(1L, s$n))
  expect_output(summary(m), "Call: mnTS")
  expect_true(is.matrix(coef(m)))
  expect_true(is.matrix(simulate(m)$Y))
})

test_that("mnGLMM fits with X = NULL and with default B arguments", {
  s <- story_setup()
  g_null <- mnGLMM(Y = s$Y[s$Tsample, ], X = NULL, V.fixed = s$V.fixed.glmm)
  expect_s3_class(g_null, "mnGLMM")

  # The B defaults are built from X after it is rebuilt inside the function.
  # Relying on lazily-evaluated formals gave the wrong dimensions.
  g_default <- mnGLMM(Y = s$Y[s$Tsample, ], X = s$X[s$Tsample, , drop = FALSE],
                      V.fixed = s$V.fixed.glmm)
  expect_s3_class(g_default, "mnGLMM")
  expect_equal(g_default$logLik, story_glmm()$logLik, tolerance = 1e-6)
})

# simulate --------------------------------------------------------------------

test_that("simulate reuses the observed sampling effort", {
  m <- story_mnTS()
  s <- story_setup()
  sim <- simulate(m, seed = 1)

  # Sample sizes used to be hard-coded to 100 regardless of the data, which
  # biased every bootstrap replicate.
  expect_equal(as.numeric(rowSums(sim$Y[s$Tsample, ])),
               as.numeric(rowSums(m$Y)))
  # Time steps with no observation carry no counts.
  expect_true(all(is.na(sim$Y[-s$Tsample, ])))
})

test_that("simulate honours size, nsim and seed", {
  m <- story_mnTS()
  s <- story_setup()

  expect_equal(unique(rowSums(simulate(m, size = 400)$Y[s$Tsample, ])), 400)
  expect_error(simulate(m, size = c(1, 2)), "one entry per element")
  expect_error(simulate(m, size = 0), "at least 1")

  expect_length(simulate(m, nsim = 3), 3)
  expect_true(is.matrix(simulate(m, nsim = 1)$Y))
  expect_error(simulate(m, nsim = 0), "positive whole number")

  expect_identical(simulate(m, seed = 42)$Y, simulate(m, seed = 42)$Y)

  # seed must not leak into the caller's RNG stream.
  set.seed(11); before <- runif(1)
  set.seed(11); invisible(simulate(m, seed = 42)); after <- runif(1)
  expect_identical(before, after)
})

# bootstrap -------------------------------------------------------------------

test_that("bootstrap returns the documented structure and counts failures correctly", {
  m <- story_mnTS()
  set.seed(5)
  invisible(capture.output(b <- bootstrap(m, reps = 2, maxit.optim = QUICK)))

  expect_named(b, c("coef.table.boot", "all_mods_pars", "all_mods"))
  expect_true(is.matrix(b$coef.table.boot))
  expect_true(all(c("boot_mean", "se", "t", "P", "upper68", "lower68",
                    "n_fail_converged", "n_total") %in%
                  colnames(b$coef.table.boot)))
  expect_equal(unique(b$coef.table.boot[, "n_total"]), 2)

  # n_fail_converged counts replicates, not the sum of bobyqa's error codes.
  codes <- b$all_mods_pars[, "opt.convergence"]
  expect_equal(unique(b$coef.table.boot[, "n_fail_converged"]),
               sum(codes != 0))
  expect_lte(unique(b$coef.table.boot[, "n_fail_converged"]), 2)

  expect_error(bootstrap(m, reps = 0), "positive whole number")
})

test_that("the deprecated boot() still works but warns", {
  # boot() was renamed to bootstrap() because it clashed with boot::boot() from
  # the recommended boot package. The old name forwards for one release cycle.
  m <- story_mnTS()
  set.seed(5)
  invisible(capture.output(
    expect_warning(b <- boot(m, reps = 1, maxit.optim = QUICK), "deprecated")
  ))
  expect_named(b, c("coef.table.boot", "all_mods_pars", "all_mods"))

  # bootstrap() is the generic; boot() must not be one, or an unregistered
  # boot.mnTS would be treated as an apparent S3 method.
  expect_true("bootstrap.mnTS" %in% as.character(utils::methods(class = "mnTS")))
  expect_false(any(grepl("^boot\\.mnTS$", as.character(utils::methods(class = "mnTS")))))
})
