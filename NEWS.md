# multinomialTS 1.0.3.9000

Unreleased development version. Given the breaking changes below, consider
releasing this as 1.1.0.

## Major changes

* `simulate.mnTS()` now draws each replicate using the **observed** multinomial
  totals rather than a hard-coded size of 100. Because `boot.mnTS()` is built on
  `simulate()`, every bootstrap replicate previously carried far more
  observation noise than the real data, which inflated bootstrapped standard
  errors and attenuated coefficients. **Bootstrap results produced with earlier
  versions should be re-run.** A `size` argument is available to override the
  default.

* `simulate.mnTS()` gains the `nsim` and `seed` arguments of the
  `stats::simulate()` generic. Both were previously accepted and silently
  ignored. With `nsim > 1` the return value is a list of `nsim` simulations.

* Rows of the simulated `Y` at time steps not in `Tsample` are now `NA` rather
  than being filled with invented counts, matching the observation design of the
  input data.

* The `REML` argument of `mnGLMM()` has been removed. It was documented but
  never reached the optimiser, so `REML = TRUE` silently returned a maximum
  likelihood fit; the underlying implementation also errored when called
  directly. Passing `REML` now raises an "unused argument" error.

* `sunfish_pollen_wide`'s truncated `P.strobu` column has been renamed to
  `P.strobus`.

## Bug fixes

* `mnTS()` now rejects a `Tsample` that is not strictly increasing. An
  out-of-order `Tsample` was previously accepted and returned a model fitted
  over the wrong time span, with no warning. `Tsample` is also checked for
  whole-number positive indices and for having one entry per row of `Y`.

* `compute.information.matrix = FALSE` now works in both `mnTS()` and
  `mnGLMM()`. Both errored with "object not found" because the branch failed to
  define values used when assembling the result.

* `X = NULL` now works in `mnTS()`. It is documented as supported, but failed
  with "Not a matrix." because RcppArmadillo cannot coerce `NULL` to a matrix.
  `summary()`, `coef()`, `simulate()` and `boot()` all handle no-covariate
  models.

* The default `B.fixed` and `B.start` of `mnGLMM()` were computed from `X` after
  `X` had been rebuilt inside the function, so they had the wrong dimensions and
  the `X = NULL` branch was unreachable. They are now built correctly.

* `mnTS()` and `mnGLMM()` reject a non-positive or non-finite
  `dispersion.fixed`/`dispersion.start`. Because dispersion is log-transformed
  internally, `dispersion.fixed = 0` produced `-Inf` and surfaced as an opaque
  "Mat::elem(): size mismatch" from Armadillo.

* `mnTS()` checks that `X` has at least `max(Tsample)` rows instead of failing
  with "Mat::col(): index out of bounds".

* `coef.mnTS()` no longer errors when the model was fitted without an
  information matrix.

* `boot.mnTS()`'s `n_fail_converged` counted the *sum of* `bobyqa`'s error
  codes, so one replicate failing with `ierr = 5` was reported as five
  failures. It now counts replicates.

* `boot.mnTS()` handles models fitted without covariates, and identifies
  successful replicates with `inherits()` instead of a list-to-string
  comparison that broke for multi-class objects.

* `print()` now dispatches for `mnTS` and `mnGLMM` objects. The methods existed
  but were absent from `NAMESPACE`, so an installed package fell through to
  `print.default`.

* The C++ likelihoods identified free parameters with `find_nonfinite`
  (NaN *and* infinities) while R used `is.na()` (NaN only), so a non-finite
  fixed parameter silently misaligned the parameter vector.

## Documentation

* All exported functions now have runnable examples.

* Corrected the documented format of `story_100_bootstraps` (a matrix, not a
  list), `sunfish_ll` (4 columns, not 2), and `story_char_matrix` (the values
  are not pre-scaled).

* `story_pollen_wide` and `sunfish_pollen_wide` no longer ship as grouped
  tibbles.

## Internal

* `mnGLMM()` computes its Hessian with the compiled likelihood, as `mnTS()`
  already did, and exposes `hessian.method.args`.

* The pure-R likelihoods have moved out of `R/` into the test suite, where they
  serve as the reference implementations the C++ is checked against.

* Added an `R-CMD-check` workflow. The macOS Intel binary is now built on an
  Intel runner and its architecture is verified; both macOS binaries previously
  came from the same Apple Silicon runner.
