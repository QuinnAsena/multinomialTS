# Open items

Maintainer notes: known issues and deferred work, carried forward from
the package audit. Not installed with the package (`.Rbuildignore`d).

Each item records **why it was deferred**, because in most cases the
reason is that the change would alter reported numbers or produce a diff
that is hard to review — not that the work is unclear.

Status at the time of writing: `R CMD check` gives 0 errors, 0 warnings
and one NOTE (`unable to verify current time`, an artefact of the
machine it was run on, not the package). Test suite: 114 assertions
passing.

------------------------------------------------------------------------

## Statistical / numerical

### 1. `mnTS()` and `mnGLMM()` use different standard-error conventions

[`mnTS()`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md)
fixes `dispersion` at its estimate and computes the Hessian over the
remaining parameters, giving **conditional** standard errors
(`R/source_mnTS.R`, the `compute.information.matrix` block).
[`mnGLMM()`](https://quinnasena.github.io/multinomialTS/reference/mnGLMM.md)
computes the full Hessian including log-dispersion and inverts it,
giving **marginal** standard errors (`R/source_mnGLMM.R`).

Inverting a submatrix of the Hessian and taking a submatrix of the
inverse are not the same thing, so the two functions report standard
errors that are not directly comparable when dispersion is estimated.
When dispersion is fixed — the default — the question does not arise,
because it is not a free parameter.

**Deferred because** aligning them changes reported standard errors for
any model that estimates dispersion. Needs a decision on which
convention is intended, ideally checked against the bootstrap (which is
the recommended way to get standard errors anyway).

### 2. `par["sigma"]` can disagree in sign with `sigma`

`sigma` enters the likelihood squared, so the optimiser can settle on a
negative value. The C++ applies
[`abs()`](https://rdrr.io/r/base/MathFun.html) to the `sigma` field it
returns but not to the `sigma` element of the `par` vector
(`src/source_mnTS_list.cpp`, `src/source_mnGLMM.cpp`), so
`object$par["sigma"]` and `object$sigma` can differ in sign.

The pure-R reference implementations in
`tests/testthat/helper-reference-impl.R` *do* apply
[`abs()`](https://rdrr.io/r/base/MathFun.html) there, so this is also a
genuine divergence between the reference and the compiled code. The
equivalency tests compare the likelihood and the `B`/`mu`/`logLik`
outputs, not `par`, so it is not currently exercised.

[`summary()`](https://rdrr.io/r/base/summary.html) and
[`coef()`](https://rdrr.io/r/stats/coef.html) exclude `sigma` from their
tables, so the practical impact is limited to anyone reading
`object$par` directly.

**Deferred because** it is cosmetic in normal use, and changing the C++
output should come with a test that actually drives sigma negative.

### 3. `mnGLMM()` does not validate `sigma.start`

[`mnTS()`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md)
rejects `sigma.start < 0` with “Sigma is less than 0”
(`R/source_mnTS.R`).
[`mnGLMM()`](https://quinnasena.github.io/multinomialTS/reference/mnGLMM.md)
has no equivalent check even though it uses `sigma` the same way.

**Deferred because** adding a new
[`stop()`](https://rdrr.io/r/base/stop.html) is a small behaviour change
that was outside the scope of the round it was found in. Low risk to
add.

### 4. Use `inv_sympd` for the symmetric positive-definite inversions

`arma::inv(FF)` in `src/source_mnTS_list.cpp` and `arma::inv(S_sub)` in
`src/source_mnGLMM.cpp` both invert matrices that are symmetric and
positive definite by construction. `arma::inv_sympd` is faster and more
numerically stable.

**Deferred because** it is a different algorithm and would perturb
results in the last few digits, which the equivalency tests assert to
machine precision. Doing it means deciding on a tolerance for those
tests.

### 5. No guard on `mu` reaching 0 or 1

`arma::pow(mu % (1 - mu), -1)` (`src/source_mnTS_list.cpp`,
`src/source_mnGLMM.cpp`) produces `Inf` if any fitted probability hits 0
or 1. The `rcond` check downstream catches the consequence and returns
the `1e+10` guard value, so the optimiser recovers, but a direct check
would give a clearer diagnostic.

------------------------------------------------------------------------

## Code quality

### 6. Run `styler::style_pkg()`

`R/` mixes tab indentation (the original code) with spaces (the newer
methods), sometimes inside a single function.

**Deferred because** a whole-package restyle produces a large
whitespace-only diff. It should be **its own commit**, with nothing else
in it, so that subsequent `git blame` and code review stay usable. Run
it, commit it alone, then confirm the test suite and `R CMD check` are
unchanged.

### 7. Vignette build cost

`vignettes/multinomialTS-vignette.qmd` takes about **30 minutes** to
render: it fits `mnGLMM`, two `mnTS` models, and a 5-replicate bootstrap
at build time. That cost is paid on every `R CMD build` and every
pkgdown deploy.

Options: reduce `maxit.optim` in the vignette, fit a smaller subset, or
precompute the fits and ship them as internal objects. Note that
shipping precomputed *bootstrap* output is what created the stale-data
problem fixed in 1.1.0 — see `data-raw/README.md` — so any
precomputation needs a regeneration script under `data-raw/`.

### 8. `boot()` deprecation shim

[`boot()`](https://quinnasena.github.io/multinomialTS/reference/boot.md)
forwards to
[`bootstrap()`](https://quinnasena.github.io/multinomialTS/reference/bootstrap.md)
with a deprecation warning (`R/source_mnTS.R`). Remove it in a future
release, along with its `man/boot.Rd` and the `Deprecated` section of
`_pkgdown.yml`.

------------------------------------------------------------------------

## Packaging / infrastructure

### 9. `build.yml`: validated for Windows and Apple Silicon, `release` job still untested

The binary-build workflow was substantially rewritten for 1.1.0
([`pkgbuild::build()`](https://pkgbuild.r-lib.org/reference/build.html)
instead of `devtools::build()`, dependencies via `setup-r-dependencies`,
a `lipo` architecture check, a tag guard on the `release` job).

A manual `workflow_dispatch` run confirmed the **windows** and
**macOS-arm64** legs build and upload artifacts correctly. The Intel leg
has since been removed (see below).

**Still unverified: the `release` job.** It is guarded to tag pushes, so
a manual dispatch never exercises it. The first `v*.*.*` tag push will
be its first run. Failure is benign — `needs: build` means a failed
build blocks publication entirely, so a partial release is not possible
— but be ready to delete the tag and release and retry:

``` sh
git push --delete origin v1.1.0 && git tag -d v1.1.0
```

#### Intel macOS: dropped, and why

`macos-13` was the last free Intel runner and **has been retired**. Jobs
requesting a retired label queue indefinitely rather than failing, which
is what happened: the Intel leg sat “Waiting for a runner” for 18 hours.

The remaining Intel labels (`macos-15-intel`, `macos-26-intel`,
`macos-*-large`) are *larger runners*, which GitHub bills even for
public repositories — “Larger runners are always charged for, even when
used by public repositories or when you have quota available from your
plan.” R packages with compiled code cannot be cross-compiled, so an
Intel binary needs a native Intel machine.

The Intel leg was therefore removed rather than paid for. Reinstating it
means either enabling billing and switching to `macos-15-intel`, or
attaching a self-hosted Intel runner.

Worth noting this removed a broken asset, not a working one: before
1.1.0 both macOS legs ran on `macos-latest`, which has been Apple
Silicon since early 2024, so every published `macOS-x86_64` binary was
really arm64.

### 10. Vignette re-building fails during `R CMD check` on Windows

`R-CMD-check` skips vignettes on the Windows job only (`build_args`
includes `--no-build-vignettes`, `args` includes `--ignore-vignettes`).
macOS and all three Linux jobs build them normally, and the pkgdown
workflow renders both articles on every push, so rendering is still
covered.

The Windows job failed about 3 minutes in — before any model fitting,
which takes ~20 minutes — reporting “Vignette re-building failed” for
both vignettes. The suspected cause is the `quarto` R package’s handling
of `system2(env = )` on Windows: it passes the environment assignment to
the CLI as an argument, so `quarto TMPDIR=C:\... -V` errors with
`Unknown command "TMPDIR=..."`. That is reproducible locally with quarto
R package 1.5.1 and CLI 1.7.30, though
[`quarto::quarto_version()`](https://quarto-dev.github.io/quarto-r/reference/quarto_version.html)
recovers from it and a direct
[`quarto::quarto_render()`](https://quarto-dev.github.io/quarto-r/reference/quarto_render.html)
succeeds — so the exact failure inside the vignette engine is not yet
pinned down.

**To diagnose properly** the job log is needed (reading it requires
admin rights on the repository). Worth re-testing after a `quarto` R
package update before investing further.

### 11. README binary URLs point at a release that must exist

`README.md` links to `v1.1.0` release assets. Those 404 until the
`v1.1.0` tag is pushed and `build.yml` publishes them.

### 12. Consider CRAN

The package now passes `R CMD check --as-cran` cleanly apart from the
environmental time NOTE. If CRAN submission is ever the goal, the
remaining work is a `cran-comments.md`, checking on r-hub/win-builder,
and confirming the 30-minute vignette build is acceptable (it very
likely is not — see item 7).
