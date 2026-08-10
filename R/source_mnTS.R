# mnTS ------------------------------------------------------------------------
#' Multinomial State-Space Model (mnTS)
#'
#' Fits a multinomial state-space model for multivariate count data,
#' allowing for latent temporal processes, covariate effects, and species
#' interactions.
#'
#' @param Y A matrix of multinomially distributed count data (e.g.,
#'     community count data).
#' @param X A matrix of covariates (predictors), which may be of mixed type.
#'     Covariates should be scaled when appropriate. Can be \code{NULL}.
#' @param Tsample A vector of row indices specifying the subset of
#'     observations in \code{Y} to treat as temporal samples.
#' @param B0.fixed A 1 by \code{ncol(Y)} matrix of species
#'     intercepts to estimate.
#' @param B0.start A matrix of starting values for \code{B0.fixed}.
#' @param B.fixed A matrix indicating which B coefficients (driver-species
#'     relationships) to estimate. Should have \code{ncol(Y)} columns and
#'     \code{ncol(X)} rows.
#' @param B.start A matrix of starting values for \code{B.fixed}. Dimensions
#'    of \code{B.start} should match \code{B.fixed}.
#' @param C.fixed A species-by-species matrix of interactions,
#'     indicating which interactions to estimate.
#' @param C.start A matrix of starting values for \code{C.fixed}.
#' @param sigma.fixed Fixed value for the overall model variance. Use
#'     \code{NA} to estimate it from the model.
#' @param sigma.start Starting value for estimating \code{sigma.fixed}.
#' @param dispersion.fixed Fixed dispersion parameter for observation-level
#'     variation. A value of 1 corresponds to no over- or under-dispersion.
#' @param dispersion.start Starting value for estimating
#'     \code{dispersion.fixed}.
#' @param V.fixed A species-by-species covariance matrix representing
#'     environmental variation.
#' @param V.start Starting values for \code{V.fixed}.
#' @param method Optimization method. Acceptable values include
#'     \code{"Nelder-Mead"}, \code{"BFGS"} (via \code{optim}),
#'     and \code{"bobyqa"} (via the \code{minqa} package).
#' @param optim.control Optional list of control parameters passed to the
#'     optimizer. See the \code{minqa} package documentation for details.
#' @param maxit.optim Maximum number of iterations for the optimizer
#'     (default is 1e+05). Increase if the optimizer needs more iterations.
#' @param compute.information.matrix Logical. If \code{TRUE}, computes the
#'     observed information matrix.
#' @param hessian.method.args A list of control parameters passed to the
#'     numerical Hessian calculator (e.g., \code{numDeriv::hessian}).
#'
#' @return An object of class \code{"mnTS"} with components including:
#'   \item{B0}{Estimated intercepts (1 x n matrix).}
#'   \item{C}{Estimated interaction matrix (n x n).}
#'   \item{B}{Estimated covariate coefficients (p x n matrix), or \code{NULL}.}
#'   \item{sigma}{Estimated process noise scale.}
#'   \item{V}{Estimated process noise covariance matrix.}
#'   \item{dispersion}{Estimated observation dispersion parameter.}
#'   \item{logLik}{Log-likelihood at the maximum.}
#'   \item{AIC}{Akaike Information Criterion.}
#'   \item{mu}{Fitted multinomial probabilities (Tmax x n matrix).}
#'   \item{y}{Fitted latent states on the log-ratio scale (Tmax x n matrix).}
#'   \item{se.y.fitted, se.mu.upper.fitted, se.mu.lower.fitted}{Kalman-filter standard errors and probability bands.}
#'   \item{par}{Estimated parameters as a named vector.}
#'   \item{se}{Approximate standard errors from the observed information matrix.}
#'   \item{inv.information.matrix}{Inverse of the observed information matrix, or \code{NA} if not invertible.}
#'
#' @seealso [mnGLMM()] for obtaining starting values, [bootstrap.mnTS()] for
#'     bootstrapped standard errors, and [simulate.mnTS()].
#'
#' @examples
#' # Tsample gives the time steps at which Y was observed: Y holds only the
#' # sampled rows, while X must cover every time step the model iterates over.
#' # This is a small, fast example (reference taxon plus two others over the
#' # first 60 time steps). See the package vignette for a full analysis.
#' data(story_pollen_matrix)
#' data(story_char_matrix)
#'
#' Y_all <- story_pollen_matrix[, 1:3]
#' sample_idx <- which(rowSums(Y_all) != 0)
#' sample_idx <- sample_idx[sample_idx <= 60]
#' Y <- Y_all[sample_idx, ]
#' X <- scale(story_char_matrix)[1:60, , drop = FALSE]
#' n <- ncol(Y)
#'
#' # NA = estimate, a number = hold fixed at that value.
#' B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
#' B.fixed <- matrix(NA, nrow = ncol(X), ncol = n)
#' B.fixed[, 1] <- 0                    # reference taxon
#' V.fixed <- matrix(NA, n, n)
#' V.fixed[1] <- 1                      # identifiability constraint
#' C.start <- 0.5 * diag(n)             # self-regulation only
#' C.fixed <- C.start
#' C.fixed[C.fixed != 0] <- NA
#'
#' \donttest{
#' ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
#'                B0.fixed = B0.fixed, B.fixed = B.fixed,
#'                C.start = C.start, C.fixed = C.fixed,
#'                V.fixed = V.fixed, dispersion.fixed = 1)
#' summary(ts_mod)
#' }
#'
#' @import numDeriv
#' @import minqa
#' @import stats
#' @export

mnTS <- function(Y, X = NULL, Tsample = 1:nrow(Y),
	B0.fixed = matrix(c(0, rep(NA, ncol(Y) - 1)), nrow = 1, ncol = ncol(Y)),
	B0.start = matrix(0, nrow = 1, ncol = ncol(Y)),
	C.fixed = diag(rep(NA, ncol(Y))), C.start = 0.01 * diag(ncol(Y)),
	B.fixed = if (is.null(X)) NULL else matrix(NA, nrow = ncol(X), ncol = ncol(Y)),
	B.start = if (is.null(X)) NULL else matrix(0, nrow = ncol(X), ncol = ncol(Y)),
	sigma.fixed = NA, sigma.start = 0.1, dispersion.fixed = 1, dispersion.start = 1, V.fixed = diag(ncol(Y)),
	V.start = diag(ncol(Y)),
	compute.information.matrix = TRUE, method = "bobyqa",
	optim.control = NULL, maxit.optim = 1e+05,
	hessian.method.args=list(eps=1e-4, d=0.0001, r=4, v=2)) {

	# check and name Y
	if (!is.matrix(Y)) {
		Y <- as.matrix(Y)
	}
	if (dim(Y)[2] > dim(Y)[1]) {
		stop("Time should run in the vertical direction.")
	}
	if (is.null(colnames(Y))) {
		colnames(Y) <- rep(NA, ncol(Y))
		for (i.col in 1:ncol(Y)) colnames(Y)[i.col] <- paste0("y", i.col)
	}
	if (any(diag(var(Y, na.rm = TRUE)) == 0)) {
		stop("The response Y (dependent variable) has no variation.")
	}
	n <- ncol(Y)

	# check Tsample
	# The filter steps forward from t = 1 to max(Tsample) and matches each row of
	# Y to a time step by position, so Tsample must be strictly increasing and
	# have exactly one entry per row of Y. Out-of-order indices used to be
	# accepted silently and produced a model fitted over the wrong time span.
	if (!is.numeric(Tsample) || anyNA(Tsample)) {
		stop("Tsample must be a numeric vector of row indices with no missing values.")
	}
	if (any(Tsample != round(Tsample)) || any(Tsample < 1)) {
		stop("Tsample must contain positive whole-number row indices.")
	}
	if (is.unsorted(Tsample, strictly = TRUE)) {
		stop("Tsample must be strictly increasing; the model steps forward through time.")
	}
	if (length(Tsample) != nrow(Y)) {
		stop("length(Tsample) (", length(Tsample), ") must equal nrow(Y) (", nrow(Y),
			"); Tsample gives the time step of each row of Y.")
	}
	Tmax <- max(Tsample)

	# check and name X
	if (!is.null(X)) {
		if (!is.matrix(X)) {
			X <- as.matrix(X)
		}
		if (is.null(colnames(X))) {
			colnames(X) <- rep(NA, ncol(X))
			for (i.col in 1:ncol(X)) colnames(X)[i.col] <- paste0("x", i.col)
		}
		if (nrow(X) < Tmax) {
			stop("X has ", nrow(X), " rows but the model iterates to time step ", Tmax,
				" (max(Tsample)); X must cover every time step.")
		}
		p <- ncol(X)
	} else {
		p <- 0
	}

	# check and name V.fixed
	if (!is.matrix(V.fixed)) {
		V.fixed <- as.matrix(V.fixed)
	}
	if (nrow(V.fixed) != n | ncol(V.fixed) != n) {
		stop("V.fixed should be diagonal with dimensions equal to the number of categories in Y.")
	}
	V.temp <- V.fixed
	V.temp[is.na(V.temp)] <- -1
	if (any(V.temp != t(V.temp))) {
		stop("V.fixed should be symmetric.")
	}
	if (is.null(colnames(V.fixed)))
		colnames(V.fixed) <- colnames(Y)
	if (is.null(rownames(V.fixed)))
		rownames(V.fixed) <- colnames(Y)

	# Check sigma
	# Sigma is squared in the calculations. This catch requested by a reviewer
	# should stop users trying to use negative inputs
	if (sigma.start < 0) {
	  stop("Sigma is less than 0")
	}

	# Check dispersion
	# dispersion is log-transformed below, so a value of 0 would put -Inf into
	# par.fixed. The C++ side treats any non-finite entry as a free parameter,
	# so that used to surface as an opaque "Mat::elem(): size mismatch".
	if (!is.na(dispersion.fixed) && (!is.finite(dispersion.fixed) || dispersion.fixed <= 0)) {
		stop("dispersion.fixed must be a positive, finite value, or NA to estimate it.")
	}
	if (!is.finite(dispersion.start) || dispersion.start <= 0) {
		stop("dispersion.start must be a positive, finite value.")
	}

	# handle V
	L.start <- chol(V.start)
	if (any(is.na(V.fixed[upper.tri(V.fixed, diag = FALSE)]))) {
		index.L.subset <- apply(V.fixed, 1, FUN = function(x) {
			any(is.na(x))
		})
		V.temp <- V.fixed
		V.temp[index.L.subset, ] <- 0
		V.temp[, index.L.subset] <- 0
		V.temp[index.L.subset, index.L.subset] <- V.start[index.L.subset, index.L.subset]
		L.fixed <- chol(V.temp)
		L.fixed[upper.tri(L.fixed, diag = TRUE) & is.na(V.fixed)] <- NA
	} else {
		if (any(is.na(diag(V.fixed)))) {
			L.fixed <- diag(n)
			diag(L.fixed)[is.na(diag(V.fixed))] <- NA
		} else {
			L.fixed <- chol(V.fixed)
		}
	}

	# log-transform dispersion
	dispersion.fixed <- log(dispersion.fixed)
	dispersion.start <- log(dispersion.start)

	# load parameters
	par.fixed <- c(B0.fixed, C.fixed, sigma.fixed, L.fixed, dispersion.fixed, B.fixed)

	B0.names <- colnames(Y)
	X.names <- colnames(X)
	C.names <- rep(NA, n * n)
	for (j in 1:n) for (i in 1:n) C.names[i + n * (j - 1)] <- paste0("sp.", colnames(Y)[i], ".", colnames(Y)[j])
	L.names <- rep(NA, n * n)
	for (j in 1:n) for (i in 1:n) L.names[i + n * (j - 1)] <- paste0("v.", colnames(Y)[i], ".", colnames(Y)[j])
	if (!is.null(X)) {
		B.names <- rep(NA, n * p)
		for (j in 1:n) for (i in 1:p) B.names[i + p * (j - 1)] <- paste0(colnames(X)[i], ".", colnames(Y)[j])
	} else {
		B.names <- NULL
	}

	names(par.fixed) <- c(B0.names, C.names, "sigma", L.names, "dispersion", B.names)

	par.start <- c(B0.start, C.start, sigma.start, L.start, dispersion.start, B.start)
	names(par.start) <- c(B0.names, C.names, "sigma", L.names, "dispersion", B.names)
	par.start <- par.start[is.na(par.fixed)]

	Y <- t(Y)
	if (!is.null(X))
		X <- t(X)

	if (!is.element(method, c("Nelder-Mead", "bobyqa", "BFGS")))
		stop("Acceptable methods are Nelder-Mead {optim},  BFGS {optim}, and bobyqa {minqa}.")

	if (is.null(optim.control)) {
		if (method == "bobyqa") {
			optim.control = list(maxfun = maxit.optim)
		} else {
			optim.control = list(maxit = maxit.optim)
		}
	}

	if (method == "Nelder-Mead") {
		opt <- optim(fn = mnTS.ml.wrapper, par = par.start, method = "Nelder-Mead", control = optim.control, par.fixed = par.fixed, Y = Y,
			X = X, Tsample = Tsample)
		if (opt$convergence != 0)
			cat("\nNelder-Mead optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$convergence
	}

	if (method == "BFGS") {
		opt <- optim(fn = mnTS.ml.wrapper, par = par.start, method = "BFGS", control = optim.control, par.fixed = par.fixed, Y = Y, X = X,
			Tsample = Tsample)
		if (opt$convergence != 0)
			cat("\nBFGS optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$convergence
	}

	if (method == "bobyqa") {
		opt <- bobyqa(fn = mnTS.ml.wrapper, par = par.start, control = optim.control,
		              par.fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
		if (opt$ierr != 0)
			cat("\nbobyqa optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$ierr
		opt$value <- opt$fval
	}

	fitted <- mnTS_ml_cpp_listout(opt$par, par_fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
	if(is.numeric(fitted)) stop("Convergence failed.")
	fitted_par <- as.vector(fitted$par)
	names(fitted_par) <- names(par.fixed)[is.na(par.fixed)]

	# Note: do not delete fitted$B to signal "no covariates". Removing the element
	# makes fitted$B fall through to partial matching and silently return
	# fitted$B0 instead. Build the output separately.
	if (!is.null(X)) {
		rownames(fitted$B) <- B0.names
		colnames(fitted$B) <- X.names
		B.out <- t(fitted$B)
	} else {
		B.out <- NULL
	}

	rownames(fitted$B0) <- B0.names
	colnames(fitted$B0) <- "(intercept)"

	rownames(fitted$C) = colnames(fitted$C) = B0.names
	rownames(fitted$V) = colnames(fitted$V) = B0.names

	# names(fitted_par) < names(par.fixed)[is.na(par.fixed)]

	if (compute.information.matrix) {
		# compute information matrix excluding sigma and dispersion
		par.with.se <- fitted_par
		par.fixed.for.se <- par.fixed


	if (!is.na(fitted_par["dispersion"])) {
  pick.pars.with.se <- !is.element(names(fitted_par), "dispersion")
  par.with.se <- fitted_par[pick.pars.with.se]
  par.fixed.for.se <- par.fixed
  par.fixed.for.se[names(fitted_par[!pick.pars.with.se])] <- fitted_par[!pick.pars.with.se]
  } else {
  par.with.se <- fitted_par
  par.fixed.for.se <- par.fixed
  }

		H <- hessian(func = mnTS.ml.wrapper, x = par.with.se, par.fixed = par.fixed.for.se, Y = Y, X = X, Tsample = Tsample, method.args = hessian.method.args)

		if (any(is.nan(H)) || rcond(H) < 1e-12) {
			inv.information.matrix <- matrix(NA, nrow(H), ncol(H))
		} else {
			inv.information.matrix <- solve(H)
		}
		se <- diag(inv.information.matrix)^0.5
		names(se) <- names(par.with.se)
		rownames(inv.information.matrix) <- names(par.with.se)
		colnames(inv.information.matrix) <- names(par.with.se)

		} else {
		  # Every one of these is referenced when the results list is built below,
		  # so they all need a value here, not just the information matrix.
		  H <- NULL
		  se <- NULL
		  par.with.se <- NULL
		  inv.information.matrix <- NULL
		}

	# Undo the log transform, but only when dispersion was actually estimated.
	# fitted_par holds free parameters only, so when dispersion is fixed there is
	# no "dispersion" element and this assignment used to *append* one with the
	# value exp(NA) = NA. That inflated npar by 1 and so AIC by 2.
	if ("dispersion" %in% names(fitted_par)) {
		fitted_par["dispersion"] <- exp(fitted_par["dispersion"])
	}

	if (!is.null(X))
		tX <- t(X)
	else tX <- NULL

	npar <- length(fitted_par)
	AIC <- -2 * fitted$logLik + 2 * npar

	results <- list(Y = t(Y), X = tX, Tsample = Tsample, B0 = t(fitted$B0), C = fitted$C,
		sigma = fitted$sigma, dispersion = fitted$dispersion,
		V = fitted$V, B = B.out, logLik = fitted$logLik,
		AIC = AIC, npar = npar, inv.information.matrix = inv.information.matrix,
		par.with.se = par.with.se, se = se, information.matrix = H,
		y = fitted$y, mu = fitted$mu, par = fitted_par, par.fixed = par.fixed,
		opt.convergence = opt.convergence, se.y.fitted = fitted$se.y.fitted,
		se.mu.upper.fitted = fitted$se.mu.upper.fitted, se.mu.lower.fitted = fitted$se.mu.lower.fitted
		)
	class(results) <- "mnTS"
	return(results)
}

# summary.mnTS --------------------------------------------------------
#' Summary Method for Multinomial State-Space Models
#'
#' Provides a summary of the fitted multinomial state-space model.
#'
#' @param object An object of class \code{"mnTS"}, typically returned by the
#'     \code{\link{mnTS}} function.
#' @param ... Additional arguments (currently unused).
#'
#' @return \code{object}, invisibly. Called for the side effect of printing
#'     the model summary to the console.
#'
#' @examples
#' # A small, fast example: the reference taxon plus two others over the first
#' # 60 time steps. See the package vignette for a full analysis.
#' data(story_pollen_matrix)
#' data(story_char_matrix)
#'
#' Y_all <- story_pollen_matrix[, 1:3]
#' sample_idx <- which(rowSums(Y_all) != 0)
#' sample_idx <- sample_idx[sample_idx <= 60]
#' Y <- Y_all[sample_idx, ]
#' X <- scale(story_char_matrix)[1:60, , drop = FALSE]
#' n <- ncol(Y)
#'
#' # NA = estimate, a number = hold fixed at that value.
#' B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
#' B.fixed <- matrix(NA, nrow = ncol(X), ncol = n)
#' B.fixed[, 1] <- 0                    # reference taxon
#' V.fixed <- matrix(NA, n, n)
#' V.fixed[1] <- 1                      # identifiability constraint
#' C.start <- 0.5 * diag(n)             # self-regulation only
#' C.fixed <- C.start
#' C.fixed[C.fixed != 0] <- NA
#'
#' \donttest{
#' ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
#'                B0.fixed = B0.fixed, B.fixed = B.fixed,
#'                C.start = C.start, C.fixed = C.fixed,
#'                V.fixed = V.fixed, dispersion.fixed = 1)
#' summary(ts_mod)
#' }
#'
#' @export

summary.mnTS <- function(object, ...) {
  mod <- object

	cat("\nCall: mnTS with Tmax = ", nrow(mod$Y), " n = ", ncol(mod$Y), "\n")

	cat("\nlogLik = ", mod$logLik, ",  AIC = ", mod$AIC, " [df = ", mod$npar, "]\n", sep = "")

	cat("\ndispersion parameter = ", mod$dispersion)

	if (!is.null(mod$inv.information.matrix)) {
	  if (!all(is.na(mod$inv.information.matrix))) {
	    cat("\n\nFitted Coefficients with approximate se\n")
	    v.pars <- names(mod$par)[grep(names(mod$par), pattern = "^v\\.")]
	    par.names.noshow.se <- c("dispersion","sigma", v.pars)
	    par.show.se <- mod$par[!is.element(names(mod$par), par.names.noshow.se)]
	    se.show.se <- mod$se[!is.element(names(mod$se), par.names.noshow.se)]
	    t.scores <- par.show.se/se.show.se
	    coef.table <- as.matrix(cbind(par.show.se, se.show.se, t.scores, 2 * pnorm(q = abs(t.scores), lower.tail = F)))
	    colnames(coef.table) <- c("Coef.", "se", "t", "P")
	    print(coef.table)
	  } else {
	    coef.table <- as.matrix(mod$par.with.se)
	    colnames(coef.table) <- "Coef."
	    print(coef.table)
	    cat("\ninformation matrix not invertible so se's not calculated")
	  }
	}

	cat("\n\nOverall model")

	cat("\n\nB0 = ", "\n")
	print(mod$B0)

	if (!is.null(mod$X)) {
		cat("\n\nB = ", "\n")
		print(mod$B)
	}

	cat("\n\nC = ", "\n")
	print(mod$C)

	cat("\nsigma = ", mod$sigma)

	cat("\n\nV = ", "\n")
	print(mod$V)

	invisible(object)
}


# print.mnTS ----------------------------------------------------------
#' Print Method for Multinomial State-Space Models
#'
#' Prints a fitted multinomial state-space model, using the same layout as
#' \code{\link{summary.mnTS}}.
#'
#' @param x An object of class \code{"mnTS"}, as returned by the
#'     \code{\link{mnTS}} function.
#' @param ... Additional arguments passed to \code{\link{summary.mnTS}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' # A small, fast example: the reference taxon plus two others over the first
#' # 60 time steps. See the package vignette for a full analysis.
#' data(story_pollen_matrix)
#' data(story_char_matrix)
#'
#' Y_all <- story_pollen_matrix[, 1:3]
#' sample_idx <- which(rowSums(Y_all) != 0)
#' sample_idx <- sample_idx[sample_idx <= 60]
#' Y <- Y_all[sample_idx, ]
#' X <- scale(story_char_matrix)[1:60, , drop = FALSE]
#' n <- ncol(Y)
#'
#' # NA = estimate, a number = hold fixed at that value.
#' B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
#' B.fixed <- matrix(NA, nrow = ncol(X), ncol = n)
#' B.fixed[, 1] <- 0                    # reference taxon
#' V.fixed <- matrix(NA, n, n)
#' V.fixed[1] <- 1                      # identifiability constraint
#' C.start <- 0.5 * diag(n)             # self-regulation only
#' C.fixed <- C.start
#' C.fixed[C.fixed != 0] <- NA
#'
#' \donttest{
#' ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
#'                B0.fixed = B0.fixed, B.fixed = B.fixed,
#'                C.start = C.start, C.fixed = C.fixed,
#'                V.fixed = V.fixed, dispersion.fixed = 1)
#' print(ts_mod)
#' }
#'
#' @export

print.mnTS <- function(x, ...) {
	summary(x, ...)
	invisible(x)
}



# coef.mnTS -----------------------------------------------------------
#' Coefficients for Multinomial State-Space Model
#'
#' Extracts estimated coefficients from a fitted multinomial state-space model.

#'
#' @param object An object of class \code{"mnTS"}, as returned by the
#'     \code{\link{mnTS}} function.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix of estimated coefficients. When an information matrix
#'     is available it has columns \code{Coef.}, \code{se}, \code{t} and
#'     \code{P}; otherwise it holds point estimates only.
#'
#' @examples
#' # A small, fast example: the reference taxon plus two others over the first
#' # 60 time steps. See the package vignette for a full analysis.
#' data(story_pollen_matrix)
#' data(story_char_matrix)
#'
#' Y_all <- story_pollen_matrix[, 1:3]
#' sample_idx <- which(rowSums(Y_all) != 0)
#' sample_idx <- sample_idx[sample_idx <= 60]
#' Y <- Y_all[sample_idx, ]
#' X <- scale(story_char_matrix)[1:60, , drop = FALSE]
#' n <- ncol(Y)
#'
#' # NA = estimate, a number = hold fixed at that value.
#' B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
#' B.fixed <- matrix(NA, nrow = ncol(X), ncol = n)
#' B.fixed[, 1] <- 0                    # reference taxon
#' V.fixed <- matrix(NA, n, n)
#' V.fixed[1] <- 1                      # identifiability constraint
#' C.start <- 0.5 * diag(n)             # self-regulation only
#' C.fixed <- C.start
#' C.fixed[C.fixed != 0] <- NA
#'
#' \donttest{
#' ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
#'                B0.fixed = B0.fixed, B.fixed = B.fixed,
#'                C.start = C.start, C.fixed = C.fixed,
#'                V.fixed = V.fixed, dispersion.fixed = 1)
#' coef(ts_mod)
#' }
#'
#' @export

coef.mnTS <- function(object, ...) {
  mod <- object

  if (!is.null(mod$inv.information.matrix)) {
    if (!all(is.na(mod$inv.information.matrix))) {
       v.pars <- names(mod$par)[grep(names(mod$par), pattern = "^v\\.")]
      par.names.noshow.se <- c("dispersion","sigma", v.pars)
      par.show.se <- mod$par[!is.element(names(mod$par), par.names.noshow.se)]
      se.show.se <- mod$se[!is.element(names(mod$se), par.names.noshow.se)]
      t.scores <- par.show.se/se.show.se
      coef.table <- as.matrix(cbind(par.show.se, se.show.se, t.scores, 2 * pnorm(q = abs(t.scores), lower.tail = F)))
      colnames(coef.table) <- c("Coef.", "se", "t", "P")
    } else {
      coef.table <- as.matrix(mod$par.with.se)
      colnames(coef.table) <- "Coef."
    }
  } else {
    # Fitted with compute.information.matrix = FALSE: point estimates only.
    coef.table <- as.matrix(mod$par)
    colnames(coef.table) <- "Coef."
  }
  return(coef.table)
}


# simulate.mnTS -------------------------------------------------------
#' Simulate Data from a Multinomial State-Space Model
#'
#' Generates simulated datasets from a fitted mnTS object.
#'
#' @param object An object of class \code{"mnTS"}, as returned by the
#'     \code{\link{mnTS}} function.
#' @param nsim Number of simulated datasets to generate (default 1).
#' @param seed Optional seed for the random number generator. Passed to
#'     \code{\link{set.seed}}; the previous RNG state is restored on exit.
#' @param size Multinomial sample size (total counts) to draw at each sampled
#'     time step. Defaults to \code{NULL}, which reuses the observed totals
#'     from the data the model was fitted to, so simulated replicates carry the
#'     same sampling effort as the data. May also be a single number, or a
#'     vector with one entry per element of \code{object$Tsample}.
#' @param ... Additional arguments (currently unused).
#'
#' @return When \code{nsim = 1}, a list with the simulated response \code{Y}
#'     and the model components used to generate it (\code{X}, \code{B0},
#'     \code{C}, \code{sigma}, \code{V}, \code{B}). When \code{nsim > 1}, an
#'     unnamed list of \code{nsim} such lists. Rows of \code{Y} at time steps
#'     not in \code{object$Tsample} are \code{NA}, matching the observation
#'     design of the input data.
#'
#' @examples
#' # A small, fast example: the reference taxon plus two others over the first
#' # 60 time steps. See the package vignette for a full analysis.
#' data(story_pollen_matrix)
#' data(story_char_matrix)
#'
#' Y_all <- story_pollen_matrix[, 1:3]
#' sample_idx <- which(rowSums(Y_all) != 0)
#' sample_idx <- sample_idx[sample_idx <= 60]
#' Y <- Y_all[sample_idx, ]
#' X <- scale(story_char_matrix)[1:60, , drop = FALSE]
#' n <- ncol(Y)
#'
#' # NA = estimate, a number = hold fixed at that value.
#' B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
#' B.fixed <- matrix(NA, nrow = ncol(X), ncol = n)
#' B.fixed[, 1] <- 0                    # reference taxon
#' V.fixed <- matrix(NA, n, n)
#' V.fixed[1] <- 1                      # identifiability constraint
#' C.start <- 0.5 * diag(n)             # self-regulation only
#' C.fixed <- C.start
#' C.fixed[C.fixed != 0] <- NA
#'
#' \donttest{
#' ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
#'                B0.fixed = B0.fixed, B.fixed = B.fixed,
#'                C.start = C.start, C.fixed = C.fixed,
#'                V.fixed = V.fixed, dispersion.fixed = 1)
#' sim <- simulate(ts_mod, seed = 1)
#' # Simulated totals match the observed sampling effort.
#' rowSums(sim$Y[sample_idx, ])[1:5]
#' rowSums(ts_mod$Y)[1:5]
#' }
#'
#' @importFrom mvtnorm rmvnorm
#' @export


simulate.mnTS <- function(object, nsim = 1, seed = NULL, size = NULL, ...) {
  mod <- object

  # Standard stats::simulate seed handling: set the seed for this call only and
  # restore whatever the caller had afterwards.
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1) {
    stop("nsim must be a single positive whole number.")
  }
  nsim <- as.integer(nsim)

  Y <- t(mod$Y)
  n <- nrow(Y)
  Tsample <- mod$Tsample
  # Iterate over the same time span as the fitted model. nrow(mod$X) was used
  # here previously, which is NULL for a model fitted without covariates.
  Tmax <- max(Tsample)

  # Multinomial sample sizes. These used to be hard-coded to 100 regardless of
  # the data, which meant bootstrap.mnTS() resampled with far less sampling effort
  # than the real counts and so overstated the bootstrap standard errors.
  # colSums(Y) is the observed total at each sampled time, i.e. the
  # sum(Y[, which(t == Tsample)]) used by the likelihood.
  size.t <- rep(NA_real_, Tmax)
  if (is.null(size)) {
    size.t[Tsample] <- colSums(Y)
  } else if (length(size) == 1L) {
    size.t[Tsample] <- size
  } else if (length(size) == length(Tsample)) {
    size.t[Tsample] <- size
  } else {
    stop("size must be a single value or have one entry per element of object$Tsample (",
         length(Tsample), ").")
  }
  if (any(!is.finite(size.t[Tsample]) | size.t[Tsample] < 1)) {
    stop("size must be finite and at least 1 at every sampled time step.")
  }

  if (!is.null(mod$X)) {
    X <- t(mod$X)
    B <- t(mod$B)
  } else {
    X <- NULL
    B <- matrix(0, n, 1)
  }

  B0 <- t(mod$B0)
  C <- mod$C
  sigma <- mod$sigma

  Se <- sigma^2 * mod$V

  sim.one <- function() {

    Y.sim <- matrix(NA, n, Tmax)

    for (t in 1:Tmax) {

      if (t == 1) {

        # Initial values at the stationary distribution

        if (!is.null(X)) {
          y <- B0 + B %*% X[, t, drop = FALSE]
        } else {
          y <- B0
        }
        y <- y + t(rmvnorm(n = 1, sigma = Se))
      } else {

        # PREDICTION EQUATIONS

        if (!is.null(X)) {
          y <- B0 + C %*% (y - B0 - B %*% X[, t - 1, drop = FALSE]) + B %*% X[, t, drop = FALSE]
        } else {
          y <- B0 + C %*% (y - B0)
        }
        y <- y + t(rmvnorm(n = 1, sigma = Se))
      }

      mu <- as.numeric(exp(y))
      mu <- mu/sum(mu)

      # Counts only exist where the data were sampled; other time steps stay NA.
      if (!is.na(size.t[t])) {
        Y.sim[, t] <- rmultinom(n = 1, size = size.t[t], prob = mu)
      }

    }
    list(Y = t(Y.sim), X = mod$X, B0 = mod$B0, C = C, sigma = abs(sigma),
         V = mod$V, B = mod$B)
  }

  # nsim = 1 returns a single simulation directly, so existing code that treats
  # the result as one dataset keeps working.
  if (nsim == 1L) {
    return(sim.one())
  }
  lapply(seq_len(nsim), function(i) sim.one())
}


# bootstrap ------------------------------------------------------------
#' Multinomial state-space bootstrap
#'
#' This is a generic bootstrap function for models.
#'
#' @param object An object for which the bootstrap method is defined.
#' @param ... Additional arguments passed to the method.
#'
#' @return The value returned by the method dispatched on \code{object};
#'     see [bootstrap.mnTS()].
#'
#' @seealso [bootstrap.mnTS()]
#'
#' @export
bootstrap <- function(object, ...) {
  UseMethod("bootstrap")
}


# boot (deprecated) ---------------------------------------------------
#' Bootstrap a model (deprecated)
#'
#' Deprecated in favour of [bootstrap()], which does the same thing under a name
#' that does not collide with \code{boot::boot()} from the recommended
#' \pkg{boot} package. This wrapper warns and forwards, and will be removed in a
#' future release.
#'
#' @param object An object of class \code{mnTS}.
#' @param ... Passed to [bootstrap()].
#'
#' @return The value returned by [bootstrap()].
#'
#' @seealso [bootstrap()], which replaces this function.
#'
#' @export
boot <- function(object, ...) {
  .Deprecated("bootstrap", package = "multinomialTS",
    msg = paste("multinomialTS::boot() is deprecated and will be removed in a",
                "future release; use bootstrap() instead. The old name clashed",
                "with boot::boot() from the recommended boot package, so",
                "whichever package was attached second masked the other."))
  bootstrap(object, ...)
}


# bootstrap.mnTS -------------------------------------------------------
#' Bootstrap method for mnTS objects
#'
#' Provides bootstrapped confidence intervals for an object of class \code{mnTS}.
#'
#' @param object An object of class \code{mnTS}.
#' @param reps Integer. Number of bootstrap replications. More is better; the
#'     package vignette discusses how to choose a value.
#' @param dispersion.fixed Numeric. Dispersion parameter (default is 1).
#' @param maxit.optim Integer. Max iterations for optimization (default 100000).
#' @param ... Additional arguments (not currently used).
#'
#' @return A list with three elements: \code{coef.table.boot}, a table of
#'     bootstrapped coefficients with standard deviations, t-scores, P-values
#'     and 68% bounds; \code{all_mods_pars}, one row of coefficients per
#'     successful replicate; and \code{all_mods}, the fitted replicate models.
#'
#' @examples
#' # A small, fast example: the reference taxon plus two others over the first
#' # 60 time steps. See the package vignette for a full analysis.
#' data(story_pollen_matrix)
#' data(story_char_matrix)
#'
#' Y_all <- story_pollen_matrix[, 1:3]
#' sample_idx <- which(rowSums(Y_all) != 0)
#' sample_idx <- sample_idx[sample_idx <= 60]
#' Y <- Y_all[sample_idx, ]
#' X <- scale(story_char_matrix)[1:60, , drop = FALSE]
#' n <- ncol(Y)
#'
#' # NA = estimate, a number = hold fixed at that value.
#' B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
#' B.fixed <- matrix(NA, nrow = ncol(X), ncol = n)
#' B.fixed[, 1] <- 0                    # reference taxon
#' V.fixed <- matrix(NA, n, n)
#' V.fixed[1] <- 1                      # identifiability constraint
#' C.start <- 0.5 * diag(n)             # self-regulation only
#' C.fixed <- C.start
#' C.fixed[C.fixed != 0] <- NA
#'
#' \donttest{
#' ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
#'                B0.fixed = B0.fixed, B.fixed = B.fixed,
#'                C.start = C.start, C.fixed = C.fixed,
#'                V.fixed = V.fixed, dispersion.fixed = 1)
#' boot_out <- bootstrap(ts_mod, reps = 3)
#' boot_out$coef.table.boot
#' }
#'
#' @export
#' @method bootstrap mnTS

bootstrap.mnTS <- function(object, reps, dispersion.fixed = 1, maxit.optim = 1e+05, ...) {
  mod <- object

  if (length(reps) != 1L || !is.finite(reps) || reps < 1) {
    stop("reps must be a single positive whole number.")
  }
  reps <- as.integer(reps)

  # simulate() returns one dataset for nsim = 1 and a list of datasets above
  # that, so normalise to a list either way.
  sims <- if (reps == 1L) list(simulate(mod)) else simulate(mod, nsim = reps)
  fit_mods <- vector(mode = "list", length = reps)

  for (j in seq_along(sims)) {

    X <- sims[[j]]$X
    Y <- sims[[j]]$Y
    n <- ncol(Y)
    Y <- Y[mod$Tsample, ]

    # With no covariates there are no B coefficients to estimate; passing a
    # non-NULL B.fixed here would add unnamed slots to the parameter vector.
    if (!is.null(X)) {
      p <- ncol(X)
      B.fixed <- matrix(NA, p, n)
      B.fixed[, 1] <- 0
      B.start <- mod$B
    } else {
      B.fixed <- NULL
      B.start <- NULL
    }

    V.fixed <- matrix(NA, n, n)
    V.fixed[1] <- 1

    C.fixed <- mod$C
    C.fixed[C.fixed != 0] <- NA

    fit_mods[[j]] <- tryCatch(
      {
        ss_mod <- mnTS(Y = Y, X = X,
                       B0.start = mod$B0,
                       B.start = B.start, B.fixed = B.fixed,
                       C.fixed = C.fixed, C.start = mod$C*0.5,
                       V.fixed = V.fixed, V.start = mod$V,
                       sigma.start = mod$sigma,
                       dispersion.fixed = dispersion.fixed,
                       Tsample = mod$Tsample, maxit.optim = maxit.optim)
      },
      error=function(cond) {
        message("Here's the original error message:")
        message(conditionMessage(cond))
        return(list(err = cond, full_mod = sims[[j]]))
      }
    )
  }
  # Filter out any models that had an error. inherits() rather than
  # class() == "mnTS", which compares a list to a string and breaks for any
  # object whose class attribute is not length 1.
  fit_mods <- fit_mods[vapply(fit_mods, inherits, logical(1), what = "mnTS")]
  if (!length(fit_mods)) {
    stop("No bootstrap replicate fitted successfully; see the error messages above.")
  }

  v.pars <- names(mod$par)[grep(names(mod$par), pattern = "^v\\.")]
  par.names.noshow.se <- c("dispersion","sigma", v.pars)

  pars <- lapply(seq_along(fit_mods), \(rep) {
    par.show.se <- fit_mods[[rep]]$par[!is.element(names(fit_mods[[rep]]$par), par.names.noshow.se)]
    par.all <- c(par.show.se, logLik = fit_mods[[rep]]$logLik, opt.convergence = fit_mods[[rep]]$opt.convergence)
  })
  all_pars <- do.call(rbind, pars)
  # drop = FALSE: with a single surviving replicate this subsetting would
  # otherwise collapse to a vector and matrixStats would reject it.
  most_pars <- all_pars[, !colnames(all_pars) %in% c("logLik", "opt.convergence"),
                        drop = FALSE]

  most_pars_mean <- matrixStats::colMeans2(most_pars, na.rm = TRUE)
  most_pars_sd <- matrixStats::colSds(most_pars, na.rm = TRUE)
  most_pars_upper_68 <- most_pars_mean + most_pars_sd
  most_pars_lower_68 <- most_pars_mean - most_pars_sd
  t_scores <- most_pars_mean / most_pars_sd
  p_vals <- 2 * pnorm(q = abs(t_scores), lower.tail = FALSE)

  # Count replicates that failed to converge. Summing the codes themselves
  # overcounts: bobyqa reports 0 on success and 1-5 for different failures, so
  # a single ierr = 5 used to be reported as five failures.
  n_fail_converged <- rep(sum(all_pars[, colnames(all_pars) %in% "opt.convergence"] != 0), length(most_pars_mean))
  n_total <- rep(nrow(all_pars), length(most_pars_mean))

  boot_bind <- cbind(
    boot_mean = most_pars_mean,
    se = most_pars_sd,
    t = t_scores,
    P = p_vals,
    upper68 = most_pars_upper_68,
    lower68 = most_pars_lower_68,
    n_fail_converged = n_fail_converged,
    n_total = n_total
  )
  return(list(coef.table.boot = boot_bind, all_mods_pars = all_pars, all_mods = fit_mods))
}
