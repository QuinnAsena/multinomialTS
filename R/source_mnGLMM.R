# Multinomial GLMM --------------------------------------------------------
#' Multinomial Generalized Linear Mixed Model (mnGLMM)
#'
#' Fits a multivariate generalized linear regression model for multinomially
#' distributed response data. This function estimates driver-species
#' relationships (B coefficients), observation-level dispersion,
#' and species covariance structures via maximum likelihood.
#'
#' @param Y A matrix of multinomially distributed count data
#'     (e.g., community count data).
#' @param X A matrix of covariates (predictors), which may be of mixed type
#'     Covariates should be scaled when appropriate. Can be \code{NULL}.
#' @param B.fixed A matrix indicating which B coefficients
#'     (driver-species relationships) to estimate. The number of columns must
#'     equal \code{ncol(Y)} (number of species), and the number of rows must
#'     equal \code{ncol(X) + 1} (intercept + covariates).
#' @param B.start A matrix of starting values for the B coefficients.
#'     Dimensions must match \code{B.fixed}.
#' @param sigma.fixed Fixed value for the overall model variance. Use \code{NA}
#'     to estimate it from the model.
#' @param sigma.start Starting value for estimating \code{sigma.fixed}
#'     (default is 0.1).
#' @param dispersion.fixed Fixed dispersion parameter to account for
#'     over- or under-dispersion. A value of 1 corresponds to no
#'     extra dispersion (pure multinomial).
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
#' @param compute.information.matrix Logical. If \code{TRUE},
#'     computes the observed information matrix.
#' @param hessian.method.args A list of control parameters passed to the
#'     numerical Hessian calculator (e.g., \code{numDeriv::hessian}).
#'
#' @return An object of class \code{"mnGLMM"} with components including:
#'   \item{B}{Estimated covariate coefficients (p x n matrix).}
#'   \item{sigma}{Estimated process noise scale.}
#'   \item{V}{Estimated species covariance matrix.}
#'   \item{dispersion}{Estimated observation dispersion parameter.}
#'   \item{mu}{Fitted multinomial probabilities (Tmax x n matrix).}
#'   \item{logLik}{Log-likelihood at the maximum.}
#'   \item{AIC}{Akaike Information Criterion.}
#'   \item{par}{Estimated parameters as a named vector.}
#'   \item{se}{Approximate standard errors from the observed information matrix.}
#'   \item{inv.information.matrix}{Inverse of the observed information matrix, or \code{NA} if not invertible.}
#'
#' @seealso [mnTS()], which uses the estimates from \code{mnGLMM} as starting
#'     values.
#'
#' @examples
#' data(story_pollen_matrix)
#' data(story_char_matrix)
#'
#' # mnGLMM cannot handle gaps, so use only the rows with pollen counts.
#' sample_idx <- which(rowSums(story_pollen_matrix) != 0)
#' Y <- story_pollen_matrix[sample_idx, ]
#' X <- scale(story_char_matrix)[sample_idx, , drop = FALSE]
#'
#' n <- ncol(Y)          # number of taxa
#' p <- ncol(X) + 1      # covariates + intercept
#'
#' # Fix the reference taxon (column 1) to 0; estimate the rest (NA = estimate).
#' B.fixed <- matrix(c(rep(0, p), rep(NA, (n - 1) * p)), p, n)
#' B.start <- matrix(c(rep(0, p), rep(0.01, (n - 1) * p)), p, n)
#'
#' glmm_mod <- mnGLMM(Y = Y, X = X, B.start = B.start, B.fixed = B.fixed,
#'                    V.fixed = diag(n))
#' summary(glmm_mod)
#'
#' @import minqa
#' @export


mnGLMM <- function(Y, X = NULL, B.fixed = NULL, B.start = NULL,
	sigma.fixed = NA, sigma.start = 0.1, dispersion.fixed = 1, dispersion.start = 1,
	V.fixed = diag(ncol(Y)), V.start = diag(ncol(Y)), method = "bobyqa", optim.control = NULL, maxit.optim = 1e+05, compute.information.matrix = TRUE,
	hessian.method.args = list(eps = 1e-04, d = 0.1)) {

	# check and name Y
	if (!is.matrix(Y)) {
		Y <- as.matrix(Y)
	}
	if (!is.null(X) && !is.null(B.fixed))
		if (!(ncol(X) == (dim(B.fixed)[1] - 1) & ncol(Y) == dim(B.fixed)[2])) {
			stop("dimensions of B.fixed should be c(ncol(X)+1, ncol(Y)).")
		}
	if (is.null(colnames(Y))) {
		colnames(Y) <- rep(NA, ncol(Y))
		for (i.col in 1:ncol(Y)) colnames(Y)[i.col] <- paste0("y", i.col)
	}
	n <- ncol(Y)

	# check and name X
	if (is.null(X)) {
		X <- matrix(1, nrow = nrow(Y), ncol = 1)
		colnames(X) <- "(intercept)"
	} else {
		if (!is.matrix(X)) {
			X <- as.matrix(X)
		}
		if (is.null(colnames(X))) {
			colnames(X) <- rep(NA, ncol(X))
			for (i.col in 1:ncol(X)) colnames(X)[i.col] <- paste0("x", i.col)
		}
		X <- cbind(matrix(1, nrow = nrow(X), ncol = 1), X)
		colnames(X)[1] <- "(intercept)"
	}
	if (nrow(X) != nrow(Y))
		stop("Variables X and Y have different lengths.")
	p <- ncol(X)

	# Defaults for B are built here rather than in the formals. X is rebuilt
	# just above (NULL becomes an intercept-only matrix; otherwise an intercept
	# column is prepended), and R evaluates default arguments lazily, so a
	# formal referring to X would see the rebuilt version and get both the
	# is.null(X) test and ncol(X) wrong.
	if (is.null(B.fixed))
		B.fixed <- matrix(c(rep(0, p), rep(NA, (n - 1) * p)), nrow = p, ncol = n)
	if (is.null(B.start))
		B.start <- matrix(0, nrow = p, ncol = n)

	# check and name V.fixed
	if (!is.matrix(V.fixed)) {
		V.fixed <- as.matrix(V.fixed)
	}
	if (nrow(V.fixed) != n | ncol(V.fixed) != n) {
		stop("V.fixed should have dimensions equal to the number of categories in Y.")
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

	# Check dispersion
	# dispersion is log-transformed below, so a value of 0 would put -Inf into
	# par.fixed, which the C++ side would treat as a free parameter.
	if (!is.na(dispersion.fixed) && (!is.finite(dispersion.fixed) || dispersion.fixed <= 0)) {
		stop("dispersion.fixed must be a positive, finite value, or NA to estimate it.")
	}
	if (!is.finite(dispersion.start) || dispersion.start <= 0) {
		stop("dispersion.start must be a positive, finite value.")
	}

	# log-transform dispersion
	dispersion.fixed <- log(dispersion.fixed)
	dispersion.start <- log(dispersion.start)

	par.fixed <- c(B.fixed, L.fixed, sigma.fixed, dispersion.fixed)

	B.names <- rep(NA, n * p)
	for (j in 1:n) for (i in 1:p) B.names[i + p * (j - 1)] <- paste0(colnames(X)[i], ".", colnames(Y)[j])
	L.names <- rep(NA, n * n)
	for (j in 1:n) for (i in 1:n) L.names[i + n * (j - 1)] <- paste0("v.", colnames(Y)[i], ".", colnames(Y)[j])

	names(par.fixed) <- c(B.names, L.names, "sigma", "dispersion")

	par.start <- c(B.start, L.start, sigma.start, dispersion.start)
	names(par.start) <- c(B.names, L.names, "sigma", "dispersion")
	par.start <- par.start[is.na(par.fixed)]

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
		opt <- optim(fn = mnGLMM.ml.wrapper, par = par.start, method = "Nelder-Mead", control = optim.control, par.fixed = par.fixed, Y = Y,
			X = X)
		if (opt$convergence != 0)
			cat("\nNelder-Mead optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$convergence
	}

	if (method == "BFGS") {
		opt <- optim(fn = mnGLMM.ml.wrapper, par = par.start, method = "BFGS", control = optim.control, par.fixed = par.fixed, Y = Y, X = X)
		if (opt$convergence != 0)
			cat("\nBFGS optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$convergence
	}

	if (method == "bobyqa") {
		opt <- bobyqa(fn = mnGLMM.ml.wrapper, par = par.start, control = optim.control, par.fixed = par.fixed, Y = Y, X = X)
		if (opt$ierr != 0)
			cat("\nbobyqa optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$ierr
		opt$value <- opt$fval
	}

	fitted <- mnGLMM_ml_cpp_listout(opt$par, par.fixed, Y, X)

	names(fitted$par) <- names(par.fixed)[is.na(par.fixed)]

	if (compute.information.matrix) {

		par.with.se <- fitted$par
		par.fixed.for.se <- par.fixed

		# Use the C++ likelihood, as mnTS does. The pure-R mnGLMM.ml is kept only
		# as the reference implementation the equivalency tests check against.
		H <- numDeriv::hessian(func = mnGLMM.ml.wrapper, x = par.with.se, par.fixed = par.fixed.for.se, Y = Y, X = X,
			method.args = hessian.method.args)
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
		# All of these are referenced when the results list is built below.
		H <- NULL
		se <- NULL
		par.with.se <- NULL
		inv.information.matrix <- NULL
		}

	fitted$par["dispersion"] <- exp(fitted$par["dispersion"])

	npar <- length(fitted$par)
	AIC <- -2 * fitted$logLik + 2 * npar

	results <- list(Y = Y, X = X, B = fitted$B, sigma = fitted$sigma, V = fitted$V, dispersion = fitted$dispersion, logLik = fitted$logLik,
		AIC = AIC, npar = npar, inv.information.matrix = inv.information.matrix, par.with.se = par.with.se, se = se,
		par = fitted$par, par.fixed = par.fixed, opt.convergence = opt.convergence)
	class(results) <- "mnGLMM"
	return(results)
}

# summary.mnGLMM --------------------------------------------------------
#' GLMM summary
#'
#' Prints the summary of a mnGLMM model.
#'
#' @param object An object of class "mnGLMM"
#' @param ... Additional arguments (currently unused).
#'
#' @return \code{object}, invisibly. Called for the side effect of printing
#'     the coefficient estimates from a fitted mnGLMM.
#'
#' @examples
#' data(story_pollen_matrix)
#' sample_idx <- which(rowSums(story_pollen_matrix) != 0)
#' n <- ncol(story_pollen_matrix)
#' glmm_mod <- mnGLMM(Y = story_pollen_matrix[sample_idx, ], V.fixed = diag(n))
#' summary(glmm_mod)
#'
#' @export

summary.mnGLMM <- function(object, ...) {
  mod <- object

	cat("\nCall: mnGLMM with Tmax = ", nrow(mod$Y), " n = ", ncol(mod$Y), "\n")

	cat("\nlogLik = ", mod$logLik, ",  AIC = ", mod$AIC, " [df = ", mod$npar, "]\n", sep = "")

	cat("\ndispersion parameter = ", mod$dispersion)

	if (!is.null(mod$inv.information.matrix)) {
		cat("\n\nFitted Coefficients with approximate se\n")
		if (!all(is.na(mod$inv.information.matrix))) {
			v.pars <- names(mod$par)[grep(names(mod$par), pattern = "^v\\.")]
			par.names.without.se <- c("dispersion","sigma", v.pars)
			par.with.se <- mod$par[!is.element(names(mod$par), par.names.without.se)]
			se <- mod$se[!is.element(names(mod$par), par.names.without.se)]

			t.scores <- par.with.se/se

			coef.table <- as.matrix(cbind(par.with.se, se, t.scores, 2 * pnorm(q = abs(t.scores), lower.tail = F)))
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

	cat("\n\nB = ", "\n")
	print(mod$B)

	cat("\nsigma = ", mod$sigma)

	cat("\n\nV = ", "\n")
	print(mod$V)

	cat("\n")

	invisible(object)
}


# print.mnGLMM --------------------------------------------------------
#' Print Method for Multinomial GLMMs
#'
#' Prints a fitted multinomial GLMM, using the same layout as
#' \code{\link{summary.mnGLMM}}.
#'
#' @param x An object of class \code{"mnGLMM"}, as returned by the
#'     \code{\link{mnGLMM}} function.
#' @param ... Additional arguments passed to \code{\link{summary.mnGLMM}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' data(story_pollen_matrix)
#' sample_idx <- which(rowSums(story_pollen_matrix) != 0)
#' n <- ncol(story_pollen_matrix)
#' glmm_mod <- mnGLMM(Y = story_pollen_matrix[sample_idx, ], V.fixed = diag(n))
#' print(glmm_mod)
#'
#' @export

print.mnGLMM <- function(x, ...) {
	summary(x, ...)
	invisible(x)
}
