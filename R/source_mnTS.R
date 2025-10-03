# Multinomial state-space ML -----------------------------------------------------
#' Multinomial state-space Maximum likelihood
#'
#' This is an internal function called by mnTS
#'

mnTS.ml <- function(par, par.fixed, Y, X = NULL, Tsample, fitted.values = FALSE) {

	parameters <- par.fixed
	parameters[is.na(par.fixed)] <- par
	n <- nrow(Y)
	Tmax <- Tsample[length(Tsample)]

	B0 <- matrix(parameters[1:n], n, 1, byrow = F)
	C <- matrix(parameters[(n + 1):(n + n * n)], n, n, byrow = F)
	sigma <- parameters[n + n * n + 1]

	L <- matrix(parameters[(n + n * n + 1 + 1):(n + n * n + 1 + n * n)], n, n)
	S <- t(L) %*% L
	#Ssd <- diag(diag(S)^-.5)
	#Se <- sigma^2 * Ssd %*% S %*% Ssd
	Ssd <- S[1, 1]
	Se <- sigma^2 * S/Ssd

	dispersion <- exp(parameters[n + n * n + 1 + n * n + 1])

	if (!is.null(X)) {
		p <- nrow(X)
		B <- matrix(parameters[(n + n * n + 1 + n * n + 1 + 1):(n + n * n + 1 + 1 + n * n + p * n)], n, p, byrow = TRUE)
	}

	logFt <- 0
	vFv <- 0

	if (fitted.values) {
		y.fitted <- matrix(NA, n, Tmax)
		mu.fitted <- matrix(NA, n, Tmax)
		se.y.fitted <- matrix(NA, n, Tmax)
		se.mu.upper.fitted <- matrix(NA, n, Tmax)
		se.mu.lower.fitted <- matrix(NA, n, Tmax)

	}

	for (t in 1:Tmax) {

		if (t == 1) {

			# Initial values at the stationary distribution

			if (!is.null(X)) {
				y <- B0 + B %*% X[, t, drop = FALSE]
			} else {
				y <- B0
			}
			if (max(abs(eigen(C)$values)) < 1) {
				PP <- solve(diag(n * n) - kronecker(C, C)) %*% matrix(Se, nrow = n * n)
				PP <- matrix(PP, n, n)
			} else {
				PP <- Se
			}

		} else {

			# PREDICTION EQUATIONS

			if (!is.null(X)) {
				y <- B0 + C %*% (y - B0 - B %*% X[, t - 1, drop = FALSE]) + B %*% X[, t, drop = FALSE]
			} else {
				y <- B0 + C %*% (y - B0)
			}
			# y <- B0 + C %*% (y - B0) + B %*% X[, t, drop = FALSE]
			PP <- C %*% PP %*% t(C) + Se
		}

		if (!any(is.na(Y[, which(t == Tsample)])) & is.element(t, Tsample)) {

			# UPDATING EQUATIONS

			size <- sum(Y[, which(t == Tsample)])

			mu <- as.numeric(exp(y))
			mu <- mu/sum(mu)

			W <- diag((mu * (1 - mu))^-1) %*% (diag(mu) - matrix(mu, n, 1) %*% matrix(mu, 1, n)) %*% diag((mu * (1 - mu))^-1)
			M <- diag(1 + mu) - (matrix(1, n, 1) %*% matrix(mu, 1, n))
			D <- diag(mu * (1 - mu))
			ZDM <- (D %*% M)[-1, ]
			FF <- ZDM %*% PP %*% t(ZDM) + dispersion * (D %*% (W/size) %*% D)[-1, -1]

			if (any(is.nan(FF)) || rcond(FF) < 1e-12)
				return(1e+10)
			invF <- solve(FF)

			v <- Y[-1, which(t == Tsample)]/size - mu[-1]

			y <- y + PP %*% t(ZDM) %*% invF %*% v
			PP <- PP - PP %*% t(ZDM) %*% invF %*% ZDM %*% PP

			logFt <- logFt + determinant(FF)$modulus[1]
			vFv <- vFv + t(v) %*% invF %*% v

		}
		if (fitted.values) {
			y.fitted[, t] <- y
			mu <- as.numeric(exp(y))
			mu <- mu/sum(mu)
			mu.fitted[, t] <- mu
			se.y.fitted[, t] <- sqrt(diag(PP))
			se.mu.upper.fitted[, t] <- exp(y.fitted[, t] + se.y.fitted[, t])/sum(as.numeric(exp(y)))
			se.mu.lower.fitted[, t] <- exp(y.fitted[, t] - se.y.fitted[, t])/sum(as.numeric(exp(y)))
		}
	}
	LL <- logFt + vFv


	if (is.complex(LL))
		LL <- 10^10

	if (!fitted.values) {
		return(LL/2)
	} else {
		B0 <- t(B0)
		colnames(B0) <- rownames(Y)
		rownames(B0) <- "mean"

		colnames(C) <- rownames(Y)
		rownames(C) <- rownames(Y)

		V <- (t(L) %*% L)
		colnames(V) <- rownames(Y)
		rownames(V) <- rownames(Y)

		if (!is.null(X)) {
			B <- t(B)
			colnames(B) <- rownames(Y)
			rownames(B) <- rownames(X)
		} else {
			B <- NULL
		}

		par <- parameters[is.na(par.fixed)]
		logLik <- -((n - 1) * length(Tsample)/2) * log(2 * pi) - LL/2

		sigma = abs(sigma)
		par[names(par) == "sigma"] <- abs(par[names(par) == "sigma"])

		return(list(B0 = B0, C = C, sigma = sigma, V = V, B = B, dispersion = dispersion, logLik = logLik, par = par, mu = t(mu.fitted),
			y = t(y.fitted), se.y.fitted = t(se.y.fitted), se.mu.upper.fitted = t(se.mu.upper.fitted), se.mu.lower.fitted = t(se.mu.lower.fitted)))
	}
}

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
#' @return An object of class \code{"mnTS"} containing estimated parameters.
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
	nn <- n - 1
	Tmax <- Tsample[length(Tsample)]

	# check and name X
	if (!is.null(X)) {
		if (!is.matrix(X)) {
			X <- as.matrix(X)
		}
		if (is.null(colnames(X))) {
			colnames(X) <- rep(NA, ncol(X))
			for (i.col in 1:ncol(X)) colnames(X)[i.col] <- paste0("x", i.col)
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
	# if (method == "bobyqa")
	# 	require(minqa)

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

	rownames(fitted$B) <- B0.names
	colnames(fitted$B) <- X.names

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

		# compute information matrix excluding V, sigma and dispersion
		par.with.se.cond <- fitted_par
		par.fixed.for.se.cond <- par.fixed
		if (!is.na(fitted_par["dispersion"])) {
			#pick.pars.with.se.cond <- !is.element(names(fitted_par), c("sigma", L.names, "dispersion"))
			pick.pars.with.se.cond <- !is.element(names(fitted_par), c(L.names, "dispersion"))
		} else {
		  #pick.pars.with.se.cond <- !is.element(names(fitted_par), c("sigma", L.names))
		  pick.pars.with.se.cond <- !is.element(names(fitted_par), c(L.names))
		}
		par.with.se.cond <- fitted_par[pick.pars.with.se.cond]
		par.fixed.for.se.cond <- par.fixed
		par.fixed.for.se.cond[names(fitted_par[!pick.pars.with.se.cond])] <- fitted_par[!pick.pars.with.se.cond]
		#	par.fixed.for.se["dispersion"] <- fitted_par["dispersion"]

		H.cond <- hessian(func = mnTS.ml.wrapper, x = par.with.se.cond, par.fixed = par.fixed.for.se.cond, Y = Y, X = X, Tsample = Tsample,
		                  method.args = hessian.method.args)

		if (any(is.nan(H.cond)) || rcond(H.cond) < 1e-12) {
		  inv.information.matrix.cond <- matrix(NA, nrow(H.cond), ncol(H.cond))
		} else {
		  inv.information.matrix.cond <- solve(H.cond)
		}
		se.cond <- diag(inv.information.matrix.cond)^0.5
		names(se.cond) <- names(par.with.se.cond)
		rownames(inv.information.matrix.cond) <- names(par.with.se.cond)
		colnames(inv.information.matrix.cond) <- names(par.with.se.cond)

		} else {
		  information.matrix <- NULL
		  information.matrix.cond <- NULL
		  inv.information.matrix <- NULL
		  inv.information.matrix.cond <- NULL
		}

	fitted_par["dispersion"] <- exp(fitted_par["dispersion"])

	if (!is.null(X))
		tX <- t(X)
	else tX <- NULL

	npar <- length(fitted_par)
	AIC <- -2 * fitted$logLik + 2 * npar

	results <- list(Y = t(Y), X = tX, Tsample = Tsample, B0 = t(fitted$B0), C = fitted$C,
		sigma = fitted$sigma, dispersion = fitted$dispersion,
		V = fitted$V, B = t(fitted$B), logLik = fitted$logLik,
		AIC = AIC, npar = npar, inv.information.matrix = inv.information.matrix,
		inv.information.matrix.cond = inv.information.matrix.cond,
		par.with.se = par.with.se, se = se, par.with.se.cond = par.with.se.cond, se.cond = se.cond,
		y = fitted$y, mu = fitted$mu, par = fitted_par, par.fixed = par.fixed,
		opt.convergence = opt.convergence, se.y.fitted = fitted$se.y.fitted,
		se.mu.upper.fitted = fitted$se.mu.upper.fitted, se.mu.lower.fitted = fitted$se.mu.lower.fitted,
		information.matrix.cond = H.cond, information.matrix = H)
	class(results) <- "mnTS"
	return(results)
}

# summary.mnTS --------------------------------------------------------
#' Summary Method for Multinomial State-Space Models
#'
#' Provides a summary of the fitted multinomial state-space model.
#'
#' @param mod An object of class \code{"mnTS"}, typically returned by the
#'     \code{\link{mnTS}} function.
#' @param ... Additional arguments (currently unused).
#'
#' @return Returns a summary list containing key model estimates.
#'     Also prints the summary to the console.
#'
#' @export

summary.mnTS <- function(mod, ...) {

	cat("\nCall: mnTS with Tmax = ", nrow(mod$Y), " n = ", ncol(mod$Y), "\n")

	cat("\nlogLik = ", mod$logLik, ",  AIC = ", mod$AIC, " [df = ", mod$npar, "]\n", sep = "")

	cat("\ndispersion parameter = ", mod$dispersion)

	if (!is.null(mod$inv.information.matrix)) {
	  if (!all(is.na(mod$inv.information.matrix))) {
	    cat("\n\nFitted Coefficients with approximate se\n")
	    v.pars <- names(mod$par)[grep(names(mod$par), pattern = "v.")]
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
	  if (!all(is.na(mod$inv.information.matrix.cond))) {
	    cat("\n\nFitted Coefficients with conditional approximate se\n")
	    t.scores.cond <- mod$par.with.se.cond/mod$se.cond
	    coef.table.cond <- as.matrix(cbind(mod$par.with.se.cond, mod$se.cond, t.scores.cond, 2 * pnorm(q = abs(t.scores.cond), lower.tail = F)))
	    colnames(coef.table.cond) <- c("Coef.", "se", "t", "P")
	    print(coef.table.cond)
	  } else {
	    coef.table.cond <- as.matrix(mod$par.with.se.cond)
	    colnames(coef.table.cond) <- "Coef."
	    print(coef.table.cond)
	    cat("\nconditional information matrix not invertible so se's not calculated")
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
}

print.mnTS <- summary.mnTS



# coef.mnTS -----------------------------------------------------------
#' Coefficients for Multinomial State-Space Model
#'
#' Extracts estimated coefficients from a fitted multinomial state-space model.

#'
#' @param mod An object of class \code{"mnTS"}, as returned by the
#'     \code{\link{mnTS}} function.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list of estimated model coefficients.
#'
#' @export

coef.mnTS <- function(mod, ...) {

  if (!is.null(mod$inv.information.matrix)) {
    if (!all(is.na(mod$inv.information.matrix))) {
       v.pars <- names(mod$par)[grep(names(mod$par), pattern = "v.")]
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
  }
  return(coef.table = coef.table)
}


# simulate.mnTS -------------------------------------------------------
#' Simulate Data from a Multinomial State-Space Model
#'
#' Generates simulated datasets from a fitted mnTS object.
#'
#' @param mod An object of class \code{"mnTS"}, as returned by the
#'     \code{\link{mnTS}} function.
#' @param ... Additional arguments passed to the simulation routine.
#'
#' @return A simulated dataset in the same structure as the response variable
#'     \code{Y} used to fit the model.
#'
#' @importFrom mvtnorm rmvnorm
#' @export


simulate.mnTS <- function(mod, ...) {

  Y <- t(mod$Y)
  n <- nrow(Y)
  Tmax <- nrow(mod$X)

  if (!is.null(mod$X)) {
    X <- t(mod$X)
    # p <- nrow(X)
    B <- t(mod$B)
  } else {
    X <- NULL
    p <- 1
    B = (matrix(0, n, p))
  }

  B0 <- t(mod$B0)
  C <- mod$C
  sigma <- mod$sigma

  # V <- mod$V
  # L <- V
  # S <- t(L) %*% L
  # S <-  mod$V
  # Ssd <- diag(diag(S)^-.5)
  # Se <- sigma^2 * Ssd %*% S %*% Ssd
  # Ssd <- S[1, 1]
  # Se <- sigma^2 * S/Ssd
  Se <- sigma^2 * mod$V

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
    size <- 100

    Y.sim[,t] <- rmultinom(n = 1, size = size, prob = mu)

  }
  return(list(Y = t(Y.sim), X = t(X), B0 = t(B0), C = C, sigma = abs(sigma),
              V = mod$V, B = t(B)))
}


# boot ------------------------------------------------------------
#' Multinomial state-space bootstrap
#'
#' This is a generic bootstrap function for models.
#'
#' @param object An object for which the bootstrap method is defined.
#' @param ... Additional arguments passed to the method.
#'
#' @export
boot <- function(object, ...) {
  UseMethod("boot")
}


# boot.mnTS -------------------------------------------------------
#' Bootstrap method for mnTS objects
#'
#' Provides bootstrapped confidence intervals for an object of class \code{mnTS}.
#'
#' @param object An object of class \code{mnTS}.
#' @param reps Integer. Number of bootstrap replications.
#' @param dispersion.fixed Numeric. Dispersion parameter (default is 1).
#' @param maxit.optim Integer. Max iterations for optimization (default 100000).
#' @param ... Additional arguments (not currently used).
#'
#' @export
#' @method boot mnTS

boot.mnTS <- function(mod, reps, dispersion.fixed = 1, maxit.optim = 1e+05, ...) {

  mods <- replicate(simulate(mod), n = reps)
  fit_mods <- vector(mode = "list", length = ncol(mods))

  for (j in seq_len(ncol(mods))) {

    X <- mods[2, j]$X
    Y <- mods[1, j]$Y
    Y <- Y[mod$Tsample, ]

    if (!is.null(X)) {
      p <- ncol(mods[2, j]$X)
    } else{
      p <- 1
    }

    n <- ncol(mods[1, j]$Y)
    B.fixed <- matrix(NA, p, n)
    B.fixed[,1] <- 0

    V.fixed <- matrix(NA, n, n)
    V.fixed[1] <- 1

    C.fixed <- mod$C
    C.fixed[C.fixed != 0] <- NA

    fit_mods[[j]] <- tryCatch(
      {
        ss_mod <- mnTS(Y = Y, X = X,
                       B0.start = mod$B0,
                       B.start = mod$B, B.fixed = B.fixed,
                       C.fixed = C.fixed, C.start = mod$C*0.5,
                       V.fixed = V.fixed, V.start = mod$V,
                       dispersion.fixed = dispersion.fixed,
                       Tsample = mod$Tsample, maxit.optim = maxit.optim)
      },
      error=function(cond) {
        message("Here's the original error message:")
        message(conditionMessage(cond))
        return(list(err = cond, full_mod = mods[1, j]))
      }
    )
  }
  # Filter out any models that had an error
  fit_mods <- fit_mods[lapply(fit_mods, class) == "mnTS"]
  # Filter out and models that failed to converge
  # fit_mods_filtered <- Filter(function(x) x$opt.convergence == 0, fit_mods)
  # n_converged <- length(fit_mods_filtered)
  # n_total <- length(fit_mods)

  v.pars <- names(mod$par)[grep(names(mod$par), pattern = "v.")]
  par.names.noshow.se <- c("dispersion","sigma", v.pars)

  pars <- lapply(seq_along(fit_mods), \(rep) {
    par.show.se <- fit_mods[[rep]]$par[!is.element(names(fit_mods[[rep]]$par), par.names.noshow.se)]
    par.all <- c(par.show.se, logLik = fit_mods[[rep]]$logLik, opt.convergence = fit_mods[[rep]]$opt.convergence)
  })
  all_pars <- do.call(rbind, pars)
  most_pars <- all_pars[, !colnames(all_pars) %in% c("logLik", "opt.convergence")]

  most_pars_mean <- matrixStats::colMeans2(most_pars, na.rm = T)
  most_pars_sd <- matrixStats::colSds(most_pars, na.rm = T)
  # most_pars_upper_68 <- matrixStats::colQuantiles(most_pars, probs = 0.84, na.rm = T)
  # most_pars_lower_68 <- matrixStats::colQuantiles(most_pars, probs = 0.16, na.rm = T)
  most_pars_upper_68 <- most_pars_mean + most_pars_sd
  most_pars_lower_68 <- most_pars_mean - most_pars_sd
  t_scores <- most_pars_mean / most_pars_sd
  p_vals <- 2 * pnorm(q = abs(t_scores), lower.tail = F)

  n_fail_converged <- rep(sum(all_pars[, colnames(all_pars) %in% "opt.convergence"]), length(most_pars_mean))
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
