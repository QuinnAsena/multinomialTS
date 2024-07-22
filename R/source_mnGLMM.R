# Multinomial GLMM ML -----------------------------------------------------
#' Multinomial GLMM Maximum likelihood
#'
#' This is a simple function that, by default, prints "Hello world". You can
#' customize the text to prin
#'
#' @param par something
#' @param par.fixed More things
#' @param Y inout matrix of
#' @param X input matrix of
#' @param fitted.values TRUE/FALSE argument indicating
#' @param REML Description
#'
#' @return This function returns a phrase to print, with or without an
#'    exclamation point added. As a side effect, this function also prints out
#'    the phrase.
#'
#'
#' @import minqa
#' @import numDeriv
#' @import stats
#'
#' @export

mnGLMM.ml <- function(par, par.fixed, Y, X, fitted.values = FALSE, REML = FALSE) {
	parameters <- par.fixed
	parameters[is.na(par.fixed)] <- par

	n <- ncol(Y)
	nn <- n - 1
	p <- ncol(X)
	Tmax <- nrow(Y)
	size <- rowSums(Y)

	B <- matrix(parameters[1:(n * p)], p, n, byrow = F)
	L <- matrix(parameters[(n * p + 1):(n * p + n * n)], n, n)
	sigma <- parameters[n * p + n * n + 1]

	S <- t(L) %*% L
	#Ssd <- diag(diag(S)^-.5)
	#P <- sigma^2 * Ssd %*% S %*% Ssd
	Ssd <- S[1, 1]
	P <- sigma^2 * S/Ssd

	dispersion <- exp(parameters[n * p + n * n + 1 + 1])

	y <- X %*% B
	mu <- exp(y)/((exp(y) %*% matrix(1, n, 1)) %*% matrix(1, 1, n))

	SS <- 0
	for (i in 1:nrow(mu)) {

		W <- diag((mu[i, ] * (1 - mu[i, ]))^-1) %*% (diag(mu[i, ]) - matrix(mu[i, ], n, 1) %*% matrix(mu[i, ], 1, n)) %*% diag((mu[i,
			] * (1 - mu[i, ]))^-1)
		M <- diag(1 + mu[i, ]) - (matrix(1, n, 1) %*% matrix(mu[i, ], 1, n))
		S <- dispersion * (W/size[i]) + M %*% P %*% t(M)
		D <- diag(mu[i, ] * (1 - mu[i, ]))
		S <- D %*% S %*% D

		v <- (Y[i, -1]/size[i] - mu[i, -1])

		if (nn > 1) {
			invS <- solve(S[-1, -1])
			if (REML) {
				SS <- SS + t(v) %*% invS %*% v + as.numeric(determinant(S[-1, -1])$modulus) + as.numeric(determinant(t(matrix(1, nn, p) %*%
					X[i, -1]) %*% invS %*% (matrix(1, nn, p) %*% X[i, -1]))$modulus)
			} else {
				SS <- SS + t(v) %*% invS %*% v + as.numeric(determinant(S[-1, -1])$modulus)
			}
		} else {
			if (REML) {
				SS <- SS + v^2/S[-1, -1] + log(S[-1, -1]) + X[i, -1]^2/S[-1, -1]
			} else {
				SS <- SS + v^2/S[-1, -1] + log(S[-1, -1])
			}
		}
	}
if (!fitted.values) {
		return(SS/2)
	} else {
		colnames(B) <- colnames(Y)
		rownames(B) <- colnames(X)

		V <- (t(L) %*% L)
		colnames(V) <- colnames(Y)
		rownames(V) <- colnames(Y)

		par <- parameters[is.na(par.fixed)]
		logLik <- -((n - 1) * Tmax/2) * log(2 * pi) - SS/2

		return(list(B = B, sigma = abs(sigma), V = V, dispersion = dispersion, logLik = logLik, par = par, mu = mu))
	}
}


# Multinomial GLMM --------------------------------------------------------
#' Multinomial GLMM
#'
#' This is a not a simple function
#' customize the text to print.
#'
#' @param Y inout matrix of
#' @param X input matrix of
#' @param fitted.values TRUE/FALSE argument indicating
#' @param REML Description
#'
#' @return This function returns a phrase to print, with or without an
#'    exclamation point added. As a side effect, this function also prints out
#'    the phrase.
#'
#'
#' @import minqa
#'
#' @export

mnGLMM <- function(Y, X = NULL, B.fixed = if (is.null(X))
	matrix(c(0, rep(NA, ncol(Y) - 1)), nrow = 1, ncol = ncol(Y))
else matrix(c(rep(0, (ncol(X) + 1) * ncol(Y)), rep(NA, (ncol(Y) - 1) * (ncol(X) + 1))), nrow = ncol(X) + 1, ncol = ncol(Y)), B.start = if (is.null(X))
	matrix(0, nrow = 1, ncol = ncol(Y))
else matrix(0, nrow = ncol(X) + 1, ncol = ncol(Y)), sigma.fixed = NA, sigma.start = 0.1, dispersion.fixed = 1, dispersion.start = 1,
	V.fixed = diag(ncol(Y)), V.start = diag(ncol(Y)), method = "bobyqa", optim.control = NULL, maxit.optim = 1e+05, REML = FALSE, compute.information.matrix = TRUE) {

	# catch cases caused when X is defined outside the function call
	# if(is.null(X)) {
# B.fixed <- matrix(c(0, rep(NA, ncol(Y) - 1)), nrow = 1, ncol = ncol(Y))
# B.start <- matrix(0, nrow = 1, ncol = ncol(Y))
# }

	# check and name Y
	if (!is.matrix(Y)) {
		Y <- as.matrix(Y)
	}
	if (!is.null(X))
		if (!(ncol(X) == (dim(B.fixed)[1] - 1) & ncol(Y) == dim(B.fixed)[2])) {
			stop("dimensions of B.fixed should be c(ncol(X)+1, ncol(Y)).")
		}
	if (is.null(colnames(Y))) {
		colnames(Y) <- rep(NA, ncol(Y))
		for (i.col in 1:ncol(Y)) colnames(Y)[i.col] <- paste0("y", i.col)
	}
	# if (any(diag(var(Y, na.rm = TRUE)) == 0)) {
	# stop("The response Y (dependent variable) has no variation.")
# }
n <- ncol(Y)
	nn <- n - 1
	Tmax <- nrow(X)

	# # V
	# if(sum(diag(V)) != n) stop("Diagonal elements of V must be one.")
#
# # expand sigma
# if(length(sigma.fixed) == 1){
#   sigma.fixed <- rep(sigma.fixed, n)
# }
# if(length(sigma.start) == 1){
#   sigma.start <- rep(sigma.start, n)
# }

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

	# log-transform dispersion
	dispersion.fixed <- log(dispersion.fixed)
	dispersion.start <- log(dispersion.start)

	par.fixed <- c(B.fixed, L.fixed, sigma.fixed, dispersion.fixed)

	B.names <- rep(NA, n * p)
	for (j in 1:n) for (i in 1:p) B.names[i + p * (j - 1)] <- paste0(colnames(X)[i], ".", colnames(Y)[j])
	L.names <- rep(NA, n * n)
	for (j in 1:n) for (i in 1:n) L.names[i + n * (j - 1)] <- paste0("v.", colnames(Y)[i], ".", colnames(Y)[j])

	#if(is.na(dispersion.fixed)) dispersion.name <- "dispersion" else dispersion.name <- NULL
	#names(par.fixed) <- c(B.names, L.names, "sigma", dispersion.name)

	names(par.fixed) <- c(B.names, L.names, "sigma", "dispersion")

	par.start <- c(B.start, L.start, sigma.start, dispersion.start)
	names(par.start) <- c(B.names, L.names, "sigma", "dispersion")
	par.start <- par.start[is.na(par.fixed)]

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
		opt <- optim(fn = mnGLMM.ml, par = par.start, method = "Nelder-Mead", control = optim.control, par.fixed = par.fixed, Y = Y,
			X = X, fitted.values = FALSE, REML = FALSE)
		if (opt$convergence != 0)
			cat("\nNelder-Mead optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$convergence
	}

	if (method == "BFGS") {
		opt <- optim(fn = mnGLMM.ml, par = par.start, method = "BFGS", control = optim.control, par.fixed = par.fixed, Y = Y, X = X,
			fitted.values = FALSE, REML = FALSE)
		if (opt$convergence != 0)
			cat("\nBFGS optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$convergence
	}

	if (method == "bobyqa") {
		opt <- bobyqa(fn = mnGLMM.ml, par = par.start, control = optim.control, par.fixed = par.fixed, Y = Y, X = X, fitted.values = FALSE,
			REML = FALSE)
		if (opt$ierr != 0)
			cat("\nbobyqa optimization failed to converge; check $opt.convergence\n")
		opt.convergence <- opt$ierr
		opt$value <- opt$fval
	}

	fitted <- mnGLMM.ml(opt$par, par.fixed = par.fixed, Y = Y, X = X, fitted.values = TRUE, REML = FALSE)

	names(fitted$par) <- names(par.fixed)[is.na(par.fixed)]

	if (compute.information.matrix) {

		par.with.se <- fitted$par
		par.fixed.for.se <- par.fixed

		H <- numDeriv::hessian(func = mnGLMM.ml, x = par.with.se, par.fixed = par.fixed.for.se, Y = Y, X = X, fitted.values = FALSE,
			REML = FALSE, method.args = list(eps = 1e-04, d = 0.1))
		if (any(is.nan(H)) || rcond(H) < 1e-12) {
			inv.information.matrix <- matrix(NA, nrow(H), ncol(H))
		} else {
			inv.information.matrix <- solve(H)
		}

		# attenpt to use a robust SE estimate
		# gradL <- NULL
		# for(i in 1:length(fitted$par)) {
		# grad.par <- fitted$par
		# grad.par[i] <- 0
		# gradL <- c(gradL, (fitted$LL + mnGLMM.ml(grad.par, par.fixed = par.fixed, Y = Y, X = X))/fitted$par[i])
		# }
		# gradL <- matrix(gradL, nrow = 1)
		# robust.cov <- inv.information.matrix %*% (t(gradL) %*% gradL) %*% inv.information.matrix
		} else {
		information.matrix <- NULL
		}
	# logLik <- logLik.glm +
	# as.numeric(-LL + pglmm.LL(0 * ss, H = H, X = X, Zt = Zt, St = St, mu = mu,
	# nested = nested, REML = REML, family = family,
	# size = size, verbose = verbose))

	fitted$par["dispersion"] <- exp(fitted$par["dispersion"])

	npar <- length(fitted$par)
	AIC <- -2 * fitted$logLik + 2 * npar

	results <- list(Y = Y, X = X, B = fitted$B, sigma = fitted$sigma, V = fitted$V, dispersion = fitted$dispersion, logLik = fitted$logLik,
		AIC = AIC, npar = npar, inv.information.matrix = inv.information.matrix, par.with.se = par.with.se, se = diag(inv.information.matrix)^0.5,
		par = fitted$par, par.fixed = par.fixed, opt.convergence = opt.convergence)
	class(results) <- "mnGLMM"
	return(results)
}

# summary.mnGLMM --------------------------------------------------------
#' GLMM summary
#'
#' This is a simple function
#' customize the text to
#'
#' @param mod inout matrix of
#' @param ... input matrix of
#'
#' @return This function returns a phrase to print, with or without an
#'    exclamation point added. As a side effect, this function also prints out
#'    the phrase.
#'
#'
#' @export

summary.mnGLMM <- function(mod, ...) {

	cat("\nCall: mnGLMM with Tmax = ", nrow(mod$Y), " n = ", ncol(mod$Y), "\n")

	# npar <- length(mod$par)
	# AIC <- -2*mod$logLik + 2*npar
	cat("\nlogLik = ", mod$logLik, ",  AIC = ", mod$AIC, " [df = ", mod$npar, "]\n", sep = "")

	cat("\ndispersion parameter = ", mod$dispersion)

	if (!is.null(mod$inv.information.matrix)) {
		cat("\n\nFitted Coefficients with approximate se\n")
		if (!all(is.na(mod$inv.information.matrix))) {
			v.pars <- names(mod$par)[grep(names(mod$par), pattern = "v.")]
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
}

print.mnGLMM <- summary.mnGLMM
