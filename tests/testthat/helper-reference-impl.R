# Reference implementations of the two likelihoods, in pure R.
#
# These were the original implementations. The package now evaluates both
# likelihoods in C++ (src/source_mnTS_list.cpp and src/source_mnGLMM.cpp), so
# they are no longer used at fit time. They are kept here as the reference the
# equivalency tests check the C++ against -- see test-cpp-equivalency.R.
#
# Keep them behaviourally identical to the C++, including the guard values
# returned for a non-invertible or non-finite likelihood.
#
# One known divergence, in the fitted.values = TRUE branch only: these functions
# apply abs() to the "sigma" element of the returned `par` vector, whereas the
# C++ returns it unmodified (it only applies abs() to the separate `sigma`
# field). sigma enters the likelihood squared, so the optimiser can settle on a
# negative value and the two then disagree in sign on par["sigma"]. The
# equivalency tests compare the likelihood and the B/mu/logLik outputs, not
# `par`, so this is not exercised.

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


	# Mirrors the C++ guard: `if (!std::isfinite(LL)) return 1e+10`. The original
	# R code tested is.complex(LL) here, which none of these operations can
	# produce, so the guard never fired.
	if (!is.finite(LL))
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


mnGLMM.ml <- function(par, par.fixed, Y, X, fitted.values = FALSE) {
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
			SS <- SS + t(v) %*% invS %*% v + as.numeric(determinant(S[-1, -1])$modulus)
		} else {
			SS <- SS + v^2/S[-1, -1] + log(S[-1, -1])
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

