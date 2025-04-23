story_pollen_matrix <- readRDS("./data/story_pollen_matrix.rds")
story_char_matrix <- readRDS("./data/story_char_matrix.rds")
Tsample <- which(rowSums(story_pollen_matrix) != 0)

story_pollen_matrix[Tsample, ]

library(multinomialTS)
# source("../multinomialTS/cpp_testing_groud/simulate_func_TEMP.R")
# library(Rcpp)
# library(RcppArmadillo)
# library(minqa)
# sourceCpp("../multinomialTS/R/source_mnTS.cpp")


X <- scale(story_char_matrix)
Y <- story_pollen_matrix
p <- ncol(X) + 1 # Number of independent variables plus intercept
n <- ncol(Y)

V.fixed = diag(n) # Covariance matrix of environmental variation in process eq
B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

glmm_mod <- multinomialTS::mnGLMM(Y = Y[Tsample, ],
                            X = X[Tsample, ,drop = F],
                            B.start = B.start, B.fixed = B.fixed,
                            V.fixed = V.fixed)
summary(glmm_mod)

B0.start <- glmm_mod$B[1, , drop = F]
B.start <- glmm_mod$B[2, , drop = F]

sigma.start <- glmm_mod$sigma

V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
V.fixed[1] = 1

V.start = V.fixed
V.start <- glmm_mod$V
# V.start <- diag(diag(V.start))

B.fixed <- matrix(NA, ncol(X), n)
B.fixed[,1] <- 0
B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)


C.start = .5 * diag(n)
C.fixed <- C.start
C.fixed[C.fixed != 0] <- NA


ss_mod <- mnTS(Y = Y[Tsample, ], X = X, Tsample = Tsample, B0.start = B0.start, B.start = B.start,
                            C.start = C.start, C.fixed = C.fixed, B0.fixed = B0.fixed,
                            V.fixed = V.fixed, V.start = V.start,
                            B.fixed = B.fixed, dispersion.fixed = 1, maxit.optim = 1e+06)
summary(ss_mod)

simulate(ss_mod)
# multinomialTS::boot.mnTS(ss_mod, 3)
coef(ss_mod)

boot(ss_mod, reps = 3)

