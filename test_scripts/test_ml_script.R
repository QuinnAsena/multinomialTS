Rcpp::sourceCpp("./R/source_mnTS_list.cpp")
Rcpp::sourceCpp("./R/source_mnTS.cpp")
library(minqa)

mnTS.ml.wrapper <- function(par, par.fixed, Y, X = NULL, Tsample) {
  retval <- mnTS_ml_cpp_list(par, par.fixed, Y, X, Tsample)
  # print(retval$LL)
  return(retval$LL)
}


opt2 <- bobyqa(fn = mnTS.ml.wrapper, par = par.start, control = optim.control, par.fixed = par.fixed, Y = Y, X = X,
       Tsample = Tsample)

opt <- bobyqa(fn = mnTS_ml_cpp, par = par.start, control = optim.control, par_fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
optlist <- bobyqa(fn = mnTS.ml.wrapper, par = par.start, control = optim.control, par.fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
optr <- bobyqa(fn = mnTS.ml, par = par.start, control = optim.control, par.fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)

test <- mnTS_ml_cpp_list(par = par.start, par_fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
test$LL

mnTS_ml_cpp(par = par.start, par_fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)

MultinomialStateSpace::multinomSS.ml(par = par.start, par.fixed = par.fixed, Y = Y, X = X, Tsample = Tsample,
                                     fitted.values = FALSE, REML = FALSE)

mnTS.ml(par = par.start, par.fixed = par.fixed, Y = Y, X = X, Tsample = Tsample,
                                     fitted.values = FALSE)




test$logLik
testr <- MultinomialStateSpace::multinomSS.ml(par = par.start, par.fixed = par.fixed, Y = Y, X = X, Tsample = 1:ncol(X),
                                              fitted.values = TRUE, REML = FALSE)
testr$logLik





fitted <- mnTS_ml_cpp_list(opt$par, par_fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
fitted$par
fitted_par <- as.vector(fitted$par)
names(fitted_par) <- names(par.fixed)[is.na(par.fixed)]

fittedr <- MultinomialStateSpace::multinomSS.ml(opt$par, par.fixed = par.fixed, Y = Y, X = X, Tsample = Tsample, fitted.values = T)
fittedr$par

fitted







# Create data
X_pulse <- c(rep(0, 30),
             rep(1, 30),
             rep(0, 30),
             rep(1, 30),
             rep(0, 30))
X_lin <- 1:150
X_lin <- (X_lin - mean(X_lin)) / sd(X_lin)

X0s <- rep(0, 150)
X1s <- rep(1, 150)

X <- t(cbind(X_pulse, X_lin, X1s, X0s))

Y <- sim_prox_func(n = 3, size = 100, Tmax = ncol(X),
                   Tstart = 0, X = X,
                   sigma = 0.2,
                   C = matrix(c(.5, 0, 0, 0, .9, 0, 0, -.3, .6), 3, 3),
                   B = NULL,
                   B0 = c(0, 0.25, -0.5),
                   V = matrix(c(1, 0.7, 0.7, 0.7, .6, .5, 0.7, .5, .6), 3, 3),
                   print_plot = T, seed = 1985)

X <- X[, 11:ncol(X), drop = F]
Y <- t(Y$Y[11:nrow(Y$Y), ])
p <- nrow(X) + 1
n <- nrow(Y)
# for uneven Tsample
Tsample <- c(1, 3, 9, 11, 14, 16, 18, 20, 21, 23, 25, 27, 28, 30, 32, 38,
             43, 48, 53, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 70, 71, 72,
             74, 76, 77, 78, 80, 82, 84, 86, 87, 88, 89, 90, 92, 93, 94, 95,
             96, 97, 100, 101, 102, 103, 104, 106, 107, 108, 109, 111, 115,
             119, 123, 126, 131, 134, 135, 137, 138)

## from cpp version
# dput(ss_cpp$par.fixed)
par.fixed <- c(y1 = 0, y2 = NA, y3 = NA, sp.y1.y1 = NA, sp.y2.y1 = 0, sp.y3.y1 = 0,
               sp.y1.y2 = 0, sp.y2.y2 = NA, sp.y3.y2 = NA, sp.y1.y3 = 0, sp.y2.y3 = NA,
               sp.y3.y3 = NA, sigma = NA, v.y1.y1 = 1, v.y2.y1 = 0, v.y3.y1 = 0,
               v.y1.y2 = NA, v.y2.y2 = NA, v.y3.y2 = 0, v.y1.y3 = NA, v.y2.y3 = NA,
               v.y3.y3 = NA, dispersion = 0, X_pulse.y1 = 0, X_lin.y1 = 0, X1s.y1 = 0,
               X0s.y1 = 0, X_pulse.y2 = NA, X_lin.y2 = NA, X1s.y2 = NA, X0s.y2 = NA,
               X_pulse.y3 = NA, X_lin.y3 = NA, X1s.y3 = NA, X0s.y3 = NA)

# dput(ss_cpp$par)
par.start <- c(y2 = -0.0356968674797618, y3 = -0.566547638732945, sp.y1.y1 = 0.246256632166609,
               sp.y2.y2 = 2.20070891772196, sp.y3.y2 = -0.699022902665956, sp.y2.y3 = 2.36270525427544,
               sp.y3.y3 = -0.343863395834362, sigma = 0.165357169702909, v.y1.y2 = -0.226761860226737,
               v.y2.y2 = -5.50590520877739e-07, v.y1.y3 = 0.141346565904594,
               v.y2.y3 = 3.09111351061685e-07, v.y3.y3 = -9.59415191068563e-10,
               X_pulse.y2 = -0.104376816254213, X_lin.y2 = 0.0235506828140919,
               X1s.y2 = -0.0672722619592519, X0s.y2 = 0.11114413652046, X_pulse.y3 = 0.0127573089763948,
               X_lin.y3 = -0.0970258339740162, X1s.y3 = -0.241755157674573,
               X0s.y3 = 0.139157671603634)



# Run the c++ ml
Rcpp::sourceCpp("./R/source_mnTS_list.cpp")
test <- mnTS_ml_cpp_list(par = par.start, par_fixed = par.fixed, Y = Y, X = X, Tsample = 1:ncol(X))
mnTS_ml_cpp(par = par.start, par_fixed = par.fixed, Y = Y, X = X, Tsample = 1:ncol(X))












## Setting up some parameters
par.start <- c(y2 = 0.229649109483497, y3 = -0.480376214942664, sp.y1.y1 = 0.5,
               sp.y2.y2 = 0.5, sp.y3.y2 = 0.01, sp.y2.y3 = 0.01, sp.y3.y3 = 0.5,
               sigma = 5.84431822064096, v.y1.y2 = 0, v.y2.y2 = 0.429198404113602,
               v.y1.y3 = 0, v.y2.y3 = 0, v.y3.y3 = 0.563289409458386, x1.y2 = 0.0747906057996346,
               x1.y3 = -0.0281864307293143)
par.fixed <- c(other = 0, Pinus = NA, Tsuga = NA, sp.other.other = NA, sp.Pinus.other = 0,
               sp.Tsuga.other = 0, sp.other.Pinus = 0, sp.Pinus.Pinus = NA,
               sp.Tsuga.Pinus = NA, sp.other.Tsuga = 0, sp.Pinus.Tsuga = NA,
               sp.Tsuga.Tsuga = NA, sigma = NA, v.other.other = 1, v.Pinus.other = 0,
               v.Tsuga.other = 0, v.other.Pinus = NA, v.Pinus.Pinus = NA, v.Tsuga.Pinus = 0,
               v.other.Tsuga = NA, v.Pinus.Tsuga = NA, v.Tsuga.Tsuga = NA, dispersion = 0,
               x1.other = 0, x1.Pinus = NA, x1.Tsuga = NA)

# Rcpp::sourceCpp("./R/source_mnTS_list.cpp")
# Rcpp::sourceCpp("./R/source_mnTS.cpp")

# Create data
X_pulse <- c(rep(0, 30),
             rep(1, 30),
             rep(0, 30),
             rep(1, 30),
             rep(0, 30))
X <- t(X_pulse)

# Y <- sim_prox_func(n = 3, size = 100, Tmax = ncol(X),
#                    Tstart = 0, X = X,
#                    sigma = 0.2,
#                    C = matrix(c(.5, 0, 0, 0, .9, 0, 0, -.3, .6), 3, 3),
#                    B = matrix(rep(0, 3), nrow = 3, ncol = 1),
#                    B0 = c(0, 0.25, -0.5),
#                    V = matrix(c(1, 0.7, 0.7, 0.7, .6, .5, 0.7, .5, .6), 3, 3),
#                    print_plot = T, seed = 1984)

X <- X[, 11:ncol(X), drop = F]
# Y <- t(Y$Y[11:nrow(Y$Y), ])
Y <- structure(c(40L, 29L, 31L, 40L, 32L, 28L, 42L, 35L, 23L, 46L,
                 33L, 21L, 41L, 39L, 20L, 39L, 30L, 31L, 40L, 25L, 35L, 45L, 26L,
                 29L, 39L, 41L, 20L, 41L, 36L, 23L, 46L, 32L, 22L, 39L, 35L, 26L,
                 39L, 33L, 28L, 44L, 39L, 17L, 41L, 39L, 20L, 32L, 44L, 24L, 36L,
                 40L, 24L, 45L, 27L, 28L, 36L, 43L, 21L, 38L, 43L, 19L, 29L, 53L,
                 18L, 42L, 48L, 10L, 34L, 51L, 15L, 30L, 52L, 18L, 32L, 41L, 27L,
                 28L, 43L, 29L, 30L, 52L, 18L, 26L, 54L, 20L, 39L, 36L, 25L, 29L,
                 46L, 25L, 30L, 46L, 24L, 30L, 48L, 22L, 31L, 46L, 23L, 40L, 40L,
                 20L, 33L, 44L, 23L, 37L, 44L, 19L, 37L, 45L, 18L, 34L, 55L, 11L,
                 39L, 37L, 24L, 37L, 48L, 15L, 35L, 40L, 25L, 37L, 35L, 28L, 32L,
                 38L, 30L, 43L, 41L, 16L, 34L, 41L, 25L, 39L, 39L, 22L, 34L, 40L,
                 26L, 36L, 46L, 18L, 42L, 44L, 14L, 35L, 48L, 17L, 36L, 44L, 20L,
                 30L, 51L, 19L, 39L, 46L, 15L, 34L, 38L, 28L, 32L, 42L, 26L, 37L,
                 49L, 14L, 38L, 42L, 20L, 36L, 47L, 17L, 41L, 39L, 20L, 32L, 51L,
                 17L, 31L, 47L, 22L, 33L, 48L, 19L, 35L, 39L, 26L, 34L, 41L, 25L,
                 39L, 36L, 25L, 31L, 43L, 26L, 38L, 42L, 20L, 42L, 34L, 24L, 34L,
                 53L, 13L, 37L, 47L, 16L, 27L, 47L, 26L, 33L, 52L, 15L, 39L, 43L,
                 18L, 36L, 45L, 19L, 26L, 51L, 23L, 25L, 58L, 17L, 32L, 48L, 20L,
                 28L, 49L, 23L, 34L, 46L, 20L, 28L, 52L, 20L, 32L, 48L, 20L, 25L,
                 43L, 32L, 36L, 47L, 17L, 35L, 48L, 17L, 36L, 49L, 15L, 38L, 44L,
                 18L, 30L, 52L, 18L, 36L, 49L, 15L, 34L, 54L, 12L, 37L, 44L, 19L,
                 32L, 45L, 23L, 32L, 45L, 23L, 31L, 51L, 18L, 22L, 60L, 18L, 33L,
                 51L, 16L, 32L, 54L, 14L, 40L, 37L, 23L, 34L, 51L, 15L, 27L, 55L,
                 18L, 27L, 55L, 18L, 32L, 44L, 24L, 39L, 45L, 16L, 30L, 53L, 17L,
                 33L, 47L, 20L, 26L, 48L, 26L, 38L, 30L, 32L, 46L, 35L, 19L, 41L,
                 35L, 24L, 32L, 47L, 21L, 33L, 44L, 23L, 33L, 32L, 35L, 34L, 44L,
                 22L, 37L, 37L, 26L, 32L, 46L, 22L, 37L, 40L, 23L, 38L, 47L, 15L,
                 36L, 46L, 18L, 27L, 51L, 22L, 31L, 42L, 27L, 30L, 46L, 24L, 29L,
                 54L, 17L, 29L, 51L, 20L, 37L, 40L, 23L, 35L, 47L, 18L, 29L, 53L,
                 18L, 26L, 50L, 24L, 28L, 51L, 21L, 28L, 58L, 14L, 28L, 53L, 19L,
                 29L, 53L, 18L, 29L, 51L, 20L, 21L, 50L, 29L, 36L, 45L, 19L, 29L,
                 53L, 18L, 32L, 48L, 20L, 39L, 42L, 19L, 31L, 53L, 16L, 35L, 43L,
                 22L, 33L, 44L, 23L, 36L, 51L, 13L), dim = c(3L, 140L))
# for uneven Tsample
Tsample <- c(1, 3, 9, 11, 14, 16, 18, 20, 21, 23, 25, 27, 28, 30, 32, 38,
             43, 48, 53, 58, 59, 60, 61, 63, 64, 65, 66, 67, 68, 70, 71, 72,
             74, 76, 77, 78, 80, 82, 84, 86, 87, 88, 89, 90, 92, 93, 94, 95,
             96, 97, 100, 101, 102, 103, 104, 106, 107, 108, 109, 111, 115,
             119, 123, 126, 131, 134, 135, 137, 138)

# This is where things happen
maxit.optim = 1e+06
# If I limit maxit.optim to 5 iterations and include a print statement for LL/2 in the R and cpp source code the first few printed values and overall results are the same
# With 1e+06 max.optim, things diverge at some point
optim.control = list(maxfun = maxit.optim)

# Run the c++ ml
# Rcpp::sourceCpp("./R/source_mnTS_list.cpp")

mnTS_ml_cpp_listout(par = par.start, par_fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
# mnTS_ml_cpp_list(par = par.start, par_fixed = par.fixed, Y = Y, X = X, Tsample = 1:ncol(X))
test <- mnTS_ml_cpp_list(par = par.start, par_fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
test$LL

mnTS_ml_cpp(par = par.start, par_fixed = par.fixed, Y = Y, X = X, Tsample = Tsample)
MultinomialStateSpace::multinomSS.ml(par = par.start, par.fixed = par.fixed, Y = Y, X = X, Tsample = Tsample,
                                     fitted.values = FALSE, REML = FALSE)


test$logLik
testr <- MultinomialStateSpace::multinomSS.ml(par = par.start, par.fixed = par.fixed, Y = Y, X = X, Tsample = 1:ncol(X),
                                              fitted.values = TRUE, REML = FALSE)
testr$logLik


mnTS_ml_cpp(par = par.start, par_fixed = par.fixed, Y = matrix(0, nrow = 3, ncol = 140), X = X, Tsample = 1:ncol(X))
mnTS_ml_cpp_list(par = par.start, par_fixed = par.fixed, Y = matrix(0, nrow = 3, ncol = 140), X = X, Tsample = 1:ncol(X))
