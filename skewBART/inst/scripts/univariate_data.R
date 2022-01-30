## Create a function for Friedman example ----

library(zeallot)

sim_fried <- function(N, P, alpha, sigma) {
  lambda <- alpha * sigma/sqrt(1+alpha^2)
  tau <- sigma/sqrt(1+alpha^2)
  X <- matrix(runif(N * P), nrow = N)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
  Z <- abs(rnorm(N, mean=0, sd=1) )
  Y <- mu + lambda * Z + rnorm(N, mean=0, sd=sqrt(tau))
  EY <- mu + lambda * sqrt(2/pi)
  return(list(X = X, Y = Y, EY = EY, mu = mu, Z=Z, tau = tau, lambda = lambda))
}

## Traning dataset : n = 250 observations, P = 5 covariates, sigma = 2, alpha = 5 ----

set.seed(12345)
c(X,Y,EY,mu,Z,tau,lambda) %<-% sim_fried(250, 5, 5, 2)

## Test dataset : n = 100 observations, P = 5 covariates, sigma = 2, alpha = 5 ----

c(test_X,test_Y,test_EY,test_mu,test_Z,test_tau,test_lambda)  %<-% sim_fried(100, 5, 5 ,2)

## Fit ----

hypers <- Hypers(X, Y)
opts <- Opts(num_burn = 10000, num_save = 10000)
fitted_skewbart <- skewBART(X, Y, test_X, hypers, opts)

## Traceplot of alpha samples and assessment of how well we recover the nonparametric function ----

par(mfrow = c(1,2))
plot(fitted_skewbart$alpha)
plot(colMeans(fitted_skewbart$y_hat_test), test_mu, pch = 2)
abline(a=0,b=1, col = 'green', lwd = 3)
