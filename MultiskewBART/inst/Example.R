## Load ----

library(MultiskewBART)
library(tidyverse)
library(zeallot)
library(MCMCpack)

## Function for generating data ----

sim_data_multi <- function(N, P, lambda, tau, rho) {
  X <- matrix(runif(N * P), nrow = N)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
  Z <- cbind(lambda[1] * abs(rnorm(N)), lambda[2] * abs(rnorm(N)))
  Sigma <- matrix(c(tau[1], sqrt(tau[1]*tau[2])*rho, sqrt(tau[1]*tau[2])*rho, tau[2]), 2, 2)
  Err <- MASS::mvrnorm(n=N, mu=c(0,0), Sigma = Sigma)
  Y <- cbind(mu, mu) + Z + Err
  EY <- rbind(mu, mu) + lambda * sqrt(2/pi)
  return( list(X = X, Y = Y, EY=EY, mu = mu, lambda = lambda, tau=tau, Z= Z, Sigma = Sigma) )
}

## Simulate dataset ----

set.seed(12345)

c(X,Y,EY,mu,lambda,tau,Z,Sigma) %<-%
  sim_data_multi(250, 5, c(2,3), c(1,1), 0.5)

c(test_X,test_Y,test_EY,test_mu,test_lambda,test_tau,test_Z,test_Sigma) %<-%
  sim_data_multi(100, 5, c(2,3), c(1,1), 0.5)

## Fit the model ----

set.seed(23948)

hypers <- Hypers(X = X, Y = Y)
opts   <- Opts(num_burn = 2500, num_save = 2500, num_print = 10)

# debug(MultiskewBART)

fitted_mskew <- MultiskewBART(X = X, Y = Y, test_X = test_X, hypers = hypers,
                              opts = opts)
MultiskewBART:::dmsn(y = Y, mu = fitted_mskew$mu[,,1],
                     Sigma = fitted_mskew$Sigma[,,1],
                     lambda = fitted_mskew$lambda[1,], give_log = TRUE)

