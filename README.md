README
================



# Bayesian Additive Regression Trees for Multivariate skewed responses

This repository contains code to implement the methodology described in
the paper “Bayesian Additive Regression Trees for Multivariate Skewed
Responses”, by Um et al (2021+).

Note that this package uses the primary functions from
[`SoftBART`](https://github.com/theodds/SoftBART) to incorporate the
SoftBART model as a component.

## Installation

This repository includes two packages: `skewBART` and `MultiskewBART`.
The packages can be installed with the `devtools` package:

``` r
library(devtools) 
devtools::install_github(repo='Seungha-Um/skewBART', subdir='/skewBART') 
devtools::install_github(repo='Seungha-Um/skewBART', subdir='/MultiskewBART', force = TRUE) 
```

## Documentation

Vignette is available at

## Status of developement

The code provided here is being actively developed and only being
provided for research purposes.

## Usage

#### skewBART model (univariate response)

The `skewBART` packages provides a nonparametric regression approach for
univariate skewed responses using Bayesian additive regression trees
(BART).

The following is an example on “Friedman’s example”:

``` r
## Load library
library(skewBART)
library(zeallot)

## Create a function for Friedman’s example
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

## Simulate dataset
## Traning dataset : n = 250 observations, P = 5 covariates, sigma = 2, alpha = 5 
set.seed(12345)
c(X,Y,EY,mu,Z,tau,lambda) %<-% sim_fried(250, 5, 5, 2)

## Test dataset : n = 100 observations, P = 5 covariates, sigma = 2, alpha = 5
c(test_X,test_Y,test_EY,test_mu,test_Z,test_tau,test_lambda)  %<-% sim_fried(100, 5, 5 ,2)
```

The function `Hypers` produces a list of the hyperparameters of the
model.

``` r
hypers <- Hypers(X, Y)
opts <- Opts(num_burn = 5000, num_save = 5000)
```

We can fit the model as

``` r
## Fit the model
fitted_skewbart <- skewBART(X, Y, test_X, hypers, opts)

## Compare the quality of the fit
rmse <- function(x,y) sqrt(mean((x-y)^2))
rmse(fitted_skewbart$y_hat_test_mean, test_mu)
```

#### MultiskewBART model (multivariate response)

The `MultiskewBART` packages provides a nonparametric regression
approach for multivariate skewed responses using Bayesian additive
regression trees (BART).

``` r
## Load library
library(MultiskewBART)
library(MASS)
library(zeallot)

## Create a function for Friedman example with multivariate framework
sim_data_multi <- function(N, P, lambda, tau, rho) {
  X <- matrix(runif(N * P), nrow = N)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5] 
  Z <- cbind(lambda[1] * abs(rnorm(N)), lambda[2] * abs(rnorm(N)))
  Sigma <- matrix(c(tau[1], sqrt(tau[1]*tau[2])*rho, sqrt(tau[1]*tau[2])*rho, tau[2]), 2, 2)
  Err <- mvrnorm(n=N, mu=c(0,0), Sigma = Sigma)
  Y <- cbind(mu, mu) + Z + Err
  EY <- rbind(mu, mu) + lambda * sqrt(2/pi)
  return( list(X = X, Y = Y, EY=EY, mu = mu, lambda = lambda, tau=tau, Z= Z, Sigma = Sigma) )
}

## Simulate dataset
## Traning dataset : n = 250 observations, P = 5 covariates, lambda = (2,3), tau = c(1,1), rho = 0.5.
set.seed(12345)
c(X,Y,EY,mu,lambda,tau,Z,Sigma) %<-% sim_data_multi(250, 5, c(2,3), c(1,1), 0.5)

## Test dataset : n = 250 observations, P = 5 covariates, lambda = (2,3), tau = c(1,1), rho = 0.5.
c(test_X,test_Y,test_EY,test_mu,test_lambda,test_tau,test_Z,test_Sigma) %<-% sim_data_multi(100, 5, c(2,3), c(1,1), 0.5)

## Create a list of the hyperparameters of the model. 
hypers <- Hypers(X = X, Y = Y)
opts <- Opts(num_burn = 5000, num_save = 5000)
```

We can fit the model as

``` r
## Fit the model
fitted_Multiskewbart <- MultiskewBART(X = X, Y = Y, test_X = test_X, hypers=hypers, opts=opts) 
## Compare the quality of the fit
rmse <- function(x,y) sqrt(mean((x-y)^2))
rmse(fitted_Multiskewbart$y_hat_test, t(test_EY))
```
