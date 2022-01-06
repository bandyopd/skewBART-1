library(skewBART)
## Create a function for Friedman example

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
## Traning dataset : n = 250 observations, P = 5 covariates, tau = 5 and lambda = 1.

set.seed(12345)
c(X,Y,EY,mu,Z,tau,lambda) %<-% sim_fried(250, 5, 5, 1)

## Test dataset : n = 100 observations, P = 5 covariates, tau = 5 and lambda = 1.

c(test_X,test_Y,test_EY,test_mu,test_Z,test_tau,test_lambda)  %<-% sim_fried(100, 5, 5 ,1)

## preprocess data 

X_norm <- quantile_normalize(X, test_X)$X
test_X_norm <- quantile_normalize(X, test_X)$test_X
Y_scaled <- scale(Y)
c(Y_dims, center_Y, scale_Y) %<-% attributes(Y_scaled)

## Create a list of the hyperparameters of the model. 
hypers <- Hypers(X = X_norm, Y = Y_scaled)

save.image(file = 'skewBART/R/univariate_data.RData')

