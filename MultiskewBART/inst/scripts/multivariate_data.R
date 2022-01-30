library(MultiskewBART)

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

## preprocess data

X_norm <- quantile_normalize(X, test_X)$X
test_X_norm <- quantile_normalize(X, test_X)$test_X
Y_scaled <- scale(Y)
center_Y <- attributes(Y_scaled)$`scaled:center`
scale_Y <- attributes(Y_scaled)$`scaled:scale`

## Create a list of the hyperparameters of the model. 

hypers <- Hypers(X = X_norm, Y = Y_scaled)
save.image(file = 'MultiskewBART/R/multivariate_data.RData')

