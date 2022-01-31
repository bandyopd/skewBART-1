#' Fit the MultiskewBART model
#'
#' Fits the MultiskewBART model of Um et al. (2021+). The model is of the form \deqn{Y_i = \mu(X_i) + \epsilon_i} where \eqn{\epsilon_i} has a multivariate skew-normal distribution.
#'
#' @param X NxP matrix of training data covariates.
#' @param Y Nxk matrix of training data response.
#' @param test_X MxP matrix of test data covariates.
#' @param hypers a list containing the hyperparameters of the model, usually constructed using the function Hypers().
#' @param opts a list containing options for running the chain, usually constructed using the function Opts().
#' @param do_skew logical, if true fit the skew-normal model; otherwise, fit a multivariate normal.
#'
#' @return Returns a list with the following components:
#' \itemize{
#'   \item mu_train: fit of the regression function to the training data for
#'         each iteration of the chain; note that the errors are _not_ mean 0,
#'         so this does not give the expected value.
#'   \item mu_test: fit of the regression function to the testing data for each
#'         iteration of the chain; note that the errors are _not_ mean 0, so
#'         this does not give the expected value.
#'   \item Sigma: posterior samples of the covariance matrix
#'   \item lambda: posterior samples of skewness parameters
#' }
#' @export
#'
MultiskewBART <- function(X, Y, test_X, hypers = NULL, opts = NULL, do_skew = TRUE, Wishart = FALSE){

  mean_Z <- Sigma_Z <- NULL

  if(is.null(hypers)) hypers <- Hypers(X,Y)
  if(is.null(opts)) opts <- Opts()
  iter <- opts$num_burn + opts$num_save
  burn <- opts$num_burn

  my_forest <- MakeForest(hypers = hypers, opts = opts)
  Lam_est <- NULL
  Sigma_out <- array(NA, c(2,2,iter))
  mu_out <- array(NA, c(nrow(Y), 2, iter))
  mu_test_out <- array(NA, c(nrow(test_X), 2, iter))
  Sigma_chain <- hypers$Sigma_hat
  l <- c(0, 0);  u <- c(Inf, Inf)
  if(Wishart) {
    S0 <- hypers$S0
    nu <- hypers$nu
   }
    
  XXtest <- quantile_normalize(X, test_X)
  X <- XXtest$X
  test_X <- XXtest$test_X

  Y_scaled <- scale(Y)
  center_Y <- attributes(Y_scaled)$`scaled:center`
  scale_Y <- attributes(Y_scaled)$`scaled:scale`
  Y <- Y_scaled

  Z <- cbind(rep(1,nrow(Y)), rep(1,nrow(Y))) * do_skew
  Lambda <- c(1,1)
  R <- Y_scaled - Z %*% diag(Lambda)

  for(j in 1:iter){
    R <- Y - Z %*% diag(Lambda)
    my_forest$set_sigma(Sigma_chain)
    mu_hat_chain <- my_forest$do_gibbs(X, R, X, 1)[,,1]
    delta <- R - mu_hat_chain
    if(Wishart) {
      B <- delta - Z %*% diag(Lambda)
      Sigma_chain <- solve(rwish( nrow(Y_scaled) + nu , solve(t(B) %*% B  + S0) ))
    }
    Sigma_chain <- update_sigma(delta, Sigma_chain, hypers)
    mu_hat_test <- my_forest$predict(test_X)
    if(do_skew) {
      c(Lambda, mean_Z, Sigma_Z) %<-% update_z_multi_2(Y, mu_hat_chain, Sigma_chain, Z)
      Lambda <- as.numeric(Lambda)
      Lam_est <- rbind(Lam_est, Lambda)
      Z <- t(sapply(1:nrow(Z), function(i) rtmvnorm(1, mu = mean_Z[i,], sigma=Sigma_Z, lb=l, ub=u)))
    }
    Sigma_out[,,j] <- Sigma_chain
    mu_out[,,j] <- mu_hat_chain
    mu_test_out[,,j] <- mu_hat_test
    if(j %% opts$num_print == 0) {
      cat("\rFinishing iteration", j, "of", iter)
    }
  }

  if(do_skew) {
    lambda <- t(t(Lam_est[(burn+1):iter,]) * scale_Y)
  } else {
    lambda <- rep(0, iter - burn)
  }
  mu <- mu_out[,,(burn+1):iter]
  mu[,1,] <- mu[,1,] * scale_Y[1] + center_Y[1]
  mu[,2,] <- mu[,2,] * scale_Y[2] + center_Y[2]
  mu_test <- mu_test_out[,,(burn+1):iter]
  mu_test[,1,] <- mu_test[,1,] * scale_Y[1] + center_Y[1]
  mu_test[,2,] <- mu_test[,2,] * scale_Y[2] + center_Y[2]
  Sigma <- Sigma_out[,,(burn+1):iter]
  for(j in 1:2) {
    for(k in 1:2) {
      Sigma[j,k,] <- Sigma[j,k,] * scale_Y[j] * scale_Y[k]
    }
  }

 
  Y_hat_train_mean <- apply(mu, c(1,2), mean)
  Y_hat_test_mean <- apply(mu_test, c(1,2), mean)
  return(list(y_hat_train = mu, y_hat_test = mu_test, 
              y_hat_train_mean = Y_hat_train_mean, y_hat_test_mean = Y_hat_test_mean,
              lambda = lambda, Sigma = Sigma))
}


