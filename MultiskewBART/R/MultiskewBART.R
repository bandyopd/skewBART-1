#' MCMC options for skewBART
#'
#' @param X NxP matrix of training data covariates.
#' @param Y Nxk matrix of training data response.
#' @param test_X MxP matrix of test data covariates
#' @param tree_num the number of trees
#' @param iter Number of MCMC iteration
#' @param burn Number of warmup iterations for the chain.
#'
#' @return Returns a list with the following components:
#' \itemize{
#'   \item y_hat_train: fit to the training data for each iteration of the chain
#'   \item y_hat_test: fit to the testing data for each iteration of the chain
#'   \item Sigma: posterior samples of the covariance matrix
#'   \item lambda: posterior samples of skewness matrix
#' }
#' @export
#'
MultiskewBART <- function(X = X_norm, Y = Y_scaled, test_X = test_X_norm, tree_num=200, iter= 5000, burn =2500){
  hypers$num_tree <- tree_num
  opts   <- Opts()
  my_forest <- MakeForest(hypers = hypers, opts = opts)
  Lam_est <- NULL
  Sigma_out <- array(NA, c(2,2,iter))
  mu_out <- array(NA, c(nrow(Y), 2, iter))
  mu_test_out <- array(NA, c(nrow(test_X), 2, iter))
  Sigma_chain <- hypers$Sigma_hat
  l <- c(0, 0);  u <- c(Inf, Inf)

  Z <- cbind(rep(1,nrow(Y)), rep(1,nrow(Y)))
  Lambda <- c(1,1)
  R <- Y_scaled - Z %*% diag(Lambda)

  for(j in 1:iter){
    R <- Y - Z %*% diag(Lambda)
    my_forest$set_sigma(Sigma_chain)
    mu_hat_chain <- my_forest$do_gibbs(X, R, X, 1)[,,1]
    delta <- Y - mu_hat_chain
    Sigma_chain <- update_sigma(delta, Sigma_chain, hypers)
    mu_hat_test <- my_forest$predict(test_X)
    c(Lambda, mean_Z, Sigma_Z) %<-% update_z_multi_2(Y, mu_hat_chain, Sigma_chain, Z)
    Lambda <- as.numeric(Lambda)
    Lam_est <- rbind(Lam_est, Lambda)
    Z <- t(sapply(1:nrow(Z), function(i) rtmvnorm(1, mu = mean_Z[i,], sigma=Sigma_Z, lb=l, ub=u)))
    Sigma_out[,,j] <- Sigma_chain
    mu_out[,,j] <- mu_hat_chain
    mu_test_out[,,j] <- mu_hat_test
    if(j %% 100 == 0) (cat("\rFinishing iteration", j, "of", iter))
  }


  EST_Lam <- colMeans(Lam_est[c(burn:iter),]) * scale_Y

  EST_mu <- t(t(apply(mu_out[,,c(burn:iter)], c(1,2), mean) %*% diag(scale_Y)) + center_Y)
  EST_train_Y <- t( t(EST_mu) + sqrt(2/pi) * EST_Lam )

  EST_t_mu <- t(t(apply(mu_test_out[,,c(burn:iter)], c(1,2), mean) %*% diag(scale_Y)) + center_Y)
  EST_test_Y <- t( t(EST_t_mu) + sqrt(2/pi) * EST_Lam )

  skew_Sig1 <- apply(Sigma_out[,,c(burn:iter)], c(1,2), mean)
  skew_Sig <- matrix(c(skew_Sig1[1,1]*scale_Y[1]^2, skew_Sig1[1,2]*scale_Y[1]*scale_Y[2],
                       skew_Sig1[2,1]*scale_Y[1]*scale_Y[2],
                       skew_Sig1[2,2]*scale_Y[2]^2),2,2)

  EST_Tau <- c(skew_Sig[1,1], skew_Sig[2,2])
  return(list(y_hat_train = EST_train_Y, y_hat_test = EST_test_Y, Sigma = skew_Sig, lambda = EST_Lam))
}


