skewBART <- function(X = X_norm, Y = Y_scaled, test_X = test_X_norm, tree_num=20, iter= 500, burn =250){
  n_train <- nrow(X)
  mu_hat <- matrix(NA, iter, n_train)
  mu_hat_test <- matrix(NA, iter, nrow(test_X))
  hypers$num_tree <- tree_num
  opts <- Opts()
  forest <- MakeForest(hypers, opts)
  EST_Tau <- EST_Lambda <-  NULL
  Z <- rep(1,n_train);  Lambda <- 1
  #like_skew <- matrix(NA, n_train, iter)

  for(i in 1:iter){
    R <- Y - Lambda*Z
    mu_hat[i,] <- forest$do_gibbs(X, R, X, 1)
    mu_hat_test[i,] <- forest$predict(test_X)
    Tau <- forest$get_params()[["sigma"]]^2
    delta <- Y - mu_hat[i,]
    Z <- rtruncnorm(1, a=0, b=Inf, mean = delta*Lambda/(Lambda^2 + Tau), sd = sqrt(Tau/(Lambda^2+Tau)))
    Lambda <- rnorm(1, mean = (t(Z) %*% delta)/(t(Z)%*%Z), sd = 1/sqrt(t(Z)%*%Z/Tau ))
    EST_Tau <- c(EST_Tau, Tau)
    EST_Lambda <- c(EST_Lambda, Lambda)
    #like_skew[,i] <- dsn(Y_scaled - mu_hat[i,], tau=0, omega= Omega, alpha = Alpha , dp=NULL, log=TRUE)
    if(i %% 100 == 0) (cat("\rFinishing iteration", i, "of",iter,"\t"))
  }

  Lam_est <- mean(EST_Lambda[-c(1:burn)]) * scale_Y
  Tau_est <- mean(EST_Tau[-c(1:burn)])  * scale_Y^2

  Y_hat_train1 <- colMeans(mu_hat[-c(1:burn),]) * scale_Y + center_Y
  Y_hat_test1 <- colMeans(mu_hat_test[-c(1:burn),]) * scale_Y + center_Y

  Y_hat_train <- Y_hat_train1 + Lam_est * sqrt(2/pi)
  Y_hat_test <- Y_hat_test1 + Lam_est * sqrt(2/pi)
  return(list(y_hat_train = Y_hat_train, y_hat_test = Y_hat_test, sigma = Tau_est, lambda = Lam_est))
}



