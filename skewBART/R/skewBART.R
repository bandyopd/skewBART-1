#' Create hyperparameter object for skewBART (Hypers function from SoftBART)
#'
#' Borrows the Hypers function from SoftBART packages to incorporate the SoftBART model as a component.
#' Creates a list which holds all the hyperparameters for use with the skewBART.
#'
#'
#' @param X NxP matrix of training data covariates.
#' @param Y Nx1 vector of training data response.
#' @param group For each column of X, gives the associated group
#' @param alpha Positive constant controlling the sparsity level
#' @param beta Parameter penalizing tree depth in the branching process prior
#' @param gamma Parameter penalizing new nodes in the branching process prior
#' @param k Related to the signal-to-noise ratio, sigma_mu = 0.5 / (sqrt(num_tree) * k). BART defaults to k = 2.
#' @param sigma_hat A prior guess at the conditional variance of Y. If not provided, this is estimated empirically by linear regression.
#' @param shape Shape parameter for gating probabilities
#' @param width Bandwidth of gating probabilities
#' @param num_tree Number of trees in the ensemble
#' @param alpha_scale Scale of the prior for alpha; if not provided, defaults to P
#' @param alpha_shape_1 Shape parameter for prior on alpha; if not provided, defaults to 0.5
#' @param alpha_shape_2 Shape parameter for prior on alpha; if not provided, defaults to 1.0
#' @param num_tree_prob Parameter for geometric prior on number of tree
#'
#' @return Returns a list containing the function arguments.
Hypers <- function(X,Y, group = NULL, alpha = 1, beta = 2, gamma = 0.95, k = 2,
                   sigma_hat = NULL, shape = 1, width = 0.1, num_tree = 20,
                   alpha_scale = NULL, alpha_shape_1 = 0.5,
                   alpha_shape_2 = 1, tau_rate = 10,
                   temperature = 1.0) {

  if(is.null(alpha_scale)) alpha_scale <- ncol(X)

  out                                  <- list()

  out$alpha                            <- alpha
  out$beta                             <- beta
  out$gamma                            <- gamma
  out$sigma_mu_hat                     <- 3.5 / (k * sqrt(num_tree))
  out$k                                <- k
  out$num_tree                         <- num_tree
  out$shape                            <- shape
  if(is.null(group)) {
    out$group                          <- 1:ncol(X) - 1
  } else {
    out$group                          <- group - 1
  }

  Y                                    <- scale(Y)
  if(is.null(sigma_hat))
    sigma_hat                          <- GetSigma(X,Y)

  out$sigma_hat                        <- sigma_hat

  out$alpha_scale                      <- alpha_scale
  out$alpha_shape_1                    <- alpha_shape_1
  out$alpha_shape_2                    <- alpha_shape_2
  out$tau_rate                         <- tau_rate
  out$temperature                      <- temperature

  return(out)

}

#' MCMC options for skewBART
#'
#' @param X NxP matrix of training data covariates.
#' @param Y Nx1 vector of training data response.
#' @param test_X MxP matrix of test data covariates
#' @param tree_num the number of trees
#' @param iter Number of MCMC iteration
#' @param burn Number of warmup iterations for the chain.
#'
#' @return Returns a list with the following components:
#' \itemize{
#'   \item y_hat_train: fit to the training data for each iteration of the chain
#'   \item y_hat_test: fit to the testing data for each iteration of the chain
#'   \item sigma: posterior samples of the error standard deviations
#'   \item lambda: posterior samples of lambda
#' }
#'
#' @export

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



