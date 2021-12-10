#' MCMC options for SoftBart
#'
#' Creates a list which provides the parameters for running the Markov chain.
#'
#' @param num_burn Number of warmup iterations for the chain.
#' @param num_thin Thinning interval for the chain.
#' @param num_save The number of samples to collect; in total, num_burn + num_save * num_thin iterations are run
#' @param num_print Interval for how often to print the chain's progress
#' @param update_sigma_mu If true, sigma_mu/k are updated, with a half-Cauchy prior on sigma_mu centered at the initial guess
#' @param update_s If true, s is updated using the Dirichlet prior.
#' @param update_alpha If true, alpha is updated using a scaled beta prime prior
#' @param update_beta If true, beta is updated using a Normal(0,2^2) prior
#' @param update_gamma If true, gamma is updated using a Uniform(0.5, 1) prior
#' @param update_tau If true, tau is updated for each tree
#' @param update_tau_mean If true, the mean of tau is updated
#'
#' @return Returns a list containing the function arguments
Opts <- function(num_burn = 2500, num_thin = 1, num_save = 2500, num_print = 100,
                 update_sigma_mu = TRUE, update_s = TRUE, update_alpha = TRUE,
                 update_beta = FALSE, update_gamma = FALSE, update_tau = TRUE,
                 update_tau_mean = FALSE, update_sigma = TRUE) {
  out <- list()
  out$num_burn        <- num_burn
  out$num_thin        <- num_thin
  out$num_save        <- num_save
  out$num_print       <- num_print
  out$update_sigma_mu <- update_sigma_mu
  out$update_s         <- update_s
  out$update_alpha    <- update_alpha
  out$update_beta     <- update_beta
  out$update_gamma    <- update_gamma
  out$update_tau      <- update_tau
  out$update_tau_mean <- update_tau_mean
  out$update_sigma    <- update_sigma
  out$update_num_tree <- FALSE

  return(out)

}

#' Fits the SoftBart model
#'
#' Runs the Markov chain for the softbart model and collects the output
#'
#' @param X NxP matrix of training data covariates
#' @param Y Nx1 vector of training data responses
#' @param X_test NxP matrix of test data covariates
#' @param hypers List of hyperparameter values obtained from Hypers function
#' @param opts List of MCMC chain settings obtained from Opts function
#'
#' @return Returns a list with the following components:
#' \itemize{
#'   \item y_hat_train: fit to the training data for each iteration of the chain
#'   \item y_hat_test: fit to the testing data for each iteration of the chain
#'   \item y_hat_train_mean: fit to the training data, averaged over iterations
#'   \item y_hat_test_mean: fit to the testing data, averaged over iterations
#'   \item sigma: posterior samples of the error standard deviations
#'   \item sigma_mu: posterior samples of sigma_mu, the standard deviation of the leaf node parameters
#'   \item s: posterior samples of s
#'   \item alpha: posterior samples of alpha
#'   \item beta: posterior samples of beta
#'   \item gamma: posterior samples of gamma
#'   \item k: posterior samples of k = 0.5 / (sqrt(num_tree) * sigma_mu)
#' }
softbart <- function(X, Y, X_test, hypers = NULL, opts = Opts()) {

  if(is.null(hypers)){
    hypers <- Hypers(X,Y)
  }

  ## Normalize Y
  Z <- scale(Y)

  ## Quantile normalize X
  n <- nrow(X)
  idx_train <- 1:n
  X_trans <- rbind(X, X_test)

  if(is.data.frame(X_trans)) {
    print("Preprocessing data frame")
    preproc_df <- preprocess_df(X_trans)
    X_trans <- preproc_df$X
    print("Using default grouping; if this is not desired, preprocess data frame manually using preprocess_df before calling.")
    hypers$group
    hypers$group <- preproc_df$group
  }

  X_trans <- quantile_normalize_bart(X_trans)
  X <- X_trans[idx_train,,drop=FALSE]
  X_test <- X_trans[-idx_train,,drop=FALSE]

  fit <- SoftBart(X,Z,X_test,
                  hypers$group,
                  hypers$alpha,
                  hypers$beta,
                  hypers$gamma,
                  hypers$sigma,
                  hypers$shape,
                  hypers$width,
                  hypers$num_tree,
                  hypers$sigma_hat,
                  hypers$k,
                  hypers$alpha_scale,
                  hypers$alpha_shape_1,
                  hypers$alpha_shape_2,
                  hypers$tau_rate,
                  hypers$num_tree_prob,
                  hypers$temperature,
                  opts$num_burn,
                  opts$num_thin,
                  opts$num_save,
                  opts$num_print,
                  opts$update_sigma_mu,
                  opts$update_s,
                  opts$update_alpha,
                  opts$update_beta,
                  opts$update_gamma,
                  opts$update_tau,
                  opts$update_tau_mean,
                  opts$update_num_tree)


  attrs <- attributes(Z)
  center_Y <- attrs[[2]]
  scale_Y <- attrs[[3]]

  fit$y_hat_train <- fit$y_hat_train * scale_Y + center_Y
  fit$y_hat_test <- fit$y_hat_test * scale_Y + center_Y
  fit$sigma <- scale_Y * fit$sigma
  fit$k <- 0.5 / (sqrt(hypers$num_tree) * fit$sigma_mu)

  fit$y_hat_train_mean <- colMeans(fit$y_hat_train)
  fit$y_hat_test_mean <- colMeans(fit$y_hat_test)

  fit$y <- Y

  class(fit) <- "softbart"

  return(fit)

}

GetSigma <- function(X,Y) {

  stopifnot(is.matrix(X) | is.data.frame(X))

  if(is.data.frame(X)) {
    X <- model.matrix(~.-1, data = X)
  }


  fit <- cv.glmnet(x = X, y = Y)
  fitted <- predict(fit, X)
  sigma_hat <- sqrt(mean((fitted - Y)^2))

  return(sigma_hat)

}

MakeForest <- function(hypers, opts) {
  mf <- Module(module = "mod_forest", PACKAGE = "skewBART")
  return(new(mf$Forest, hypers, opts))
}
