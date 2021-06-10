Hypers <- function(X, Y, group = NULL, alpha = 1, beta = 2, gamma = 0.95, k = 2,
                   num_tree = 50, temperature = 1.0) {


  J <- ncol(Y)
  P <- ncol(X)

  if(is.null(group)) group <- Matrix(diag(P), sparse = TRUE)
  if(J != 2) stop("CURRENTLY ONLY SUPPORT TWO DIMENSIONAL RESPONSE")

  Sigma_mu_hat <- diag(2) * (3 / k / sqrt(num_tree))^2

  ## Get Sigma hyperparameters
  fit_lm_Y_1 <- lm(Y[,1] ~ X)
  fit_lm_Y_2 <- lm(Y[,2] ~ X + Y[,1])
  fit_lm_Y_2_x <- lm(Y[,2] ~ X)

  sigma_1_hat <- summary(fit_lm_Y_1)$sigma
  sigma_2_hat <- summary(fit_lm_Y_2)$sigma
  r_hat       <- sigma_2_hat / sigma_1_hat
  Sigma_hat   <- cov(cbind(residuals(fit_lm_Y_1), residuals(fit_lm_Y_2_x)))

   

  out <- list(alpha = alpha, beta = beta, gamma = gamma,
              Sigma_mu_hat = Sigma_mu_hat, k = k, num_tree = num_tree,
              Sigma_hat = Sigma_hat, temperature = temperature,
              sigma_1_hat = sigma_1_hat, sigma_2_hat = sigma_2_hat,
              r_hat = r_hat, group = group)
  
  return(out)

}
