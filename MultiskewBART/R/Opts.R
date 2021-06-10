Opts <- function(num_burn = 5000, num_thin = 1, num_save = 5000,
                 num_print = 100, update_Sigma_mu = FALSE, update_Sigma = FALSE,
                 update_s = FALSE, update_alpha = FALSE) {


  return(list(num_burn = num_burn, num_thin = num_thin, num_save = num_save,
              num_print = num_print, update_Sigma_mu = update_Sigma_mu,
              update_Sigma = update_Sigma, update_s = update_s,
              update_alpha = update_alpha))

}
