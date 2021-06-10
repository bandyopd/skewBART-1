#include "hypers.h"

using namespace Rcpp;
using namespace arma;

// void Hypers::UpdateSigma(const arma::vec& residuals) {
//   sigma = update_sigma_halfcauchy(residuals, sigma_hat, sigma, temperature);
// }

// void Hypers::UpdateSigmaMu(const arma::vec& means) {
//   sigma_mu = update_sigma_halfcauchy(means, sigma_mu_hat, sigma_mu, temperature);
// }


Hypers::Hypers(Rcpp::List hypers)
{
  alpha        = hypers["alpha"];
  beta         = hypers["beta"];
  gamma        = hypers["gamma"];
  Sigma_hat    = as<mat>(hypers["Sigma_hat"]);
  Sigma        = Sigma_hat;
  Sigma_mu_hat = as<mat>(hypers["Sigma_mu_hat"]);
  Sigma_mu     = Sigma_mu_hat;
  A            = inv(Sigma);
  B            = inv(Sigma_mu);
  temperature  = hypers["temperature"];
  num_tree    = hypers["num_tree"];

  sp_mat tmp = as<sp_mat>(hypers["group"]);
  group = zeros<sp_umat>(tmp.n_rows, tmp.n_cols);
  for(sp_mat::iterator it = tmp.begin(); it != tmp.end(); ++it) {
    group(it.row(), it.col()) = std::round(*it); 
  }
  group_size = sum(group, 1);
  p_group_var = tmp;
  for(sp_umat::iterator it = group.begin(); it != group.end(); ++it) {
    p_group_var(it.row(), it.col()) = 1.0 / group_size(it.row());
  }


  // Initialize s
  s = 1.0 / group.n_rows * ones<vec>(group.n_rows);
  logs = log(s);
}

int Hypers::SampleVar() const {
  uvec group_and_var = zeros<uvec>(2);
  group_and_var(0) = sample_class(s);
  group_and_var(1) = sample_class_row(p_group_var, group_and_var(0));
  return group_and_var(1);
}





/*Note: Because the shape of the Dirichlet will mostly be small, we sample from
  the Dirichlet distribution by sampling log-gamma random variables using the
  technique of Liu, Martin, and Syring (2017+) and normalizing using the
  log-sum-exp trick */
// void UpdateS(std::vector<Node*>& forest, Hypers& hypers) {

//   // Get shape vector
//   vec shape_up = hypers.alpha / ((double)hypers.s.size()) * ones<vec>(hypers.s.size());
//   shape_up = shape_up + get_var_counts(forest, hypers);

//   // Sample unnormalized s on the log scale
//   for(int i = 0; i < shape_up.size(); i++) {
//     hypers.logZ(i) = rlgam(shape_up(i));
//   }
//   // Normalize s on the log scale, then exponentiate
//   hypers.logs = hypers.logZ - log_sum_exp(hypers.logZ);
//   hypers.s = exp(hypers.logs);
//   hypers.logZ = hypers.logs + rlgam(hypers.alpha);

// }

// // NOTE: the log-likelihood here is -n Gam(alpha/n) + alpha * mean_log_Z + (shape - 1) * log(alpha) - rate * alpha
// void Hypers::UpdateAlpha() {


//   // Get the Gamma approximation

//   double n = logZ.size();
//   double R = mean(logZ); mean_log_Z = R;
//   double alpha_hat = exp(log_sum_exp(logZ));
//   a_hat = alpha_shape_1 + alpha_hat * alpha_hat * Rf_trigamma(alpha_hat / n) / n;
//   b_hat = 1.0 / alpha_scale + (a_hat - alpha_shape_1) / alpha_hat +
//     Rf_digamma(alpha_hat / n) - R;
//   int M = 10;
//   for(int i = 0; i < M; i++) {
//     alpha_hat = a_hat / b_hat;
//     a_hat = alpha_shape_1 + alpha_hat * alpha_hat * Rf_trigamma(alpha_hat / n) / n;
//     b_hat = 1.0 / alpha_scale + (a_hat - alpha_shape_1) / alpha_hat +
//       Rf_digamma(alpha_hat / n) - R;
//   }
//   double A = a_hat * .75;
//   double B = b_hat * .75;

//   // double n = logZ.size();
//   // double R = sum(logZ);
//   // double alpha_hat = exp(log_sum_exp(logZ)) / n;
//   // a_hat = 1.0 + alpha_hat * alpha_hat * n * Rf_trigamma(alpha_hat);
//   // b_hat = (a_hat - 1.0) / alpha_hat + n * Rf_digamma(alpha_hat) - R;
//   // int M = 10;
//   // for(int i = 0; i < M; i++) {
//   //   alpha_hat = a_hat / b_hat;
//   //   a_hat = 1.0 + alpha_hat * alpha_hat * n * Rf_trigamma(alpha_hat);
//   //   b_hat = (a_hat - 1.0) / alpha_hat + n * Rf_digamma(alpha_hat) - R;
//   // }
//   // a_hat = a_hat / 1.3;
//   // b_hat = b_hat / 1.3;

//   // Sample from the gamma approximation
//   double alpha_prop = R::rgamma(A, 1.0 / B);


//   // Compute logliks
//   double loglik_new = - n * R::lgammafn(alpha_prop / n) + alpha_prop * R +
//     (alpha_shape_1 - 1.0) * log(alpha_prop) - alpha_prop / alpha_scale +
//     R::dgamma(alpha, A, 1.0 / B, 1);
//   double loglik_old = -n * R::lgammafn(alpha / n) + alpha * R +
//     (alpha_shape_1 - 1.0) * log(alpha) - alpha / alpha_scale +
//     R::dgamma(alpha_prop, A, 1.0 / B, 1);

//   // Accept or reject
//   if(log(unif_rand()) < loglik_new - loglik_old) {
//     alpha = alpha_prop;
//   }

//   // arma::vec logliks = zeros<vec>(rho_propose.size());
//   // rho_loglik loglik;
//   // loglik.mean_log_s = mean(logs);
//   // loglik.p = (double)s.size();
//   // loglik.alpha_scale = alpha_scale;
//   // loglik.alpha_shape_1 = alpha_shape_1;
//   // loglik.alpha_shape_2 = alpha_shape_2;

//   // for(int i = 0; i < rho_propose.size(); i++) {
//   //   logliks(i) = loglik(rho_propose(i));
//   // }

//   // logliks = exp(logliks - log_sum_exp(logliks));
//   // double rho_up = rho_propose(sample_class(logliks));
//   // alpha = rho_to_alpha(rho_up, alpha_scale);

// }

// // void Hypers::UpdateAlpha() {

// //   double rho = alpha_to_rho(alpha, alpha_scale);
// //   double psi = mean(log(s));
// //   double p = (double)s.size();

// //   double loglik = alpha * psi + Rf_lgammafn(alpha) - p * Rf_lgammafn(alpha / p) +
// //     logpdf_beta(rho, alpha_shape_1, alpha_shape_2);

// //   // 50 MH proposals
// //   for(int i = 0; i < 50; i++) {
// //     double rho_propose = Rf_rbeta(alpha_shape_1, alpha_shape_2);
// //     double alpha_propose = rho_to_alpha(rho_propose, alpha_scale);

// //     double loglik_propose = alpha_propose * psi + Rf_lgammafn(alpha_propose) -
// //       p * Rf_lgammafn(alpha_propose/p) +
// //       logpdf_beta(rho_propose, alpha_shape_1, alpha_shape_2);

// //     if(log(unif_rand()) < loglik_propose - loglik) {
// //       alpha = alpha_propose;
// //       rho = rho_propose;
// //       loglik = loglik_propose;
// //     }
// //   }
// // }

// // void Hypers::UpdateAlpha() {

// //   rho_loglik loglik;
// //   loglik.mean_log_s = mean(logs);
// //   loglik.p = (double)s.size();
// //   loglik.alpha_scale = alpha_scale;
// //   loglik.alpha_shape_1 = alpha_shape_1;
// //   loglik.alpha_shape_2 = alpha_shape_2;

// //   double rho = alpha_to_rho(alpha, alpha_scale);
// //   rho = slice_sampler(rho, loglik, 0.1, 0.0 + exp(-10.0), 1.0);
// //   alpha = rho_to_alpha(rho, alpha_scale);
// // }


// // double loglik_data(const arma::vec& Y, const arma::vec& Y_hat, const Hypers& hypers) {
// //   vec res = Y - Y_hat;
// //   double out = -0.5 * Y.size() * log(M_2_PI * pow(hypers.sigma,2.0)) -
// //     dot(res, res) * 0.5 / pow(hypers.sigma,2.0);
// //   return out;
// // }

// arma::vec loglik_data(const arma::vec& Y, const arma::vec& Y_hat, const Hypers& hypers) {
//   vec res = Y - Y_hat;
//   vec out = zeros<vec>(Y.size());
//   for(int i = 0; i < Y.size(); i++) {
//     out(i) = -0.5 * log(M_2_PI * pow(hypers.sigma,2)) - 0.5 * pow(res(i) / hypers.sigma, 2);
//   }
//   return out;
// }
