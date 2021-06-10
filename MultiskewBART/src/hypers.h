#ifndef HYPERS_H
#define HYPERS_H

#include <RcppArmadillo.h>
#include "math.h"

class Hypers {

  public:

  double alpha;
  double beta;
  double gamma;
  arma::mat Sigma;
  arma::mat Sigma_hat;
  arma::mat Sigma_mu;
  arma::mat Sigma_mu_hat;
  arma::mat A;
  arma::mat B;

  double temperature;
  int num_tree;
  arma::vec s;
  arma::vec logs;
  arma::sp_umat group;
  arma::uvec group_size;
  arma::sp_mat p_group_var;

  // Updates
  // void UpdateSigma(const arma::vec& residuals);
  // void UpdateSigmaMu(const arma::vec& means);
  /* void UpdateAlpha(); */

  // Construct/Destructure
  Hypers(Rcpp::List hypers);

  // Utilities
  int SampleVar() const;

};












/* arma::vec loglik_data(const arma::vec& Y, const arma::vec& Y_hat, const Hypers& hypers); */
/* double alpha_to_rho(double alpha, double scale); */
/* double rho_to_alpha(double rho, double scale); */
/* double forest_loglik(std::vector<Node*>& forest, double gamma, double beta); */
/* double tree_loglik(Node* node, int node_depth, double gamma, double beta); */
/* void UpdateS(std::vector<Node*>& forest, Hypers& hypers); */


#endif
