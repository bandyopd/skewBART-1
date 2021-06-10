#ifndef FOREST_H
#define FOREST_H

struct Forest;

#include <RcppArmadillo.h>
#include "node.h"
#include "hypers.h"
#include "opts.h"
#include "mcmc.h"


class Forest {

 private:

  Nodevec trees;
  Hypers hypers;
  Opts opts;

 public:

  Forest(Rcpp::List hypers, Rcpp::List opts);
  ~Forest();
  void IterateGibbs(arma::mat& Y_hat, const arma::mat& X, const arma::mat& Y);
  arma::cube do_gibbs(const arma::mat& X,
                     const arma::mat& Y,
                     const arma::mat& X_test,
                     int num_iter);

  // Getters and Setters, from R
  void set_s(const arma::vec& s);
  void set_sigma(const arma::mat& Sigma);
  arma::vec get_s();
  Rcpp::List get_params();

  // Prediction interfact
  arma::mat predict_mat(const arma::mat& X_test);

  // Var counts
  arma::uvec get_counts();

};








#endif
