#include "forest.h"

using namespace arma;
using namespace Rcpp;

Forest::Forest(Rcpp::List hypers, Rcpp::List opts) : hypers(hypers), opts(opts) {
  trees.resize(this->hypers.num_tree);
  for(int i = 0; i < this->hypers.num_tree; i++) {
    trees[i] = new Node(this->hypers);
  }
}


Forest::~Forest() {
  for(int i = 0; i < trees.size(); i++) delete trees[i];
}

void Forest::IterateGibbs(arma::mat& Y_hat,
                          const arma::mat& X,
                          const arma::mat& Y) {


  TreeBackfit(trees, Y_hat, hypers, X, Y, opts);
  arma::mat res = Y - Y_hat;
  // arma::mat means = get_means(trees);
  // if(opts.update_sigma) hypers.UpdateSigma(res);
  // if(opts.update_sigma_mu) hypers.UpdateSigmaMu(means);
  // if(opts.update_s) UpdateS(trees, hypers);
  // if(opts.update_alpha) hypers.UpdateAlpha();
  // if(opts.update_num_tree) update_num_tree(forest, hypers, opts, Y, Y - Y_hat, X);

  Rcpp::checkUserInterrupt();

}

arma::cube Forest::do_gibbs(const arma::mat& X,
                            const arma::mat& Y,
                            const arma::mat& X_test,
                            int num_iter) {

  mat Y_hat = predict(trees, X, hypers);
  cube Y_out = zeros<cube>(X_test.n_rows, Y.n_cols, num_iter);

  for(int i = 0; i < num_iter; i++) {
    // Rcout << "Iterate";
    IterateGibbs(Y_hat, X, Y);
    // Rcout << "Predict";
    Y_out.slice(i) = predict(trees, X_test, hypers);
    if((i+1) % opts.num_print == 0) {
      Rcout << "\rFinishing iteration " << i+1 << "\t\t\t";
    }
  }

  return Y_out;

}

void Forest::set_s(const arma::vec& s_) {
  hypers.s = s_;
  hypers.logs = log(s_);
}

arma::vec Forest::get_s() {
  return hypers.s;
}

void Forest::set_sigma(const arma::mat& Sigma_) {
  hypers.Sigma = Sigma_;
  hypers.A = inv(Sigma_);
}

Rcpp::List Forest::get_params() {

  List out;
  out["alpha"] = hypers.alpha;
  out["sigma"] = hypers.Sigma;
  out["sigma_mu"] = hypers.Sigma_mu;

  return out;
}

arma::mat Forest::predict_mat(const arma::mat& X_test) {
  return predict(trees, X_test, hypers);
}

arma::uvec Forest::get_counts() {
  return get_var_counts(trees, hypers);
}

RCPP_MODULE(mod_forest) {

  class_<Forest>("Forest")

    .constructor<Rcpp::List, Rcpp::List>()
    .method("do_gibbs", &Forest::do_gibbs)
    .method("get_s", &Forest::get_s)
    .method("set_sigma", &Forest::set_sigma)
    // .method("get_params", &Forest::get_params)
    .method("predict", &Forest::predict_mat)
    .method("get_counts", &Forest::get_counts)
    .method("set_s", &Forest::set_s)
    ;

}
