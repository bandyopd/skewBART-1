#ifndef MY_MATH_H
#define MY_MATH_H


#include <RcppArmadillo.h>


int sample_class(const arma::vec& probs);
int sample_class_row(const arma::sp_mat& probs, int row);
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& Precision);
double logit(double x);
double expit(double x);
double activation(double x, double c, double tau);
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& Precision);
double rlgam(double shape);
double log_sum_exp(const arma::vec& x);
double cauchy_jacobian(double tau, double sigma_hat);
double logpdf_beta(double x, double a, double b);

#endif
