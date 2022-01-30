// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rmvnorm
arma::vec rmvnorm(const arma::vec& mean, const arma::mat& Precision);
RcppExport SEXP _skewBART_rmvnorm(SEXP meanSEXP, SEXP PrecisionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Precision(PrecisionSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm(mean, Precision));
    return rcpp_result_gen;
END_RCPP
}
// rlgam
double rlgam(double shape);
RcppExport SEXP _skewBART_rlgam(SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(rlgam(shape));
    return rcpp_result_gen;
END_RCPP
}
// update_b
arma::vec update_b(const arma::vec& Z, const arma::vec& b_new, const arma::vec& mu_new, const arma::vec& b_old, const arma::vec& mu_old, const arma::uvec& cluster_idx, double sigma);
RcppExport SEXP _skewBART_update_b(SEXP ZSEXP, SEXP b_newSEXP, SEXP mu_newSEXP, SEXP b_oldSEXP, SEXP mu_oldSEXP, SEXP cluster_idxSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b_new(b_newSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_new(mu_newSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b_old(b_oldSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_old(mu_oldSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type cluster_idx(cluster_idxSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(update_b(Z, b_new, mu_new, b_old, mu_old, cluster_idx, sigma));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_mod_forest();

static const R_CallMethodDef CallEntries[] = {
    {"_skewBART_rmvnorm", (DL_FUNC) &_skewBART_rmvnorm, 2},
    {"_skewBART_rlgam", (DL_FUNC) &_skewBART_rlgam, 1},
    {"_skewBART_update_b", (DL_FUNC) &_skewBART_update_b, 7},
    {"_rcpp_module_boot_mod_forest", (DL_FUNC) &_rcpp_module_boot_mod_forest, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_skewBART(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
