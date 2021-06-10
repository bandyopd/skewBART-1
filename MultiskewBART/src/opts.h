#ifndef OPTS_H
#define OPTS_H

#include <RcppArmadillo.h>
#include "math.h"

struct Opts;

struct Opts {

  int num_burn;
  int num_thin;
  int num_save;
  int num_print;

  bool update_Sigma_mu;
  bool update_Sigma;
  bool update_s;
  bool update_alpha;

  Opts(Rcpp::List opts);

};


#endif
