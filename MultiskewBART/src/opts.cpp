#include "opts.h"

using namespace arma;
using namespace Rcpp;


Opts::Opts(Rcpp::List opts) {
  num_burn        = opts["num_burn"];
  num_thin        = opts["num_thin"];
  num_save        = opts["num_save"];
  num_print       = opts["num_print"];
  update_Sigma_mu = opts["update_Sigma_mu"];
  update_Sigma    = opts["update_Sigma"];
  update_s        = opts["update_s"];
  update_alpha    = opts["update_alpha"];
}
