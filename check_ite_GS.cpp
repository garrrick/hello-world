#include "fun.h"
#include "iteration.h"
#include "ite_GS.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List check_ite_GS(std::string dist_name,
                  std::string model_name,
                  List mydata,
                  NumericMatrix par1,
                  NumericMatrix par2,
                  double cu_from,
                  double cu_to,
                  double ll_initial,
                  double approx_weight,
                  double tolerance,
                  int stopCount,
                  bool b_logical,
                  NumericVector b_index) {
  NumericVector obs = as<NumericVector>(mydata["obs"]);
  NumericVector cen1 = as<NumericVector>(mydata["cen1"]);
  NumericVector cen2 = as<NumericVector>(mydata["cen2"]);
  NumericVector cen3 = as<NumericVector>(mydata["cen3"]);
  NumericVector cen4 = as<NumericVector>(mydata["cen4"]);
  NumericVector gh = as<NumericVector>(mydata["gh"]);
  NumericMatrix p1x = as<NumericMatrix>(mydata["p1x"]);
  NumericMatrix p2x = as<NumericMatrix>(mydata["p2x"]);
  NumericMatrix cxx = as<NumericMatrix>(mydata["cxx"]);
  surv W;
  W.base(dist_name,
         b_logical,
         approx_weight,
         ll_initial,
         obs,
         b_index,
         cen1,
         cen2,
         cen3,
         cen4,
         gh,
         p1x,
         p2x,
         cxx);
  List output = ite_GS(model_name, par1, par2, cu_from, cu_to, &W, tolerance, stopCount);
  return output;
}
