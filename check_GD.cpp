#include "fun.h"
#include "iteration.h"
#include "GD.h"
#include "ite_GD.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List check_GD(std::string dist_name,
              std::string model_name,
              List mydata,
              NumericMatrix par1,
              NumericMatrix par2,
              NumericMatrix cu,
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
         FALSE,
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
  List output = GD(model_name, par1, par2, cu, &W, tolerance, stopCount);
  return output;
}
