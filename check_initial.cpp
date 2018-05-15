#include "fun.h"
#include "initial.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List check_initial(std::string dist_name,
                   std::string model_name,
                   List mydata,
                   NumericMatrix y0,
                   NumericMatrix y1,
                   NumericMatrix y2,
                   double ll_initial,
                   double approx_weight,
                   double tolerance,
                   int stopCount) {
  NumericVector obs = as<NumericVector>(mydata["obs"]);
  NumericVector cen1 = as<NumericVector>(mydata["cen1"]);
  NumericVector cen2 = as<NumericVector>(mydata["cen2"]);
  NumericVector cen3 = as<NumericVector>(mydata["cen3"]);
  NumericVector cen4 = as<NumericVector>(mydata["cen4"]);
  NumericVector gh = as<NumericVector>(mydata["gh"]);
  NumericMatrix p1x = as<NumericMatrix>(mydata["p1x"]);
  NumericMatrix p2x = as<NumericMatrix>(mydata["p2x"]);
  NumericMatrix cxx = as<NumericMatrix>(mydata["cxx"]);
  NumericVector b_index(obs.size());
  for(int i = 0; i < obs.size(); i++){
    b_index[i] = 1;
  }
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
  mat cx = as<mat>(cxx);
  vec ct = as<vec>(obs);
  List output = initial(model_name, 
                        y0, 
                        y1, 
                        y2,
                        cx,
                        ct,
                        cen3, 
                        &W, 
                        tolerance, 
                        stopCount);
  return output;
}
