#include "fun.h"
#include "cen_est.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List check_cen_est(std::string dist_name,
                   std::string model_name,
                   List mydata,
                   NumericMatrix y0,
                   NumericMatrix y1,
                   double ll_initial,
                   double approx_weight,
                   double tolerance,
                   int stopCount) {
  List output;
  NumericVector obs = as<NumericVector>(mydata["obs"]);
  int ssize = obs.size();
  NumericVector cen1 = as<NumericVector>(mydata["cen1"]);
  NumericVector cen2 = as<NumericVector>(mydata["cen2"]);
  NumericVector cen3 = as<NumericVector>(mydata["cen3"]);
  NumericVector cen4 = as<NumericVector>(mydata["cen4"]);
  NumericVector gh = as<NumericVector>(mydata["gh"]);
  NumericMatrix p1x = as<NumericMatrix>(mydata["p1x"]);
  NumericMatrix p2x = as<NumericMatrix>(mydata["p2x"]);
  NumericMatrix cxx = as<NumericMatrix>(mydata["cxx"]);
  int p1x_size = p1x.ncol();
  int p2x_size = p2x.ncol();
  int cxx_size = cxx.ncol();
  int all_para_size = p1x_size+p2x_size+cxx_size-2;
  NumericVector b_index(ssize);
  NumericVector cen5(ssize);
  NumericMatrix cen1x(ssize, 1), cen2x(ssize, all_para_size);
  for(int i = 0; i < ssize; i++){
    b_index[i] = 1;
    cen5[i] = 1-cen1[i];
    cen1x(i, 0) = 1;
    for (int j = 0; j < all_para_size; j++) {
      if (j >= 0 && j < p1x_size) {
        cen2x(i, j) = p1x(i, j);
      } else if (j >= p1x_size && j < (p1x_size+p2x_size-1)) {
        int jnd = j-p1x_size+1;
        cen2x(i, j) = p2x(i, jnd);
      } else {
        int jnd = j-p1x_size-p2x_size+2;
        cen2x(i, j) = cxx(i, jnd);
      }
    }
  }
  surv C;
  C.base(dist_name,
         FALSE,
         approx_weight,
         ll_initial,
         obs,
         b_index,
         cen5,
         cen2,
         cen3,
         cen4,
         gh,
         cen1x,
         cen2x,
         cxx);
  List outc = cen_est(model_name, 
                      y0, 
                      y1, 
                      &C, 
                      tolerance, 
                      stopCount);
  NumericVector cest = as<NumericVector>(outc["est"]);
  output["cen_est"] = cest;
  return outc;
}
