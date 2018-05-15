#include "fun.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gra(std::string dist_name,
                  List mydata,
                  NumericVector bp1,
                  NumericVector bp2,
                  NumericVector cu,
                  NumericVector bindex,
                  double ll_initial,
                  double approx_weight) {
  NumericVector obs = as<NumericVector>(mydata["obs"]);
  NumericVector cen1 = as<NumericVector>(mydata["cen1"]);
  NumericVector cen2 = as<NumericVector>(mydata["cen2"]);
  NumericVector cen3 = as<NumericVector>(mydata["cen3"]);
  NumericVector cen4 = as<NumericVector>(mydata["cen4"]);
  NumericVector gh = as<NumericVector>(mydata["gh"]);
  NumericMatrix p1x = as<NumericMatrix>(mydata["p1x"]);
  NumericMatrix p2x = as<NumericMatrix>(mydata["p2x"]);
  NumericMatrix cxx = as<NumericMatrix>(mydata["cxx"]);
  surv tes;
  tes.base(dist_name,
           FALSE,
           approx_weight,
           ll_initial,
           obs,
           bindex,
           cen1,
           cen2,
           cen3,
           cen4,
           gh,
           p1x,
           p2x,
           cxx);
  NumericVector yy(p1x.ncol()+p2x.ncol());
  for (int i = 0; i < p1x.ncol()+p2x.ncol(); i++) {
    if (i < p1x.ncol()) {
      yy[i] = bp1[i];
    } else {
      int jnd = i-p1x.ncol();
      yy[i] = bp2[jnd];
    }
  }
  tes.input_model(yy);
  NumericVector ics = tes.ite_curetime_score(cu);
  return(ics);
}
