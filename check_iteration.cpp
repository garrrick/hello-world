#include "fun.h"
#include "iteration.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List check_iteration(std::string dist_name,
                     List mydata,
                     NumericMatrix para,
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
  surv cc;
  cc.base(dist_name,
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
  int inn = para.nrow();
  int para_size = p1x.ncol()+p2x.ncol();
  NumericMatrix est(inn, para_size), gra(inn, para_size);
  LogicalVector conv(inn);
  IntegerVector coun(inn);
  NumericVector li(inn), fd_norm(inn), fd_length(inn);
  StringVector c_status;
  List steps;
  for (int i = 0; i < inn; i++) {
    NumericVector parai(para_size);
    NumericVector cui(cxx.ncol());
    for (int j = 0; j < para.ncol(); j++) {
      if (j < para_size) {
        parai[j] = para(i,j);
      } else {
        int jnd = j - para_size;
        cui[jnd] = para(i,j);
      }
    }
    forwarder_input_cure(&cc, cui);
    List outi = iteration(parai, &forwarder_ite_model_llh, &forwarder_ite_model_score, &cc, tolerance, stopCount);
    est(i,_) = as<NumericVector>(outi["est"]);
    gra(i,_) = as<NumericVector>(outi["gradient"]);
    conv[i] = as<bool>(outi["converge"]);
    coun[i] = as<int>(outi["final_count"]);
    li[i] = as<double>(outi["likelihood"]);
    std::string cs = as<std::string>(outi["converge_status"]);
    c_status.push_back(cs);
    NumericVector stepr = as<NumericVector>(outi["step_size"]);
    steps.push_back(stepr);
    fd_norm[i] = as<double>(outi["feasible_direction_norm"]);
    fd_length[i] = as<double>(outi["feasible_direction_length"]);
  }
  List out;
  out["est"] = est;
  out["gradient"] = gra;
  out["conv"] = conv;
  out["coun"] = coun;
  out["li"] = li;
  out["converge_status"] = c_status;
  out["feasible_direction_norm"] = fd_norm;
  out["feasible_direction_length"] = fd_length;
  out["step"] = steps;
  return out;
}
