#include "fun.h"
#include "iteration_curetime.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List check_iteration_curetime(std::string dist_name,
                              List mydata,
                              NumericMatrix mp,
                              arma::mat para,
                              arma::mat excess_covariate,
                              arma::vec excess_event,
                              double ll_initial,
                              double approx_weight,
                              double tolerance,
                              int stopCount,
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
  int inn = para.n_rows;
  int sample_size = obs.size();
  mat est(inn, para.n_cols), gra(inn, para.n_cols);
  LogicalVector conv(inn);
  IntegerVector coun(inn), ss(inn);
  NumericVector li(inn), fd_norm(inn), fd_length(inn);
  StringVector conv_status;
  List steps, Armijo_count, phase, ac;
  for (int i = 0; i < inn; i++) {
    surv cc;
    cc.base(dist_name,
            TRUE,
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
    forwarder_input_model(&cc, mp(i,_));
    rowvec parai = para.row(i);
    List outi = iteration_curetime(parai.t(), &forwarder_ite_curetime_llh, &forwarder_ite_curetime_score, excess_covariate, excess_event, &cc, tolerance, stopCount);
    est(i, span::all) = as<rowvec>(outi["est"]);
    gra(i, span::all) = as<rowvec>(outi["gradient"]);
    conv[i] = as<bool>(outi["converge"]);
    coun[i] = as<int>(outi["final_count"]);
    li[i] = as<double>(outi["likelihood"]);
    std::string cons = as<std::string>(outi["converge_status"]);
    conv_status.push_back(cons);
    NumericVector stepr = as<NumericVector>(outi["step_size"]);
    steps.push_back(stepr);
    NumericVector Ac = as<NumericVector>(outi["Armijo_count"]);
    Armijo_count.push_back(Ac);
    ss[i] = as<int>(outi["smallest_step_index"]);
    fd_norm[i] = as<double>(outi["feasible_direction_norm"]);
    fd_length[i] = as<double>(outi["feasible_direction_length"]);
    StringVector pr = as<StringVector>(outi["phase"]);
    phase.push_back(pr);
    NumericMatrix acc = as<NumericMatrix>(outi["active_set"]);
    ac.push_back(acc);
  }
  List out;
  out["est"] = est;
  out["gradient"] = gra;
  out["conv"] = conv;
  out["coun"] = coun;
  out["conv_status"] = conv_status;
  out["li"] = li;
  out["feasible_direction_norm"] = fd_norm;
  out["feasible_direction_length"] = fd_length;
  out["step"] = steps;
  out["phase"] = phase;
  out["active_set"] = ac;
  out["Armijo_count"] = Armijo_count;
  return out;
}
