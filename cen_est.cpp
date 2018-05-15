#include "fun.h"
#include "iteration.h"
#include "cen_est.h"
#include <Rcpp.h>

using namespace Rcpp;

List cen_est(std::string mtype,
             NumericMatrix mi1,
             NumericMatrix mi2,
             void* class_address,
             double tole,
             int stop_count) {
  List output;
  int inn = mi1.nrow();
  int mp1_size = mi1.ncol();
  int mp2_size = mi2.ncol();
  int model_size = mp1_size+mp2_size;
  NumericMatrix est(inn, model_size);
  NumericVector l(inn);
  LogicalVector converge(inn);
  IntegerVector final_count(inn);
  LogicalVector cand(inn);
  List phase;
  for (int i = 0; i < inn; ++i) {
    NumericVector ini_mp1 = mi1(i,_);
    NumericVector ini_mp2 = mi2(i,_);
    NumericVector model_initial(model_size);
    for (int j = 0; j < model_size; j++) {
      if (j < mp1_size) {
        model_initial[j] = ini_mp1[j];
      } else {
        int jnd = j-mp1_size;
        model_initial[j] = ini_mp2[jnd];
      }
    }
    List model_out = iteration(model_initial, &forwarder_model_llh, &forwarder_model_score, class_address, tole, stop_count);
    NumericVector model_est = as<NumericVector>(model_out["est"]);
    bool model_converge = as<bool>(model_out["converge"]);
    int model_count = as<int>(model_out["final_count"]);
    est(i,_) = model_est;
    converge[i] = model_converge;
    final_count[i] = model_count;
    l[i] = forwarder_model_llh(class_address, model_est);

    NumericVector estl = model_est;
    estl.push_back(l[i]);
    LogicalVector check_na_estl = !is_na(estl);
    LogicalVector check_nan_estl = !is_nan(estl);
    LogicalVector check_inf_estl = !is_infinite(estl);
    LogicalVector check_estl = check_na_estl & check_nan_estl & check_inf_estl;
    bool check_all_estl = is_true(all(check_estl));
    cand[i] = check_all_estl & model_converge;
  }
  output["phase"] = phase;
  DataFrame summary = DataFrame::create(Named("est") = est,
                                        Named("converge") = converge,
                                        Named("count") = final_count,
                                        Named("likelihood") = l,
                                        Named("candidate") = cand);
  output["summary"] = summary;
  NumericVector candl, out(model_size);
  double maxl;
  if (is_true(any(cand))) {
    for (int i = 0; i < inn; i++) {
      if (cand[i]) {
        candl.push_back(l[i]);
      }
    }
    maxl = max(candl);
    for (int i = 0; i < inn; i++) {
      if (l[i] == maxl) {
        out = est(i,_);
      }
    }
  } else {
    for (int i = 0; i < model_size; i++) {
      if (i < mp1_size) {
        out[i] = mi1(0,i);
      } else if ((i >= mp1_size) && (i < model_size)) {
        int jnd = i-mp1_size;
        out[i] = mi2(0,jnd);
      }
    }
  }
  output["est"] = out;
  return output;
}
