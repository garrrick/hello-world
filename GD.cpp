#include "fun.h"
#include "iteration.h"
#include "GD.h"
#include <Rcpp.h>

using namespace Rcpp;

List GD(std::string mtype,
        NumericMatrix y0,
        NumericMatrix y1,
        NumericMatrix y2,
        void* class_address,
        double tole,
        int stop_count) {
  List output;
  int inn = y0.nrow();
  int mp1_size = y0.ncol();
  int mp2_size = y1.ncol();
  int cure_size = y2.ncol();
  int model_size = mp1_size+mp2_size;
  int allpara_size = model_size+cure_size;
  NumericMatrix est(inn, allpara_size);
  LogicalVector converge(inn);
  IntegerVector count(inn);
  NumericVector l(inn), li(inn);
  LogicalVector cand1(inn), cand2(inn);
  for (int i = 0; i < inn; ++i) {
    NumericVector ini_mp1 = y0(i,_);
    NumericVector ini_mp2 = y1(i,_);
    NumericVector ini_cu = y2(i,_);
    NumericVector all_initial(allpara_size);
    for (int j = 0; j < allpara_size; j++) {
      if (j < mp1_size) {
        all_initial[j] = ini_mp1[j];
      } else if ((j >= mp1_size)&(j < model_size)) {
        int jnd = j-mp1_size;
        all_initial[j] = ini_mp2[jnd];
      } else {
        int jnd = j-model_size;
        all_initial[j] = ini_cu[jnd];
      }
    }
    List all_out;
    if (mtype.compare("mixture") == 0) {
      all_out = iteration(all_initial, &forwarder_mixture_llh, &forwarder_mixture_score, class_address, tole, stop_count);
    } else if (mtype.compare("non-mixture") == 0) {
      all_out = iteration(all_initial, &forwarder_nmixture_llh, &forwarder_nmixture_score, class_address, tole, stop_count);
    } else if (mtype.compare("time") == 0) {
      all_out = iteration(all_initial, &forwarder_llh, &forwarder_score, class_address, tole, stop_count);
    } else {
      stop("Choose mixture, non-mixture, or time as model type");
    }
    NumericVector all_est = as<NumericVector>(all_out["est"]);
    bool all_converge = as<bool>(all_out["converge"]);
    int all_count = as<int>(all_out["final_count"]);
    double sl = as<double>(all_out["likelihood"]);
    est(i,_) = all_est;
    converge[i] = all_converge;
    count[i] = all_count;
    l[i] = sl;

    NumericVector estl = all_est;
    if (mtype.compare("time") == 0) {
      li[i] = forwarder_llh_indicator(class_address, all_est);
    } else {
      li[i] = sl;
    }
    estl.push_back(li[i]);
    LogicalVector check_na_estl = !is_na(estl);
    LogicalVector check_nan_estl = !is_nan(estl);
    LogicalVector check_inf_estl = !is_infinite(estl);
    LogicalVector check_estl = check_na_estl & check_nan_estl & check_inf_estl;
    bool check_all_estl = is_true(all(check_estl));
    cand1[i] = check_all_estl & all_converge;
    cand2[i] = check_all_estl & !all_converge;
  }
  DataFrame summary = DataFrame::create(Named("est") = est,
                                        Named("converge") = converge,
                                        Named("count") = count,
                                        Named("smoothed_likelihood") = l,
                                        Named("likelihood") = li,
                                        Named("candidate1") = cand1,
                                        Named("candidate2") = cand2);
  output["summary"] = summary;
  bool cand1conv = is_true(any(cand1));
  output["converge"] = cand1conv;
  NumericVector candl, out(allpara_size);
  double maxl;
  if (is_true(any(cand1))) {
    for (int i = 0; i < inn; i++) {
      if (cand1[i]) {
        candl.push_back(li[i]);
      }
    }
    maxl = max(candl);
    for (int i = 0; i < inn; i++) {
      if (li[i] == maxl) {
        out = est(i,_);
      }
    }
  } else if (is_true(any(cand2))) {
    for (int i = 0; i < inn; i++) {
      if (cand2[i]) {
        candl.push_back(li[i]);
      }
    }
    maxl = max(candl);
    for (int i = 0; i < inn; i++) {
      if (li[i] == maxl) {
        out = est(i,_);
      }
    }
  } else {
    for (int i = 0; i < allpara_size; i++) {
      if (i < mp1_size) {
        out[i] = y0(0,i);
      } else if ((i >= mp1_size)&(i < model_size)) {
        int jnd = i-mp1_size;
        out[i] = y1(0,jnd);
      } else {
        int jnd = i-model_size;
        out[i] = y2(0,jnd);
      }
    }
  }
  output["est"] = out;
  if (mtype.compare("mixture") == 0) {
    output["loglikelihood"] = forwarder_mixture_llh(class_address, out);
    output["score"] = forwarder_mixture_score(class_address, out);
  } else if (mtype.compare("non-mixture") == 0) {
    output["loglikelihood"] = forwarder_nmixture_llh(class_address, out);
    output["score"] = forwarder_nmixture_score(class_address, out);
  } else if (mtype.compare("time") == 0) {
    output["loglikelihood"] = forwarder_llh_indicator(class_address, out);
    output["score"] = forwarder_score(class_address, out);
  }
  return output;
}
