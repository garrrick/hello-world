#include "fun.h"
#include "iteration.h"
#include "iteration_curetime.h"
#include "initial.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

List initial(std::string mtype,
             NumericMatrix mi1,
             NumericMatrix mi2,
             NumericMatrix ci,
             arma::mat X,
             arma::vec T,
             NumericVector ce,
             void* class_address,
             double tole,
             int stop_count) {
  List output;
  int inn = mi1.nrow();
  int mp1_size = mi1.ncol();
  int mp2_size = mi2.ncol();
  int cure_size = ci.ncol();
  int model_size = mp1_size+mp2_size;
  int allpara_size = model_size+cure_size;
  NumericMatrix est(inn, allpara_size);
  NumericVector l(inn), li(inn);
  LogicalVector Mconverge(inn), Cconverge(inn);
  IntegerVector Mfinal_count(inn), Cfinal_count(inn);
  LogicalVector cand1(inn), cand2(inn), cand3(inn);
  List phase;
  for (int i = 0; i < inn; ++i) {
    NumericVector ini_mp1 = mi1(i,_);
    NumericVector ini_mp2 = mi2(i,_);
    NumericVector ini_cu = ci(i,_);
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
    forwarder_input_model(class_address, model_est);
    NumericVector cure_initial = ini_cu;
    List cure_out;
    if (mtype.compare("mixture") == 0) {
      cure_out = iteration(cure_initial, &forwarder_ite_curerate_mixture_llh, &forwarder_ite_curerate_mixture_score, class_address, tole, stop_count);
    } else if (mtype.compare("non-mixture") == 0) {
      cure_out = iteration(cure_initial, &forwarder_ite_curerate_nmixture_llh, &forwarder_ite_curerate_nmixture_score, class_address, tole, stop_count);
    } else if (mtype.compare("time") == 0) {
      int ce_no = sum(ce);
      if (ce_no == 0) {
        cure_out = iteration(cure_initial, &forwarder_curetime_llh, &forwarder_curetime_score, class_address, tole, stop_count);
      } else {
        int ec_count = 0;
        uvec excess_covariate_index;
        for (uword i = 0; i < ce.size(); i++) {
          if (ce[i] == 1) {
            ec_count++;
            excess_covariate_index.resize(ec_count);
            excess_covariate_index(ec_count-1) = i;
          }
        }
        mat excess_covariate = X.rows(excess_covariate_index);
        vec excess_T = T.elem(excess_covariate_index);
        vec cure_init = as<vec>(cure_initial);
        cure_out = iteration_curetime(cure_init, &forwarder_curetime_llh, &forwarder_curetime_score, excess_covariate, excess_T, class_address, tole, stop_count);
        StringVector pr = as<StringVector>(cure_out["phase"]);
        phase.push_back(pr);
      }
    } else {
      stop("Choose mixture, non-mixture, or time as model type");
    }
    NumericVector cure_est = as<NumericVector>(cure_out["est"]);
    bool cure_converge = as<bool>(cure_out["converge"]);
    int cure_count = as<int>(cure_out["final_count"]);
    double sl = as<double>(cure_out["likelihood"]);
    NumericVector all_est(allpara_size);
    for (int j = 0; j < allpara_size; j++) {
      if (j < model_size) {
        all_est[j] = model_est[j];
      } else {
        int jnd = j-model_size;
        all_est[j] = cure_est[jnd];
      }
    }
    est(i,_) = all_est;
    Mconverge[i] = model_converge;
    Cconverge[i] = cure_converge;
    Mfinal_count[i] = model_count;
    Cfinal_count[i] = cure_count;
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
    cand1[i] = check_all_estl & model_converge & cure_converge;
    cand2[i] = check_all_estl & ((!model_converge & cure_converge) | (model_converge & !cure_converge));
    cand3[i] = check_all_estl & !model_converge & !cure_converge;
  }
  output["phase"] = phase;
  DataFrame summary = DataFrame::create(Named("est") = est,
                                        Named("model_converge") = Mconverge,
                                        Named("cure_converge") = Cconverge,
                                        Named("model_count") = Mfinal_count,
                                        Named("cure_count") = Cfinal_count,
                                        Named("smoothed_likelihood") = l,
                                        Named("likelihood") = li,
                                        Named("candidate1") = cand1,
                                        Named("candidate2") = cand2,
                                        Named("candidate3") = cand3);
  output["summary"] = summary;
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
  } else if (is_true(any(cand3))) {
    for (int i = 0; i < inn; i++) {
      if (cand3[i]) {
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
        out[i] = mi1(0,i);
      } else if ((i >= mp1_size)&&(i < model_size)) {
        int jnd = i-mp1_size;
        out[i] = mi2(0,jnd);
      } else {
        int jnd = i-model_size;
        out[i] = ci(0,jnd);
      }
    }
  }
  output["est"] = out;
  return output;
}
