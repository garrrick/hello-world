#include "fun.h"
#include "iteration.h"
#include "iteration_curetime.h"
#include "ite_GD.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

List ite_GD(std::string mtype,
            NumericMatrix y0,
            NumericMatrix y1,
            NumericMatrix y2,
            arma::mat X,
            arma::vec T,
            NumericVector ce,
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
  LogicalVector Mconverge(inn), Cconverge(inn);
  IntegerVector Mfinal_count(inn), Cfinal_count(inn);
  NumericVector l(inn), li(inn);
  LogicalVector cand1(inn), cand2(inn), cand3(inn);
  IntegerVector Count(inn);
  List mfc, cfc, ac;
  for (int i = 0; i < inn; ++i) {
    NumericVector a0(allpara_size);
    for (int j = 0; j < allpara_size; j++) {
      if (j < mp1_size) {
        a0[j] = y0(i,j);
      } else if ((j >= mp1_size) && (j < model_size)) {
        int jnd = j-mp1_size;
        a0[j] = y1(i,jnd);
      } else {
        int jnd = j-model_size;
        a0[j] = y2(i,jnd);
      }
    }
    LogicalVector check_a0_na = !is_na(a0);
    LogicalVector check_a0_nan = !is_nan(a0);
    LogicalVector check_a0_inf = !is_infinite(a0);
    LogicalVector check_a0 = check_a0_na & check_a0_nan & check_a0_inf;
    NumericVector l_value;
    if (mtype.compare("mixture") == 0) {
      l_value[0] = forwarder_mixture_llh(class_address, a0);
    } else if (mtype.compare("non-mixture") == 0) {
      l_value[0] = forwarder_nmixture_llh(class_address, a0);
    } else if (mtype.compare("time") == 0) {
      l_value[0] = forwarder_llh(class_address, a0);
    } else {
      stop("Choose mixture, non-mixture, or time as model type");
    }
    LogicalVector check_a0_l_na = !is_na(l_value);
    LogicalVector check_a0_l_nan = !is_nan(l_value);
    LogicalVector check_a0_l_inf = !is_infinite(l_value);
    LogicalVector check_a0_l = check_a0_l_na & check_a0_l_nan & check_a0_l_inf;
    NumericVector s_values(allpara_size);
    if (mtype.compare("mixture") == 0) {
      s_values = forwarder_mixture_score(class_address, a0);
    } else if (mtype.compare("non-mixture") == 0) {
      s_values = forwarder_nmixture_score(class_address, a0);
    } else if (mtype.compare("time") == 0) {
      s_values = forwarder_score(class_address, a0);
    }
    LogicalVector check_a0_s_na = !is_na(s_values);
    LogicalVector check_a0_s_nan = !is_nan(s_values);
    LogicalVector check_a0_s_inf = !is_infinite(s_values);
    LogicalVector check_a0_s = check_a0_s_na & check_a0_s_nan & check_a0_s_inf;
    bool b_check_a0 = is_true(all(check_a0)) && is_true(all(check_a0_l)) && is_true(all(check_a0_s));
    Count[i] = 0;
    Mfinal_count[i] = 0;
    Cfinal_count[i] = 0;
    IntegerVector mcount, ccount;
    bool StopFlag = TRUE;
    if (b_check_a0) {
      while (StopFlag) {
        Count[i]++;
        NumericVector model_initial(model_size), cure_initial(cure_size);
        for (int j = 0; j < allpara_size; j++) {
          if (j < model_size) {
            model_initial[j] = a0[j];
          } else {
            int jnd = j - model_size;
            cure_initial[jnd] = a0[j];
          }
        }
        forwarder_input_cure(class_address, cure_initial);
        List model_out;
        if (mtype.compare("mixture") == 0) {
          model_out = iteration(model_initial, &forwarder_ite_model_mixture_llh, &forwarder_ite_model_mixture_score, class_address, tole, stop_count);
        } else if (mtype.compare("non-mixture") == 0) {
          model_out = iteration(model_initial, &forwarder_ite_model_nmixture_llh, &forwarder_ite_model_nmixture_score, class_address, tole, stop_count);
        } else if (mtype.compare("time") == 0) {
          model_out = iteration(model_initial, &forwarder_ite_model_llh, &forwarder_ite_model_score, class_address, tole, stop_count);
        }
        NumericVector model_est = as<NumericVector>(model_out["est"]);
        bool model_converge = as<bool>(model_out["converge"]);
        int model_count = as<int>(model_out["final_count"]);
        mcount.push_back(model_count);
        Mfinal_count[i]++;
        forwarder_input_model(class_address, model_est);
        List cure_out;
        if (mtype.compare("mixture") == 0) {
          cure_out = iteration(cure_initial, &forwarder_ite_curerate_mixture_llh, &forwarder_ite_curerate_mixture_score, class_address, tole, stop_count);
        } else if (mtype.compare("non-mixture") == 0) {
          cure_out = iteration(cure_initial, &forwarder_ite_curerate_nmixture_llh, &forwarder_ite_curerate_nmixture_score, class_address, tole, stop_count);
        } else if (mtype.compare("time") == 0) {
          int ce_no = sum(ce);
          if (ce_no == 0) {
            cure_out = iteration(cure_initial, &forwarder_ite_curetime_llh, &forwarder_ite_curetime_score, class_address, tole, stop_count);
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
            cure_out = iteration_curetime(cure_init, &forwarder_ite_curetime_llh, &forwarder_ite_curetime_score, excess_covariate, excess_T, class_address, tole, stop_count);
          }
        }
        NumericVector cure_est = as<NumericVector>(cure_out["est"]);
        bool cure_converge = as<bool>(cure_out["converge"]);
        int cure_count = as<int>(cure_out["final_count"]);
        ccount.push_back(cure_count);
        Cfinal_count[i]++;
        NumericVector a1(allpara_size);
        for (int j = 0; j < allpara_size; j++) {
          if (j < model_size) {
            a1[j] = model_est[j];
          } else {
            int jnd = j-model_size;
            a1[j] = cure_est[jnd];
          }
        }
        bool b_check_a1 = model_converge && cure_converge;
        if (b_check_a1) {
          bool check_stop = (Count[i] == stop_count);
          double norm2 = 0, dn = 0;
          for(int j = 0; j < allpara_size; j++){
            norm2 += pow(a0[j]-a1[j], 2);
            dn += pow(a0[j], 2);
          }
          double norm;
          if (dn != 0) {
            norm = sqrt(norm2)/sqrt(dn);
          } else {
            norm = sqrt(norm2)/(sqrt(dn)+.0000001);
          }
          bool check_converge_a1 = (norm <= tole);
          if (check_converge_a1 || check_stop) {
            est(i,_) = a1;
            Mconverge[i] = model_converge;
            Cconverge[i] = cure_converge;
            mfc.push_back(mcount);
            cfc.push_back(ccount);
            if (mtype.compare("mixture") == 0) {
              l[i] = forwarder_mixture_llh(class_address, a1);
              li[i] = l[i];
            } else if (mtype.compare("non-mixture") == 0) {
              l[i] = forwarder_nmixture_llh(class_address, a1);
              li[i] = l[i];
            } else if (mtype.compare("time") == 0) {
              l[i] = forwarder_llh(class_address, a1);
              li[i] = forwarder_llh_indicator(class_address, a1);
            }
            NumericMatrix acc = as<NumericMatrix>(cure_out["active_set"]);
            ac.push_back(acc);
            StopFlag = FALSE;
          } else {
            a0 = a1;
          }
        } else {
          est(i,_) = a0;
          Mconverge[i] = FALSE;
          Cconverge[i] = FALSE;
          mfc.push_back(mcount);
          cfc.push_back(ccount);
          if (mtype.compare("mixture") == 0) {
            l[i] = forwarder_mixture_llh(class_address, a0);
            li[i] = l[i];
          } else if (mtype.compare("non-mixture") == 0) {
            l[i] = forwarder_nmixture_llh(class_address, a0);
            li[i] = l[i];
          } else if (mtype.compare("time") == 0) {
            l[i] = forwarder_llh(class_address, a0);
            li[i] = forwarder_llh_indicator(class_address, a0);
          }
          NumericMatrix acc = as<NumericMatrix>(cure_out["active_set"]);
          ac.push_back(acc);
          StopFlag = FALSE;
        }
      }
      NumericVector estl = est(i,_);
      estl.push_back(li[i]);
      LogicalVector check_na_estl = !is_na(estl);
      LogicalVector check_nan_estl = !is_nan(estl);
      LogicalVector check_inf_estl = !is_infinite(estl);
      LogicalVector check_estl = check_na_estl & check_nan_estl & check_inf_estl;
      bool check_all_estl = is_true(all(check_estl));
      cand1[i] = check_all_estl && Mconverge[i] && Cconverge[i];
      cand2[i] = check_all_estl && ((!Mconverge[i] && Cconverge[i]) || (Mconverge[i] && !Cconverge[i]));
      cand3[i] = check_all_estl && !Mconverge[i] && !Cconverge[i];
    } else {
      est(i,_) = a0;
      Mconverge[i] = FALSE;
      Cconverge[i] = FALSE;
      mfc.push_back(0);
      cfc.push_back(0);
      if (mtype.compare("mixture") == 0) {
        l[i] = forwarder_mixture_llh(class_address, a0);
        li[i] = l[i];
      } else if (mtype.compare("non-mixture") == 0) {
        l[i] = forwarder_nmixture_llh(class_address, a0);
        li[i] = l[i];
      } else if (mtype.compare("time") == 0) {
        l[i] = forwarder_llh(class_address, a0);
        li[i] = forwarder_llh_indicator(class_address, a0);
      }
      cand1[i] = FALSE;
      cand2[i] = FALSE;
      cand3[i] = FALSE;
    }
  }
  DataFrame summary = DataFrame::create(Named("est") = est,
                                        Named("model_converge") = Mconverge,
                                        Named("cure_converge") = Cconverge,
                                        Named("model_count") = Mfinal_count,
                                        Named("cure_count") = Cfinal_count,
                                        Named("count") = Count,
                                        Named("smoothed_likelihood") = l,
                                        Named("likelihood") = li,
                                        Named("candidate1") = cand1,
                                        Named("candidate2") = cand2,
                                        Named("candidate3") = cand3);
  output["active_set"] = ac;
  output["summary"] = summary;
  output["model_iterations_history"] = mfc;
  output["cure_iterations_history"] = cfc;
  bool cand1conv = is_true(any(cand1));
  output["converge"] = cand1conv;
  bool noconv = is_true(all(cand3));
  output["no_converge"] = noconv;
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
        out[i] = y0(0,i);
      } else if ((i >= mp1_size) && (i < model_size)) {
        int jnd = i-mp1_size;
        out[i] = y1(0,jnd);
      } else {
        int jnd = i-model_size;
        out[i] = y2(0,jnd);
      }
    }
  }
  NumericVector mout(model_size), cout(cure_size);
  for (int j = 0; j < allpara_size; j++) {
    if (j < model_size) {
      mout[j] = out[j];
    } else {
      int jnd = j-model_size;
      cout[jnd] = out[j];
    }
  }
  output["est"] = out;
  if (mtype.compare("mixture") == 0) {
    output["loglikelihood"] = forwarder_mixture_llh(class_address, out);
    output["score"] = forwarder_mixture_score(class_address, out);
    output["ite_model_score"] = forwarder_ite_model_mixture_score(class_address, mout);
    output["ite_curerate_score"] = forwarder_ite_curerate_mixture_score(class_address, cout);
  } else if (mtype.compare("non-mixture") == 0) {
    output["loglikelihood"] = forwarder_nmixture_llh(class_address, out);
    output["score"] = forwarder_nmixture_score(class_address, out);
    output["ite_model_score"] = forwarder_ite_model_nmixture_score(class_address, mout);
    output["ite_curerate_score"] = forwarder_ite_curerate_nmixture_score(class_address, cout);
  } else if (mtype.compare("time") == 0) {
    output["loglikelihood"] = forwarder_llh_indicator(class_address, out);
    output["sloglikelihood"] = forwarder_llh(class_address, out);
    output["score"] = forwarder_score(class_address, out);
    output["ite_model_score"] = forwarder_ite_model_score(class_address, mout);
    output["ite_curetime_score"] = forwarder_ite_curetime_score(class_address, cout);
  }
  return output;
}
