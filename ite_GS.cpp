#include "fun.h"
#include "iteration.h"
#include "GS.h"
#include "ite_GS.h"
#include <Rcpp.h>

using namespace Rcpp;

List ite_GS(std::string mtype,
            NumericMatrix par1,
            NumericMatrix par2,
            double cu_from,
            double cu_to,
            void* class_address,
            double tole,
            int stop_count) {
  List output;
  int inn = par1.nrow();
  int mp1_size = par1.ncol();
  int mp2_size = par2.ncol();
  int model_size = mp1_size+mp2_size;
  int cure_size = 1;
  int allpara_size = model_size+cure_size;
  NumericMatrix est(inn, allpara_size);
  LogicalVector Mconverge(inn);
  IntegerVector Mfinal_count(inn), Cfinal_count(inn);
  NumericVector l(inn), li(inn);
  LogicalVector cand1(inn), cand2(inn);
  IntegerVector Count(inn);
  for (int i = 0; i < inn; ++i) {
    NumericVector a0(allpara_size);
    for (int j = 0; j < allpara_size; j++) {
      if (j < mp1_size) {
        a0[j] = par1(i,j);
      } else if ((j >= mp1_size) & (j < model_size)) {
        int jnd = j-mp1_size;
        a0[j] = par2(i,jnd);
      } else {
        a0[j] = (cu_from+cu_to)/2;
      }
    }
    LogicalVector check_a0_na = !is_na(a0);
    LogicalVector check_a0_nan = !is_nan(a0);
    LogicalVector check_a0_inf = !is_infinite(a0);
    LogicalVector check_a0 = check_a0_na & check_a0_nan & check_a0_inf;
    NumericVector l_value(1);
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
    bool b_check_a0 = is_true(all(check_a0)) & is_true(all(check_a0_l)) & is_true(all(check_a0_s));
    Count[i] = 0;
    if (b_check_a0) {
      while (TRUE) {
        NumericVector model_initial(model_size);
        for (int j = 0; j < model_size; j++) {
          model_initial[j] = a0[j];
        }
        forwarder_input_model(class_address, model_initial);
        List cure_out;
        if (mtype.compare("mixture") == 0) {
          cure_out = GS(cu_from, cu_to, &forwarder_grid_mixture_llh, class_address, stop_count);
        } else if (mtype.compare("non-mixture") == 0) {
          cure_out = GS(cu_from, cu_to, &forwarder_grid_nmixture_llh, class_address, stop_count);
        } else if (mtype.compare("time") == 0) {
          cure_out = GS(cu_from, cu_to, &forwarder_grid_llh, class_address, stop_count);
        }
        NumericVector cure_est(1);
        cure_est[0] = as<double>(cure_out["est"]);
        int cure_count = as<int>(cure_out["count"]);
        forwarder_input_cure(class_address, cure_est);
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


        NumericVector a1(allpara_size);
        for (int j = 0; j < allpara_size; j++) {
          if (j < model_size) {
            a1[j] = model_est[j];
          } else {
            a1[j] = cure_est[0];
          }
        }
        LogicalVector check_a1_na = !is_na(a1);
        LogicalVector check_a1_nan = !is_nan(a1);
        LogicalVector check_a1_inf = !is_infinite(a1);
        LogicalVector check_a1 = check_a1_na & check_a1_nan & check_a1_inf;
        NumericVector l_value(1);
        if (mtype.compare("mixture") == 0) {
          l_value[0] = forwarder_mixture_llh(class_address, a1);
        } else if (mtype.compare("non-mixture") == 0) {
          l_value[0] = forwarder_nmixture_llh(class_address, a1);
        } else if (mtype.compare("time") == 0) {
          l_value[0] = forwarder_llh(class_address, a1);
        }
        LogicalVector check_a1_l_na = !is_na(l_value);
        LogicalVector check_a1_l_nan = !is_nan(l_value);
        LogicalVector check_a1_l_inf = !is_infinite(l_value);
        LogicalVector check_a1_l = check_a1_l_na & check_a1_l_nan & check_a1_l_inf;
        NumericVector s_values(allpara_size);
        if (mtype.compare("mixture") == 0) {
          s_values = forwarder_mixture_score(class_address, a1);
        } else if (mtype.compare("non-mixture") == 0) {
          s_values = forwarder_nmixture_score(class_address, a1);
        } else if (mtype.compare("time") == 0) {
          s_values = forwarder_score(class_address, a1);
        }
        LogicalVector check_a1_s_na = !is_na(s_values);
        LogicalVector check_a1_s_nan = !is_nan(s_values);
        LogicalVector check_a1_s_inf = !is_infinite(s_values);
        LogicalVector check_a1_s = check_a1_s_na & check_a1_s_nan & check_a1_s_inf;
        bool b_check_a1 = is_true(all(check_a1)) & is_true(all(check_a1_l)) & is_true(all(check_a1_s));
        Count[i] += 1;
        if (b_check_a1) {
          bool check_stop = (Count[i] == stop_count);
          double norm2 = 0, dn = 0;
          for(int j = 0; j < allpara_size; j++){
            norm2 += pow(a0[j]-a1[j], 2);
            dn += pow(a0[j], 2);
          }
          double norm = sqrt(norm2)/(sqrt(dn)+.0001);
          bool check_converge_a1 = (norm <= tole);
          if (check_converge_a1 || check_stop) {
            est(i,_) = a1;
            Mconverge[i] = model_converge;
            Mfinal_count[i] = model_count;
            Cfinal_count[i] = cure_count;
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
            break;
          } else {
            a0 = a1;
          }
        } else {
          est(i,_) = a0;
          Mconverge[i] = FALSE;
          Mfinal_count[i] = model_count;
          Cfinal_count[i] = cure_count;
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
          break;
        }
      }
      NumericVector estl = est(i,_);
      estl.push_back(li[i]);
      LogicalVector check_na_estl = !is_na(estl);
      LogicalVector check_nan_estl = !is_nan(estl);
      LogicalVector check_inf_estl = !is_infinite(estl);
      LogicalVector check_estl = check_na_estl & check_nan_estl & check_inf_estl;
      bool check_all_estl = is_true(all(check_estl));
      cand1[i] = check_all_estl & Mconverge[i];
      cand2[i] = check_all_estl & !Mconverge[i];
    } else {
      est(i,_) = a0;
      Mconverge[i] = FALSE;
      Mfinal_count[i] = 0;
      Cfinal_count[i] = 0;
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
    }
  }
  DataFrame summary = DataFrame::create(Named("est") = est,
                                        Named("model_converge") = Mconverge,
                                        Named("model_count") = Mfinal_count,
                                        Named("cure_count") = Cfinal_count,
                                        Named("count") = Count,
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
        out[i] = par1(0,i);
      } else if ((i >= mp1_size)&&(i < model_size)) {
        int jnd = i-mp1_size;
        out[i] = par2(0,jnd);
      } else {
        out[i] = 1;
      }
    }
  }
  NumericVector mout(model_size), cout(cure_size);
  for (int j = 0; j < allpara_size; j++) {
    if (j < model_size) {
      mout[j] = out[j];
    } else {
      int jnd = j - model_size;
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
    output["score"] = forwarder_score(class_address, out);
    output["ite_model_score"] = forwarder_ite_model_score(class_address, mout);
    output["ite_curetime_score"] = forwarder_ite_curetime_score(class_address, cout);
  }
  return output;
}
