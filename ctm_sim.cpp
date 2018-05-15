#include "fun.h"
#include "iteration.h"
#include "initial.h"
#include "GD.h"
#include "ite_GD.h"
#include "GS.h"
#include "ite_GS.h"
#include "cen_est.h"
#include "parametric_bootstrap_sample.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List ctm_sim(std::string dist_name,
             std::string model_name,
             std::string cen_model_name,
             std::string est_method,
             int boo_times,
             List mydata,
             NumericMatrix par1,
             NumericMatrix par2,
             NumericMatrix cu,
             NumericVector cenp1,
             NumericVector cenp2,
             double genp1,
             double genp2,
             int iic,
             int eic,
             int bic,
             double ll_initial,
             NumericVector approx_weight,
             double tolerance,
             int stopCount) {
  List output;
  NumericVector obs = as<NumericVector>(mydata["obs"]);
  vec ct = as<vec>(obs);
  NumericVector cen1 = as<NumericVector>(mydata["cen1"]);
  NumericVector cen2 = as<NumericVector>(mydata["cen2"]);
  NumericVector cen3 = as<NumericVector>(mydata["cen3"]);
  NumericVector cen4 = as<NumericVector>(mydata["cen4"]);
  NumericVector gh = as<NumericVector>(mydata["gh"]);
  NumericMatrix p1x = as<NumericMatrix>(mydata["p1x"]);
  NumericMatrix p2x = as<NumericMatrix>(mydata["p2x"]);
  NumericMatrix cxx = as<NumericMatrix>(mydata["cxx"]);
  mat cx = as<mat>(cxx);
  int sample_size = obs.size();
  int par1_size = par1.ncol();
  int par2_size = par2.ncol();
  int cu_size = cu.ncol();
  int model_size = par1_size+par2_size;
  int allpara_size = model_size+cu_size;
  int a_size = approx_weight.size();
  double cen1_rate = mean(cen1);
  double cen2_rate = mean(cen2);
  double cen3_rate = mean(cen3);
  double cen4_rate = mean(cen4);
  output["censoring_rate"] = cen1_rate;
  output["To_rate"] = cen2_rate;
  output["Te_rate"] = cen3_rate;
  output["T_rate"] = cen4_rate;

  NumericMatrix initial_parameter(a_size, allpara_size), est_parameter(a_size, allpara_size);
  NumericVector alike(a_size);
  LogicalVector noconv(a_size);
  for (int i = 0; i < a_size; i++) {
    double aw = approx_weight[i];
    NumericVector b_ind(sample_size);
    for (int j = 0; j < sample_size; j++) {
      b_ind[j] = 1;
    }
    surv W;
    W.base(dist_name,
           FALSE,
           aw,
           ll_initial,
           obs,
           b_ind,
           cen1,
           cen2,
           cen3,
           cen4,
           gh,
           p1x,
           p2x,
           cxx);
    List initial_out = initial(model_name, par1, par2, cu, cx, ct, cen3, &W, tolerance, stopCount);
    NumericVector ip = as<NumericVector>(initial_out["est"]);
    initial_parameter(i,_) = ip;
    NumericVector init_par1(par1_size), init_par2(par2_size), init_cu(cu_size);
    for (int j = 0; j < allpara_size; j++) {
      if (j < par1_size) {
        init_par1[j] = ip[j];
      } else if ((j >= par1_size) && (j < model_size)) {
        int jnd = j - par1_size;
        init_par2[jnd] = ip[j];
      } else {
        int jnd = j - model_size;
        init_cu[jnd] = ip[j];
      }
    }

    NumericMatrix e_par1(iic, par1_size), e_par2(iic, par2_size), e_cu(iic, cu_size);
    for (int j = 0; j < eic; j++) {
      if (j == 0) {
        e_par1(j,_) = init_par1;
        e_par2(j,_) = init_par2;
        e_cu(j,_) = init_cu;
      } else if (j == 1) {
        e_par1(j,_) = par1(0,_);
        e_par2(j,_) = par2(0,_);
        e_cu(j,_) = cu(0,_);
      } else {
        for (int ii = 0; ii < par1_size; ii++) {
          NumericVector mipar1 = runif(1, -1, 1);
          e_par1(j,ii) = par1(0,ii)+mipar1[0];
        }
        for (int ii = 0; ii < par2_size; ii++) {
          NumericVector mipar2 = runif(1, -3, 3);
          e_par2(j,ii) = par2(0,ii)+mipar2[0];
        }
        for (int ii = 0; ii < cu_size; ii++) {
          e_cu(j,ii) = init_cu[ii];
        }
      }
    }
    List est_out;
    if (est_method.compare("ite") == 0) {
      if (cu_size == 1) {
        est_out = ite_GS(model_name, e_par1, e_par2, -5, 5, &W, tolerance, stopCount);
      } else {
        est_out = ite_GD(model_name, e_par1, e_par2, e_cu, cx, ct, cen3, &W, tolerance, stopCount);
      }
    } else if (est_method.compare("GD") == 0) {
      est_out = GD(model_name, e_par1, e_par2, e_cu, &W, tolerance, stopCount);
    }
    NumericVector ep = as<NumericVector>(est_out["est"]);
    est_parameter(i,_) = ep;
    bool nconv = as<bool>(est_out["no_converge"]);
    noconv[i] = nconv;
    double a_l = as<double>(est_out["loglikelihood"]);
    alike[i] = a_l;
  }
  int aindex = which_max(alike);
  NumericVector out_est = est_parameter(aindex,_);
  output["est"] = out_est;
  double out_a = approx_weight[aindex];
  output["approx_weight"] = out_a;
  bool out_noconv = noconv[aindex];
  output["no_converge"] = out_noconv;

  NumericVector est_par1(par1_size), est_par2(par2_size), est_cu(cu_size);
  for (int j = 0; j < allpara_size; j++) {
    if (j < par1_size) {
      est_par1[j] = out_est[j];
    } else if ((j >= par1_size) && (j < model_size)) {
      int jnd = j - par1_size;
      est_par2[jnd] = out_est[j];
    } else {
      int jnd = j - model_size;
      est_cu[jnd] = out_est[j];
    }
  }

  int cen_para_size = allpara_size-2;
  NumericMatrix c_par1(bic, 1), c_par2(bic, cen_para_size);
  for (int j = 0; j < bic; j++) {
    if (j == 0) {
      c_par1(j,0) = 0;
      for (int ii = 0; ii < cen_para_size; ii++) {
        c_par2(j,ii) = 0;
      }
    } else {
      double cii1 = R::runif(-1, 1);
      c_par1(j,0) = cii1;
      NumericVector cii2 = runif(cen_para_size, -1, 1);
      c_par2(j,_) = cii2;
    }
  }
  
  NumericVector b_ind = runif(sample_size, 0, 1);
  NumericVector cen5(sample_size);
  NumericMatrix cen1x(sample_size, 1), cen2x(sample_size, cen_para_size);
  for(int i = 0; i < sample_size; i++){
    cen5[i] = 1-cen1[i];
    cen1x(i,0) = 1;
    for (int j = 0; j < cen_para_size; j++) {
      if (j >= 0 && j < par1_size) {
        cen2x(i, j) = p1x(i, j);
      } else if (j >= par1_size && j < (model_size-1)) {
        int jnd = j-par1_size+1;
        cen2x(i, j) = p2x(i, jnd);
      } else {
        int jnd = j-model_size+2;
        cen2x(i, j) = cxx(i, jnd);
      }
    }
  }
  surv C;
  C.base(dist_name,
         FALSE,
         out_a,
         ll_initial,
         obs,
         b_ind,
         cen5,
         cen2,
         cen3,
         cen4,
         gh,
         cen1x,
         cen2x,
         cxx);
  List ecen = cen_est(cen_model_name, c_par1, c_par2, &C, tolerance, stopCount);
  NumericVector bc = as<NumericVector>(ecen["est"]);
  int c1_size = cen1x.ncol();
  int c2_size = cen2x.ncol();
  NumericVector c1(c1_size);
  NumericVector c2(c2_size);
  for (int i = 0; i < bc.size(); i++) {
    if (i < c1_size) {
      c1[i] = bc[i];
    } else {
      int jnd = i-c1_size;
      c2[jnd] = bc[i];
    }
  }
  NumericMatrix b_par1(bic, par1_size), b_par2(bic, par2_size), b_cu(bic, cu_size);
  for (int j = 0; j < bic; j++) {
    if (j == 0) {
      b_par1(j,_) = est_par1;
      b_par2(j,_) = est_par2;
      b_cu(j,_) = est_cu;
    } else {
      for (int ii = 0; ii < par1_size; ii++) {
        NumericVector mipar1 = runif(1, -1, 1);
        b_par1(j,ii) = est_par1[ii]+mipar1[0];
      }
      for (int ii = 0; ii < par2_size; ii++) {
        NumericVector mipar2 = runif(1, -3, 3);
        b_par2(j,ii) = est_par2[ii]+mipar2[0];
      }
      for (int ii = 0; ii < cu_size; ii++) {
        b_cu(j,ii) = est_cu[ii];
      }
    }
  }
  /*testing for censoring estimation*/
  NumericMatrix boo_parameter(boo_times, allpara_size);
  NumericMatrix testx(sample_size, 1);
  for (int i = 0; i < sample_size; i++) {
    testx(i,0) = 1;
  }
  /*testing end*/
  double c4_of_c234 = cen4_rate/(1-cen1_rate);
  NumericMatrix bcen_percentage(boo_times, 4);
  for (int i = 0; i < boo_times; i++) {
    bool StopFlag = TRUE;
    int count = 0;
    while (StopFlag) {
      count++;
      List bdata = parametric_bootstrap_sample(dist_name, out_est, p1x, p2x, cxx, testx, testx, cenp1, cenp2, genp1, genp2, c4_of_c234);
      NumericVector bobs = as<NumericVector>(bdata["bobs"]);
      NumericVector bcen1 = as<NumericVector>(bdata["bcen1"]);
      NumericVector bcen2 = as<NumericVector>(bdata["bcen2"]);
      NumericVector bcen3 = as<NumericVector>(bdata["bcen3"]);
      NumericVector bcen4 = as<NumericVector>(bdata["bcen4"]);
      NumericVector bgh = as<NumericVector>(bdata["bgh"]);
      NumericVector b_ind = runif(sample_size, 0, 1);
      surv B;
      B.base(dist_name,
             FALSE,
             out_a,
             ll_initial,
             bobs,
             b_ind,
             bcen1,
             bcen2,
             bcen3,
             bcen4,
             bgh,
             p1x,
             p2x,
             cxx);
      List boo_out;
      if (est_method.compare("ite") == 0) {
        if (cu_size == 1) {
          boo_out = ite_GS(model_name, b_par1, b_par2, -5, 5, &B, tolerance, stopCount);
        } else {
          vec bt = as<vec>(bobs);
          boo_out = ite_GD(model_name, b_par1, b_par2, b_cu, cx, bt, bcen3, &B, tolerance, stopCount);
        }
      } else if (est_method.compare("GD") == 0) {
        boo_out = GD(model_name, b_par1, b_par2, b_cu, &B, tolerance, stopCount);
      }
      bool bconv = as<bool>(boo_out["converge"]);
      if (bconv) {
        StopFlag = FALSE;
        double mean_bcen1 = mean(bcen1);
        double mean_bcen2 = mean(bcen2);
        double mean_bcen3 = mean(bcen3);
        double mean_bcen4 = mean(bcen4);
        bcen_percentage(i,0) = mean_bcen1;
        bcen_percentage(i,1) = mean_bcen2;
        bcen_percentage(i,2) = mean_bcen3;
        bcen_percentage(i,3) = mean_bcen4;
        NumericVector bp = as<NumericVector>(boo_out["est"]);
        boo_parameter(i,_) = bp;
      } else if (count == stopCount) {
        StopFlag = FALSE;
        double mean_bcen1 = mean(bcen1);
        double mean_bcen2 = mean(bcen2);
        double mean_bcen3 = mean(bcen3);
        double mean_bcen4 = mean(bcen4);
        bcen_percentage(i,0) = mean_bcen1;
        bcen_percentage(i,1) = mean_bcen2;
        bcen_percentage(i,2) = mean_bcen3;
        bcen_percentage(i,3) = mean_bcen4;
        NumericVector bp = as<NumericVector>(boo_out["est"]);
        boo_parameter(i,_) = bp;
        warning("no converge in %i-th bootstrapping");
      }
    }
  }

  output["bootstrap_cen_rate"] = bcen_percentage;
  output["boo_est"] = boo_parameter;
  NumericVector out_sd(allpara_size);
  for(int i = 0; i < allpara_size; i++) {
    NumericVector para = boo_parameter(_,i);
    double para_sd = sd(para);
    out_sd[i] = para_sd;
  }
  output["sd"] = out_sd;

  return output;
}
