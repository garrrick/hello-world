#include "fun.h"
#include "iteration.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List check_fun_h(std::string dist_name,
                 List mydata,
                 NumericVector beta_par1,
                 NumericVector beta_par2,
                 NumericVector beta_curetime,
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
  List out;
  NumericVector ee(obs.size());
  for(int i = 0; i < obs.size(); i++){
    ee[i] = 1;
  }
  surv tes;
  tes.base(dist_name,
           FALSE,
           approx_weight,
           ll_initial,
           obs,
           ee,
           cen1,
           cen2,
           cen3,
           cen4,
           gh,
           p1x,
           p2x,
           cxx);
  NumericVector yy(beta_par1.size()+beta_par2.size());
  for (int i = 0; i < beta_par1.size()+beta_par2.size(); i++) {
    if (i < beta_par1.size()) {
      yy[i] = 100;
    } else {
      int jnd = i-beta_par1.size();
      yy[i] = -100;
    }
  }
  double ml = tes.model_llh(yy);
  out["ml"] = ml;
  NumericVector ms = tes.model_score(yy);
  out["ms"] = ms;
  tes.input_model(yy);
  double cl = tes.curetime_llh(beta_curetime);
  out["cl"] = cl;
  NumericVector cs = tes.curetime_score(beta_curetime);
  out["cs"] = cs;
  double icl = tes.ite_curetime_llh(beta_curetime);
  out["icl"] = icl;
  NumericVector ics = tes.ite_curetime_score(beta_curetime);
  out["ics"] = ics;
  tes.input_cure(beta_curetime);
  double iml = tes.ite_model_llh(yy);
  out["iml"] = iml;
  NumericVector ims = tes.ite_model_score(yy);
  out["ims"] = ims;
  tes.input_model(yy);
  double icml = tes.ite_curerate_mixture_llh(beta_curetime);
  out["icml"] = icml;
  NumericVector icms = tes.ite_curerate_mixture_score(beta_curetime);
  out["icms"] = icms;
  double incl = tes.ite_curerate_nmixture_llh(beta_curetime);
  out["incl"] = incl;
  NumericVector incs = tes.ite_curerate_nmixture_score(beta_curetime);
  out["incs"] = incs;
  int coun = beta_par1.size()+beta_par2.size()+beta_curetime.size();
  NumericVector zz(coun);
  for (int i = 0; i < coun; i++) {
    if (i < beta_par1.size()) {
      zz[i] = beta_par1[i];
    } else if ((i >= beta_par1.size())&(i < beta_par1.size()+beta_par2.size())) {
      int jnd = i-beta_par1.size();
      zz[i] = beta_par2[jnd];
    } else {
      int jnd = i-beta_par1.size()-beta_par2.size();
      zz[i] = beta_curetime[jnd];
    }
  }
  double l = tes.llh(zz);
  out["l"] = l;
  NumericVector s = tes.score(zz);
  out["s"] = s;
  double nml = tes.nmixture_llh(zz);
  out["nml"] = nml;
  NumericVector nms = tes.nmixture_score(zz);
  out["nms"] = nms;
  NumericVector ww(201);
  for(int i = 0; i < 201; i++){
    ww[i] = i*.1-10;
  }
  int qq = beta_curetime.size()-1;
  NumericVector c2(qq);
  for(int i = 0; i < qq; i++){
    c2[i] = beta_curetime[i+1];
  }
  return out;
}
