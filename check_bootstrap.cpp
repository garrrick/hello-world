#include "fun.h"
#include "iteration.h"
#include "iteration_curetime.h"
#include "initial.h"
#include "GD.h"
#include "ite_GD.h"
#include "GS.h"
#include "ite_GS.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List check_bootstrap(std::string dist_name,
                     std::string model_name,
                     int boo_times,
                     List mydata,
                     NumericMatrix par1,
                     NumericMatrix par2,
                     NumericMatrix cu,
                     double genp1,
                     double genp2,
                     double cenp1,
                     double cenp2,
                     double ll_initial,
                     double approx_weight,
                     double tolerance,
                     int stopCount) {
  List output;
  NumericVector obs = as<NumericVector>(mydata["obs"]);
  NumericVector cen1 = as<NumericVector>(mydata["cen1"]);
  NumericVector cen2 = as<NumericVector>(mydata["cen2"]);
  NumericVector cen3 = as<NumericVector>(mydata["cen3"]);
  NumericVector cen4 = as<NumericVector>(mydata["cen4"]);
  NumericVector gh = as<NumericVector>(mydata["gh"]);
  NumericMatrix p1x = as<NumericMatrix>(mydata["p1x"]);
  NumericMatrix p2x = as<NumericMatrix>(mydata["p2x"]);
  NumericMatrix cxx = as<NumericMatrix>(mydata["cxx"]);
  int sample_size = obs.size();
  int par1_size = par1.ncol();
  int par2_size = par2.ncol();
  int cu_size = cu.ncol();
  int model_size = par1_size+par2_size;
  int allpara_size = model_size+cu_size;
  double cen1_rate = mean(cen1);
  double cen2_rate = mean(cen2);
  double cen3_rate = mean(cen3);
  double cen4_rate = mean(cen4);
  output["censoring_rate"] = cen1_rate;
  output["To_rate"] = cen2_rate;
  output["Te_rate"] = cen3_rate;
  output["T_rate"] = cen4_rate;

  NumericMatrix boo_parameter(boo_times, allpara_size), boo_index(boo_times, sample_size);
  NumericVector boo_l(boo_times);
  LogicalVector boo_noconv(boo_times);
  List mfc, cfc;
  vec ot = as<vec>(obs);
  vec c1 = as<vec>(cen1);
  vec c2 = as<vec>(cen2);
  vec c3 = as<vec>(cen3);
  vec c4 = as<vec>(cen4);
  vec ghaz = as<vec>(gh);
  mat p1 = as<mat>(p1x);
  mat p2 = as<mat>(p2x);
  mat cx = as<mat>(cxx);
  for (int i = 0; i < boo_times; i++) {
    NumericVector b_ind = rbeta(sample_size, sqrt(2)-1, 1);
    NumericVector booindex_1 = runif(sample_size);
    NumericVector booindex_2(sample_size);
    for (int j = 0; j < sample_size; j++) {
      booindex_2[j] = booindex_1[j]*(sample_size-1);
    }
    NumericVector booindex_3 = round(booindex_2, 0);
    uvec bindex = as<uvec>(booindex_3);
    vec bt = ot.elem(bindex);
    vec b1 = c1.elem(bindex);
    vec b2 = c2.elem(bindex);
    vec b3 = c3.elem(bindex);
    vec b4 = c4.elem(bindex);
    vec bghaz = ghaz.elem(bindex);
    mat bp1 = p1.rows(bindex);
    mat bp2 = p2.rows(bindex);
    mat bcx = cx.rows(bindex);
    NumericVector bobs = wrap(bt);
    NumericVector bcen1 = wrap(b1);
    NumericVector bcen2 = wrap(b2);
    NumericVector bcen3 = wrap(b3);
    NumericVector bcen4 = wrap(b4);
    NumericVector bgh = wrap(bghaz);
    NumericMatrix bp1x = wrap(bp1);
    NumericMatrix bp2x = wrap(bp2);
    NumericMatrix bcxx = wrap(bcx);
    surv B;
    B.base(dist_name,
           FALSE,
           approx_weight,
           ll_initial,
           bobs,
           b_ind,
           bcen1,
           bcen2,
           bcen3,
           bcen4,
           bgh,
           bp1x,
           bp2x,
           bcxx);
    List boo_out;
    if (cu_size == 1) {
      boo_out = ite_GS(model_name, par1, par2, -5, 5, &B, tolerance, stopCount);
    } else {
      boo_out = ite_GD(model_name, par1, par2, cu, bcx, bt, bcen3, &B, tolerance, stopCount);
    }
    NumericVector bp = as<NumericVector>(boo_out["est"]);
    boo_parameter(i,_) = bp;
    boo_index(i,_) = b_ind;
    boo_l[i] = B.llh(bp);
    List boo_mfc = as<List>(boo_out["model_iterations_history"]);
    mfc.push_back(boo_mfc);
    List boo_cfc = as<List>(boo_out["cure_iterations_history"]);
    cfc.push_back(boo_cfc);
    bool noconv = as<bool>(boo_out["no_converge"]);
    boo_noconv[i] = noconv;
  }

  output["boo_est"] = boo_parameter;
  NumericVector out_sd(allpara_size);
  for(int i = 0; i < allpara_size; i++) {
    NumericVector para = boo_parameter(_,i);
    double para_sd = sd(para);
    out_sd[i] = para_sd;
  }
  output["sd"] = out_sd;
  output["bootstrap_index"] = boo_index;
  output["l"] = boo_l;
  output["mfc"] = mfc;
  output["cfc"] = cfc;
  output["no_converge"] = boo_noconv;
  return output;
}
