#include "parametric_bootstrap_sample.h"
#include <Rcpp.h>

using namespace Rcpp;

List parametric_bootstrap_sample(std::string dist,
                                 NumericVector y,
                                 NumericMatrix p1x,
                                 NumericMatrix p2x,
                                 NumericMatrix cxx,
                                 NumericMatrix cen1x,
                                 NumericMatrix cen2x,
                                 NumericVector cp1,
                                 NumericVector cp2,
                                 double gp1,
                                 double gp2,
                                 double c4_of_c234) {
  int sam_size = cxx.nrow();
  int allpara_size = y.size();
  int mp1_size = p1x.ncol();
  int mp2_size = p2x.ncol();
  int cure_size = cxx.ncol();
  int cp1_size = cen1x.ncol();
  int cp2_size = cen2x.ncol();
  int cen_size = cp1_size+cp2_size;

  NumericVector mp1(mp1_size), mp2(mp2_size), mc(cure_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else if ((i >= mp1_size) && (i < mp1_size+mp2_size)) {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    } else {
      int jnd = i-mp1_size-mp2_size;
      mc[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += p1x(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += p2x(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cxx(i,j)*mc[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector c1(sam_size), c2(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double cc1 = 0, cc2 = 0;
    for (int j = 0; j < cp1_size; j++) {
      cc1 += cen1x(i,j)*cp1[j];
    }
    for (int j = 0; j < cp2_size; j++) {
      cc2 += cen2x(i,j)*cp2[j];
    }
    c1[i] = cc1;
    c2[i] = cc2;
  }
  
  NumericVector ep1(sam_size), ep2(sam_size), ec(sam_size);
  NumericVector be1(sam_size), be2(sam_size), bc1(sam_size);
  LogicalVector mix1(sam_size), mix2(sam_size);
  if (dist.compare("Weibull") == 0) {
    be2 = rweibull(sam_size, gp1, gp2);
    for (int i = 0; i < sam_size; i++) {
      ep1[i] = exp(m1[i]);
      ep2[i] = exp(m2[i]);
      ec[i] = exp(m3[i]);
      be1[i] = R::rweibull(ep1[i], ep2[i]);
      mix1[i] = be1[i] <= ec[i];
      mix2[i] = be1[i] < be2[i];
      double bcp1 = exp(c1[i]);
      double bcp2 = exp(c2[i]);
      bc1[i] = R::rweibull(bcp1, bcp2);
    }
  } else if (dist.compare("log-normal") == 0) {
    be2 = rlnorm(sam_size, gp1, gp2);
    for (int i = 0; i < sam_size; i++) {
      ep1[i] = m1[i];
      ep2[i] = m2[i];
      ec[i] = exp(m3[i]);
      be1[i] = R::rlnorm(ep1[i], ep2[i]);
      mix1[i] = be1[i] <= ec[i];
      mix2[i] = be1[i] < be2[i];
      double bcp1 = c1[i];
      double bcp2 = c2[i];
      bc1[i] = R::rlnorm(bcp1, bcp2);
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double u2 = R::runif(0,1);
      be2[i] = gp2*pow((1-u2)/u2,-1/gp1);
      ep1[i] = exp(m1[i]);
      ep2[i] = exp(m2[i]);
      ec[i] = exp(m3[i]);
      double u1 = R::runif(0,1);
      be1[i] = ep2[i]*pow((1-u1)/u1,-1/ep1[i]);
      mix1[i] = be1[i] <= ec[i];
      mix2[i] = be1[i] < be2[i];
      double bcp1 = exp(c1[i]);
      double bcp2 = exp(c2[i]);
      double u3 = R::runif(0,1);
      bc1[i] = bcp2*pow((1-u3)/u3,-1/bcp1);
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector be12 = ifelse(mix2, be1, be2);
  NumericVector be3 = ifelse(mix1, be12, be2);

  NumericVector bcen1(sam_size), bcen2(sam_size), bcen3(sam_size), bcen4(sam_size), bobs(sam_size);
  for (int i = 0; i < sam_size; i++) {
    if (bc1[i] < be3[i]) {
      bcen1[i] = 1;
    } else {
      bcen1[i] = 0;
    }
    if (bcen1[i] == 0) {
      bcen4[i] = R::rbinom(1, c4_of_c234);
    } else {
      bcen4[i] = 0;
    }
    if (bcen4[i] == 0 && bcen1[i] == 0) {
      if (be3[i] == be2[i]) {
        bcen2[i] = 1;
      } else {
        bcen3[i] = 1;
      }
    } else {
      bcen2[i] = 0;
      bcen3[i] = 0;
    }
    bobs[i] = bcen1[i]*bc1[i]+bcen2[i]*be2[i]+bcen3[i]*be1[i]+bcen4[i]*be3[i];
  }

  NumericVector bgcdf(sam_size), bgpdf(sam_size), bgh(sam_size);
  if (dist.compare("Weibull") == 0) {
    bgcdf = pweibull(bobs, gp1, gp2);
    bgpdf = dweibull(bobs, gp1, gp2);
  } else if (dist.compare("log-normal") == 0) {
    bgcdf = plnorm(bobs, gp1, gp2);
    bgpdf = dlnorm(bobs, gp1, gp2);
  } else if (dist.compare("log-logistic") == 0){
    for (int i = 0; i < sam_size; i++) {
      bgcdf[i] = 1-1/(1+pow(bobs[i]/gp2, gp1));
      bgpdf[i] = ((gp1/gp2)*pow(bobs[i]/gp2, gp1-1))/pow(1+pow(bobs[i]/gp2, gp1), 2);
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }
  for (int i = 0; i < sam_size; i++) {
    bgh[i] = bgpdf[i]/(1-bgcdf[i]);
  }

  List out;
  out["bobs"] = bobs;
  out["bcen1"] = bcen1;
  out["bcen2"] = bcen2;
  out["bcen3"] = bcen3;
  out["bcen4"] = bcen4;
  out["bgh"] = bgh;

  return out;
}
