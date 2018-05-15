#ifndef FUN_H
#define FUN_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

class surv {
public:
  void base(std::string,
            bool b_logical,
            double a,
            double s,
            NumericVector lo,
            NumericVector b_index,
            NumericVector cs1,
            NumericVector cs2,
            NumericVector cs3,
            NumericVector cs4,
            NumericVector gh,
            NumericMatrix mp1c,
            NumericMatrix mp2c,
            NumericMatrix cc);
  void input_model(NumericVector);
  void input_cure(NumericVector);
  double model_llh(NumericVector);
  NumericVector model_score(NumericVector);
  double curetime_llh(NumericVector);
  NumericVector curetime_score(NumericVector);
  double llh(NumericVector);
  double llh_indicator(NumericVector);
  NumericVector grid_llh(NumericVector);
  NumericVector grid_mixture_llh(NumericVector);
  NumericVector grid_nmixture_llh(NumericVector);
  double ite_curetime2_llh(NumericVector);
  NumericVector ite_curetime2_score(NumericVector);
  NumericVector score(NumericVector);
  double ite_model_llh(NumericVector);
  NumericVector ite_model_score(NumericVector);
  double ite_curetime_llh(NumericVector);
  NumericVector ite_curetime_score(NumericVector);
  double mixture_llh(NumericVector);
  double nmixture_llh(NumericVector);
  NumericVector mixture_score(NumericVector);
  NumericVector nmixture_score(NumericVector);
  double ite_model_mixture_llh(NumericVector);
  double ite_model_nmixture_llh(NumericVector);
  NumericVector ite_model_mixture_score(NumericVector);
  NumericVector ite_model_nmixture_score(NumericVector);
  double ite_curerate_mixture_llh(NumericVector);
  NumericVector ite_curerate_mixture_score(NumericVector);
  double ite_curerate_nmixture_llh(NumericVector);
  NumericVector ite_curerate_nmixture_score(NumericVector);
  std::string dist;
  bool bl;
  double appro_weight;
  double smooth_para;
  NumericVector bi;
  NumericVector last_obs;
  NumericVector cen_status1;
  NumericVector cen_status2;
  NumericVector cen_status3;
  NumericVector cen_status4;
  NumericVector general_hazard;
  NumericMatrix mp1_covariate;
  NumericMatrix mp2_covariate;
  NumericMatrix cure_covariate;
  NumericVector mp1;
  NumericVector mp2;
  NumericVector cure;
};

inline void surv::base (std::string dis,
                        bool b_logical,
                        double a,
                        double s,
                        NumericVector lo,
                        NumericVector b_index,
                        NumericVector cs1,
                        NumericVector cs2,
                        NumericVector cs3,
                        NumericVector cs4,
                        NumericVector gh,
                        NumericMatrix mp1c,
                        NumericMatrix mp2c,
                        NumericMatrix cc) {
  dist = dis;
  bl = b_logical;
  appro_weight = a;
  smooth_para = s;
  last_obs = lo;
  bi = b_index;
  cen_status1 = cs1;
  cen_status2 = cs2;
  cen_status3 = cs3;
  cen_status4 = cs4;
  general_hazard = gh;
  mp1_covariate = mp1c;
  mp2_covariate = mp2c;
  cure_covariate = cc;
}

inline void surv::input_model(NumericVector y) {
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  NumericVector p1(mp1_size), p2(mp2_size);
  for(int i = 0; i < y.size(); i++) {
    if (i < mp1_size) {
      p1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      p2[jnd] = y[i];
    }
  }
  mp1 = p1;
  mp2 = p2;
}

inline void surv::input_cure(NumericVector y) {
  cure = y;
}

inline double surv::model_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();

  NumericVector par1(mp1_size), par2(mp2_size);
  for(int i = 0; i < y.size(); i++) {
    if (i < mp1_size) {
      par1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      par2[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*par1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*par2[j];
    }
    m1[i] = p1;
    m2[i] = p2;
  }

  NumericVector sp(sam_size), hp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double fp = fn[i]/(g2*last_obs[i]);
      hp[i] = fp/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double inner_mp1 = 0, inner_mp2 = 0;
  for (int i = 0; i < mp1_size; i++) {
    inner_mp1 += pow(par1[i], 2);
  }
  for (int i = 0; i < mp2_size; i++) {
    inner_mp2 += pow(par2[i], 2);
  }

  double penalty = 0.5*smooth_para*(inner_mp1+inner_mp2);

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    individual_loglikelihood[i] = (1-cen_status1[i])*log(hp[i])+log(sp[i]);
  }
  double out = mean(individual_loglikelihood)-penalty;
  return out;
}

inline NumericVector surv::model_score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();

  NumericVector par1(mp1_size), par2(mp2_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      par1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      par2[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*par1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*par2[j];
    }
    m1[i] = p1;
    m2[i] = p2;
  }

  NumericMatrix d_mp1(sam_size, mp1_size), d_mp2(sam_size, mp2_size);
  if ((dist.compare("Weibull") == 0) || (dist.compare("log-logistic") == 0)) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = g1*mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = g2*mp2_covariate(i,j);
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    for (int i = 0; i < sam_size; i++) {
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = mp2_covariate(i,j);
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector sp(sam_size);
  NumericVector hp(sam_size);
  NumericVector fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = hp[i]*sp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1sp(i,j) = -sp2*log(sp1)*sp[i]*d_mp1(i,j);
          d1fp(i,j) = (pow(g1,-1)+log(sp1)-sp2*log(sp1))*fp[i]*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          d1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          d1fp(i,j) = (last_obs[i]*hp[i]-g1)*pow(g2,-1)*fp[i]*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1fp(i,j) = (fp[i]*standard[i]/g2)*d_mp1(i,j);
          d1sp(i,j) = (fn[i]/g2)*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          d1fp(i,j) = (fp[i]/g2)*(pow(standard[i], 2)-1)*d_mp2(i,jnd);
          d1sp(i,j) = ((fn[i]*standard[i])/g2)*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1fp(i,j) = fp[i]*(pow(g1,-1)+log(sp1)-2*last_obs[i]*pow(g1,-1)*log(sp1)*hp[i])*d_mp1(i,j);
          d1sp(i,j) = -last_obs[i]*pow(g1,-1)*log(sp1)*fp[i]*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          d1fp(i,j) = pow(g2,-1)*fp[i]*(-g1+2*last_obs[i]*hp[i])*d_mp1(i,jnd);
          d1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
        }
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector d1p(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    d1p[i] = smooth_para*y[i];
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      individual_score(i,j) = (1-cen_status1[i])*pow(hp[i],-1)*d1hp(i,j) + pow(sp[i],-1)*d1sp(i,j) - d1p[j];
    }
  }

  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = mean(individual_score(_,i));
  }

  return out;
}

inline double surv::curetime_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = y.size();

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*y[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size), spc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = pow(1+spc2, -1);
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double approx_sigma = pow(sam_size, -1/appro_weight);
  NumericVector ind_approx(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    ind_approx[i] = pow(1+exp(-(gc-last_obs[i])*pow(approx_sigma,-1)), -1);
  }

  double inner_cure = 0;
  for (int i = 0; i < cure_size; i++) {
    inner_cure += pow(y[i], 2);
  }
  double penalty = 0.5*smooth_para*inner_cure;

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double th = general_hazard[i]+hp[i]*ind_approx[i];
    double ts = (sp[i]-spc[i])*ind_approx[i]+spc[i];
    double tf = th*ts;
    double fe = fp[i];
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = loglike1+loglike2+loglike3+loglike4-penalty;
  }

  double out = mean(individual_loglikelihood);
  return out;
}

inline NumericVector surv::curetime_score(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = y.size();
  int allpara_size = mp1_size+mp2_size+cure_size;

  NumericVector allpara(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      allpara[i] = mp1[i];
    } else if ((i >= mp1_size) && (i < mp1_size+mp2_size)) {
      int jnd = i-mp1_size;
      allpara[i] = mp2[jnd];
    } else {
      int jnd = i-mp1_size-mp2_size;
      allpara[i] = y[jnd];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*y[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_cr(sam_size, cure_size);

  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    for (int j = 0; j < cure_size; j++) {
      d_cr(i,j) = gc*cure_covariate(i,j);
    }
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  NumericVector spc(sam_size), fpc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = hp[i]*sp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
      fpc[i] = hp1*pow(spc1, g1-1)*spc[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double gc = exp(m3[i]);
      fp[i] = fn[i]/(g2*last_obs[i]);
      fpc[i] = fnc[i]/(g2*gc);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = pow(1+spc2, -1);
      double hpc2 = pow(spc1, g1-1);
      double hpc = hp1*hpc2*spc[i];
      fpc[i] = spc[i]*hpc;
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double approx_sigma = pow(sam_size, -1/appro_weight);
  NumericVector ind_approx(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    ind_approx[i] = pow(1+exp(-(gc-last_obs[i])*pow(approx_sigma,-1)), -1);
  }

  NumericMatrix d1i(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < cure_size; j++) {
      d1i(i,j) = pow(approx_sigma,-1)*ind_approx[i]*(1-ind_approx[i])*d_cr(i,j);
    }
  }

  NumericVector d1p(cure_size);
  for (int i = 0; i < cure_size; i++) {
    d1p[i] = smooth_para*y[i];
  }

  NumericMatrix d1spc(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < cure_size; j++) {
      d1spc(i,j) = -fpc[i]*d_cr(i,j);
    }
  }

  NumericMatrix individual_score(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < cure_size; j++) {
      double th = general_hazard[i]+hp[i]*ind_approx[i];
      double ts = (sp[i]-spc[i])*ind_approx[i]+spc[i];
      double tdh = pow(th,-1)*hp[i]*d1i(i,j);
      double tds = pow(ts,-1)*(-d1spc(i,j)*ind_approx[i]+(sp[i]-spc[i])*d1i(i,j)+d1spc(i,j));
      double tdf = tds+tdh;
      individual_score(i,j) =
        cen_status1[i]*tds+
        cen_status2[i]*tds+
        cen_status4[i]*tdf-d1p[j];
    }
  }

  NumericVector out(cure_size);
  for (int i = 0; i < cure_size; i++) {
    out[i] = mean(individual_score(_,i));
  }

  return out;
}

inline double surv::llh(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size), cr(cure_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else if ((i >= mp1_size) && (i < mp1_size+mp2_size)) {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    } else {
      int jnd = i-mp1_size-mp2_size;
      cr[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cr[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  NumericVector spc(sam_size), fpc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
      fpc[i] = hp1*pow(spc1, g1-1)*spc[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double gc = exp(m3[i]);
      fp[i] = fn[i]/(g2*last_obs[i]);
      fpc[i] = fnc[i]/(g2*gc);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = pow(1+spc2, -1);
      double hpc2 = pow(spc1, g1-1);
      double hpc = hp1*hpc2*spc[i];
      fpc[i] = spc[i]*hpc;
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double approx_sigma = pow(sam_size, -1/appro_weight);
  NumericVector ind_approx(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    ind_approx[i] = pow(1+exp(-(gc-last_obs[i])*pow(approx_sigma,-1)),-1);
  }

  double inner_cure = 0;
  for (int i = 0; i < cure_size; i++) {
    inner_cure += pow(cr[i], 2);
  }
  double penalty = 0.5*smooth_para*inner_cure;

  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double th = general_hazard[i]+hp[i]*ind_approx[i];
    double ts = (sp[i]-spc[i])*ind_approx[i]+spc[i];
    double tf = th*ts;
    double fe = fp[i];
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4-penalty);
  }

  double sum_bi = sum(bi);
  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline double surv::llh_indicator(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size), cr(cure_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else if ((i >= mp1_size) && (i < mp1_size+mp2_size)) {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    } else {
      int jnd = i-mp1_size-mp2_size;
      cr[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cr[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  NumericVector spc(sam_size), fpc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
      fpc[i] = hp1*pow(spc1, g1-1)*spc[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double gc = exp(m3[i]);
      fp[i] = fn[i]/(g2*last_obs[i]);
      fpc[i] = fnc[i]/(g2*gc);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = pow(1+spc2, -1);
      double hpc2 = pow(spc1, g1-1);
      double hpc = hp1*hpc2*spc[i];
      fpc[i] = spc[i]*hpc;
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector ind(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    if (gc - last_obs[i] > 0) {
      ind[i] = 1;
    } else {
      ind[i] = 0;
    }
  }

  double inner_cure = 0;
  for (int i = 0; i < cure_size; i++) {
    inner_cure += pow(cr[i], 2);
  }
  double penalty = 0.5*smooth_para*inner_cure;

  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double th = general_hazard[i]+hp[i]*ind[i];
    double ts = (sp[i]-spc[i])*ind[i]+spc[i];
    double tf = th*ts;
    double fe = fp[i];
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4-penalty);
  }

  double sum_bi = sum(bi);
  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline NumericVector surv::score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size), cr(cure_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else if ((i >= mp1_size) && (i < mp1_size+mp2_size)) {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    } else {
      int jnd = i-mp1_size-mp2_size;
      cr[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cr[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_mp1(sam_size, mp1_size);
  NumericMatrix d_mp2(sam_size, mp2_size);
  NumericMatrix d_cr(sam_size, cure_size);
  if ((dist.compare("Weibull") == 0) || (dist.compare("log-logistic") == 0)) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = g1*mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = g2*mp2_covariate(i,j);
      }
      for (int j = 0; j < cure_size; j++) {
        d_cr(i,j) = gc*cure_covariate(i,j);
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double gc = exp(m3[i]);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = mp2_covariate(i,j);
      }
      for (int j = 0; j < cure_size; j++) {
        d_cr(i,j) = gc*cure_covariate(i,j);
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  NumericVector spc(sam_size), fpc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
      fpc[i] = hp1*pow(spc1, g1-1)*spc[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double gc = exp(m3[i]);
      fp[i] = fn[i]/(g2*last_obs[i]);
      fpc[i] = fnc[i]/(g2*gc);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      double hpc2 = pow(spc1, g1-1);
      spc[i] = pow(1+spc2, -1);
      double hpc = hp1*hpc2*spc[i];
      fpc[i] = hpc*spc[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double approx_sigma = pow(sam_size, -1/appro_weight);
  NumericVector ind_approx(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    ind_approx[i] = pow(1+exp(-(gc-last_obs[i])*pow(approx_sigma,-1)), -1);
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  NumericMatrix d1spc(sam_size, allpara_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1sp(i,j) = -sp2*log(sp1)*sp[i]*d_mp1(i,j);
          d1fp(i,j) = (pow(g1,-1)+log(sp1)-sp2*log(sp1))*fp[i]*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = -spc2*log(spc1)*spc[i]*d_mp1(i,j);
        } else if ((j >= mp1_size) && (j < mp1_size+mp2_size)) {
          int jnd = j-mp1_size;
          d1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          d1fp(i,j) = (last_obs[i]*hp[i]-g1)*pow(g2,-1)*fp[i]*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = spc1*fpc[i]*d_mp2(i,jnd);
        } else {
          int jnd = j-mp1_size-mp2_size;
          d1sp(i,j) = 0;
          d1fp(i,j) = 0;
          d1hp(i,j) = 0;
          d1spc(i,j) = -fpc[i]*d_cr(i,jnd);
        }
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1fp(i,j) = (fp[i]*standard[i]/g2)*d_mp1(i,j);
          d1sp(i,j) = (fn[i]/g2)*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = (fnc[i]/g2)*d_mp1(i,j);
        } else if ((j >= mp1_size) && (j < mp1_size+mp2_size)) {
          int jnd = j-mp1_size;
          d1fp(i,j) = (fp[i]/g2)*(pow(standard[i], 2)-1)*d_mp2(i,jnd);
          d1sp(i,j) = ((fn[i]*standard[i])/g2)*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = ((fnc[i]*standardc[i])/g2)*d_mp2(i,jnd);
        } else {
          int jnd = j-mp1_size-mp2_size;
          d1sp(i,j) = 0;
          d1fp(i,j) = 0;
          d1hp(i,j) = 0;
          d1spc(i,j) = -fpc[i]*d_cr(i,jnd);
        }
      }
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double spc1 = gc/g2;
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1fp(i,j) = fp[i]*(pow(g1,-1)+log(sp1)-2*last_obs[i]*pow(g1,-1)*log(sp1)*hp[i])*d_mp1(i,j);
          d1sp(i,j) = -last_obs[i]*pow(g1,-1)*log(sp1)*fp[i]*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = -gc*pow(g1,-1)*log(spc1)*fpc[i]*d_mp1(i,j);
        } else if ((j >= mp1_size) && (j < mp1_size+mp2_size)) {
          int jnd = j-mp1_size;
          d1fp(i,j) = pow(g2,-1)*fp[i]*(-g1+2*last_obs[i]*hp[i])*d_mp1(i,jnd);
          d1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = spc1*fpc[i]*d_mp1(i,jnd);
        } else {
          int jnd = j-mp1_size-mp2_size;
          d1sp(i,j) = 0;
          d1fp(i,j) = 0;
          d1hp(i,j) = 0;
          d1spc(i,j) = -fpc[i]*d_cr(i,jnd);
        }
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix d1i(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      if (j < mp1_size+mp2_size) {
        d1i(i,j) = 0;
      } else {
        int jnd = j-mp1_size-mp2_size;
        d1i(i,j) = pow(approx_sigma,-1)*ind_approx[i]*(1-ind_approx[i])*d_cr(i,jnd);
      }
    }
  }

  NumericVector d1p(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    if (i >= mp1_size+mp2_size) {
      int jnd = i-mp1_size-mp2_size;
      d1p[i] = smooth_para*cr[jnd];
    } else {
      d1p[i] = 0;
    }
  }

  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      double th = general_hazard[i]+hp[i]*ind_approx[i];
      double ts = (sp[i]-spc[i])*ind_approx[i]+spc[i];
      double fe = fp[i];
      double tdh = pow(th,-1)*(hp[i]*d1i(i,j)+d1hp(i,j)*ind_approx[i]);
      double tds = pow(ts,-1)*((d1sp(i,j)-d1spc(i,j))*ind_approx[i]+(sp[i]-spc[i])*d1i(i,j)+d1spc(i,j));
      double tdf = tdh+tds;
      double dfe = pow(fe,-1)*d1fp(i,j);
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status3[i]*dfe+
            cen_status4[i]*tdf-d1p[j]);
    }
  }

  double sum_bi = sum(bi);
  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}

inline NumericVector surv::grid_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();
  int cure_range = y.size();
  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }
  NumericVector out(cure_range);
  for (int i = 0; i < cure_range; i++) {
    NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
    for (int ii = 0; ii < sam_size; ii++) {
      double p1 = 0, p2 = 0, pc = 0;
      for (int j = 0; j < mp1_size; j++) {
        p1 += mp1_covariate(ii,j)*mp1[j];
      }
      for (int j = 0; j < mp2_size; j++) {
        p2 += mp2_covariate(ii,j)*mp2[j];
      }
      for (int j = 0; j < cure_size; j++) {
        pc += cure_covariate(ii,j)*y[i];
      }
      m1[ii] = p1;
      m2[ii] = p2;
      m3[ii] = pc;
    }

    NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
    NumericVector spc(sam_size), fpc(sam_size);
    if (dist.compare("Weibull") == 0) {
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = exp(m1[ii]);
        double g2 = exp(m2[ii]);
        double gc = exp(m3[ii]);
        double sp1 = last_obs[ii]/g2;
        double sp2 = pow(sp1, g1);
        sp[ii] = exp(-sp2);
        double hp1 = g1/g2;
        double hp2 = pow(sp1, g1-1);
        hp[ii] = hp1*hp2;
        fp[ii] = sp[ii]*hp[ii];
        double spc1 = gc/g2;
        double spc2 = pow(spc1, g1);
        spc[ii] = exp(-spc2);
        fpc[ii] = hp1*pow(spc1, g1-1)*spc[ii];
      }
    } else if (dist.compare("log-normal") == 0) {
      NumericVector standard(sam_size), standardc(sam_size);
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = m1[ii];
        double g2 = m2[ii];
        double gc = exp(m3[ii]);
        standard[ii] = (log(last_obs[ii])-g1)/g2;
        standardc[ii] = (log(gc)-g1)/g2;
      }
      sp = 1-pnorm(standard);
      spc = 1-pnorm(standardc);
      NumericVector fn = dnorm(standard);
      NumericVector fnc = dnorm(standardc);
      for (int ii = 0; ii < sam_size; ii++) {
        double g2 = m2[ii];
        double gc = exp(m3[ii]);
        fp[ii] = fn[ii]/(g2*last_obs[ii]);
        fpc[ii] = fnc[ii]/(g2*gc);
        hp[ii] = fp[ii]/sp[ii];
      }
    } else if (dist.compare("log-logistic") == 0) {
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = exp(m1[ii]);
        double g2 = exp(m2[ii]);
        double gc = exp(m3[ii]);
        double sp1 = last_obs[ii]/g2;
        double sp2 = pow(sp1, g1);
        sp[ii] = pow(1+sp2, -1);
        double hp1 = g1/g2;
        double hp2 = pow(sp1, g1-1);
        hp[ii] = hp1*hp2*sp[ii];
        fp[ii] = sp[ii]*hp[ii];
        double spc1 = gc/g2;
        double spc2 = pow(spc1, g1);
        spc[ii] = pow(1+spc2, -1);
        double hpc2 = pow(spc1, g1-1);
        double hpc = hp1*hpc2*spc[ii];
        fpc[ii] = spc[ii]*hpc;
      }
    } else {
      stop("Weibull, log-normal, or log-logistic");
    }

    NumericVector ind(sam_size);
    for (int ii = 0; ii < sam_size; ii++) {
      double gc = exp(m3[ii]);
      if (gc - last_obs[ii] > 0) {
        ind[ii] = 1;
      } else {
        ind[ii] = 0;
      }
    }
    NumericVector individual_loglikelihood(sam_size);
    for (int ii = 0; ii < sam_size; ii++) {
      double th = general_hazard[ii]+hp[ii]*ind[ii];
      double ts = (sp[ii]-spc[ii])*ind[ii]+spc[ii];
      double tf = th*ts;
      double fe = fp[ii]*ind[ii];
      double loglike1 = log(pow(ts, cen_status1[ii]));
      double loglike2 = log(pow(ts, cen_status2[ii]));
      double loglike3;
      if (fe == 0) {
        loglike3 = 0;
      } else {
        loglike3 = log(pow(fe, cen_status3[ii]));
      }
      double loglike4 = log(pow(tf, cen_status4[ii]));
      individual_loglikelihood[ii] = bi[ii]*(loglike1+loglike2+loglike3+loglike4);
    }
    out[i] = sum(individual_loglikelihood)/sum_bi;
  }
  return out;
}

inline NumericVector surv::grid_mixture_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();
  int cure_range = y.size();

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }
  NumericVector out(cure_range);
  for (int i = 0; i < cure_range; i++) {
    NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
    for (int ii = 0; ii < sam_size; ii++) {
      double p1 = 0, p2 = 0, pc = 0;
      for (int j = 0; j < mp1_size; j++) {
        p1 += mp1_covariate(ii,j)*mp1[j];
      }
      for (int j = 0; j < mp2_size; j++) {
        p2 += mp2_covariate(ii,j)*mp2[j];
      }
      for (int j = 0; j < cure_size; j++) {
        pc += cure_covariate(ii,j)*y[i];
      }
      m1[ii] = p1;
      m2[ii] = p2;
      m3[ii] = pc;
    }

    NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
    if (dist.compare("Weibull") == 0) {
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = exp(m1[ii]);
        double g2 = exp(m2[ii]);
        double sp1 = last_obs[ii]/g2;
        double sp2 = pow(sp1, g1);
        sp[ii] = exp(-sp2);
        double hp1 = g1/g2;
        double hp2 = pow(sp1, g1-1);
        hp[ii] = hp1*hp2;
        fp[ii] = sp[ii]*hp[ii];
      }
    } else if (dist.compare("log-normal") == 0) {
      NumericVector standard(sam_size), standardc(sam_size);
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = m1[ii];
        double g2 = m2[ii];
        standard[ii] = (log(last_obs[ii])-g1)/g2;
      }
      sp = 1-pnorm(standard);
      NumericVector fn = dnorm(standard);
      for (int ii = 0; ii < sam_size; ii++) {
        double g2 = m2[ii];
        fp[ii] = fn[ii]/(g2*last_obs[ii]);
        hp[ii] = fp[ii]/sp[ii];
      }
    } else if (dist.compare("log-logistic") == 0) {
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = exp(m1[ii]);
        double g2 = exp(m2[ii]);
        double sp1 = last_obs[ii]/g2;
        double sp2 = pow(sp1, g1);
        sp[ii] = pow(1+sp2, -1);
        double hp1 = g1/g2;
        double hp2 = pow(sp1, g1-1);
        hp[ii] = hp1*hp2*sp[ii];
        fp[ii] = sp[ii]*hp[ii];
      }
    } else {
      stop("Weibull, log-normal, or log-logistic");
    }

    NumericVector individual_loglikelihood(sam_size);
    for (int ii = 0; ii < sam_size; ii++) {
      double pipi = pow(1+exp(-m3[ii]), -1);
      double se = pipi+(1-pipi)*sp[ii];
      double fe = (1-pipi)*fp[ii];
      double he = fe*pow(se, -1);
      double th = general_hazard[ii]+he;
      double ts = se;
      double tf = th*ts;
      double loglike1 = log(pow(ts, cen_status1[ii]));
      double loglike2 = log(pow(ts, cen_status2[ii]));
      double loglike3 = log(pow(fe, cen_status3[ii]));
      double loglike4 = log(pow(tf, cen_status4[ii]));
      individual_loglikelihood[ii] = bi[ii]*(loglike1+loglike2+loglike3+loglike4);
    }
    out[i] = sum(individual_loglikelihood)/sum_bi;
  }
  return out;
}

inline NumericVector surv::grid_nmixture_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();
  int cure_range = y.size();

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }
  NumericVector out(cure_range);
  for (int i = 0; i < cure_range; i++) {
    NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
    for (int ii = 0; ii < sam_size; ii++) {
      double p1 = 0, p2 = 0, pc = 0;
      for (int j = 0; j < mp1_size; j++) {
        p1 += mp1_covariate(ii,j)*mp1[j];
      }
      for (int j = 0; j < mp2_size; j++) {
        p2 += mp2_covariate(ii,j)*mp2[j];
      }
      for (int j = 0; j < cure_size; j++) {
        pc += cure_covariate(ii,j)*y[i];
      }
      m1[ii] = p1;
      m2[ii] = p2;
      m3[ii] = pc;
    }

    NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
    NumericVector spc(sam_size), fpc(sam_size);
    if (dist.compare("Weibull") == 0) {
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = exp(m1[ii]);
        double g2 = exp(m2[ii]);
        double sp1 = last_obs[ii]/g2;
        double sp2 = pow(sp1, g1);
        sp[ii] = exp(-sp2);
        double hp1 = g1/g2;
        double hp2 = pow(sp1, g1-1);
        hp[ii] = hp1*hp2;
        fp[ii] = sp[ii]*hp[ii];
      }
    } else if (dist.compare("log-normal") == 0) {
      NumericVector standard(sam_size), standardc(sam_size);
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = m1[ii];
        double g2 = m2[ii];
        standard[ii] = (log(last_obs[ii])-g1)/g2;
      }
      sp = 1-pnorm(standard);
      NumericVector fn = dnorm(standard);
      for (int ii = 0; ii < sam_size; ii++) {
        double g2 = m2[ii];
        fp[ii] = fn[ii]/(g2*last_obs[ii]);
        hp[ii] = fp[ii]/sp[ii];
      }
    } else if (dist.compare("log-logistic") == 0) {
      for (int ii = 0; ii < sam_size; ii++) {
        double g1 = exp(m1[ii]);
        double g2 = exp(m2[ii]);
        double sp1 = last_obs[ii]/g2;
        double sp2 = pow(sp1, g1);
        sp[ii] = pow(1+sp2, -1);
        double hp1 = g1/g2;
        double hp2 = pow(sp1, g1-1);
        hp[ii] = hp1*hp2*sp[ii];
        fp[ii] = sp[ii]*hp[ii];
      }
    } else {
      stop("Weibull, log-normal, or log-logistic");
    }

    NumericVector individual_loglikelihood(sam_size);
    for (int ii = 0; ii < sam_size; ii++) {
      double pipi = pow(1+exp(-m3[ii]), -1);
      double se = pow(pipi, 1-sp[ii]);
      double he = -fp[ii]*log(pipi);
      double fe = se*he;
      double th = general_hazard[ii]+he;
      double ts = se;
      double tf = th*ts;
      double loglike1 = log(pow(ts, cen_status1[ii]));
      double loglike2 = log(pow(ts, cen_status2[ii]));
      double loglike3 = log(pow(fe, cen_status3[ii]));
      double loglike4 = log(pow(tf, cen_status4[ii]));
      individual_loglikelihood[ii] = bi[ii]*(loglike1+loglike2+loglike3+loglike4);
    }
    out[i] = sum(individual_loglikelihood)/sum_bi;
  }
  return out;
}

inline double surv::ite_model_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cure[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  NumericVector spc(sam_size), fpc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
      fpc[i] = hp1*pow(spc1, g1-1)*spc[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double gc = exp(m3[i]);
      fp[i] = fn[i]/(g2*last_obs[i]);
      fpc[i] = fnc[i]/(g2*gc);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = pow(1+spc2, -1);
      double hpc2 = pow(spc1, g1-1);
      double hpc = hp1*hpc2*spc[i];
      fpc[i] = spc[i]*hpc;
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector ind(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    if (gc - last_obs[i] > 0) {
      ind[i] = 1;
    } else {
      ind[i] = 0;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double th = general_hazard[i]+hp[i]*ind[i];
    double ts = (sp[i]-spc[i])*ind[i]+spc[i];
    double tf = th*ts;
    double fe = fp[i];
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4);
  }

  double sum_bi = sum(bi);
  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline NumericVector surv::ite_model_score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cure[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_mp1(sam_size, mp1_size);
  NumericMatrix d_mp2(sam_size, mp2_size);
  if ((dist.compare("Weibull") == 0) || (dist.compare("log-logistic") == 0)) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = g1*mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = g2*mp2_covariate(i,j);
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    for (int i = 0; i < sam_size; i++) {
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = mp2_covariate(i,j);
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  NumericVector spc(sam_size), fpc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
      fpc[i] = hp1*pow(spc1, g1-1)*spc[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double gc = exp(m3[i]);
      fp[i] = fn[i]/(g2*last_obs[i]);
      fpc[i] = fnc[i]/(g2*gc);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = pow(1+spc2, -1);
      double hpc2 = pow(spc1, g1-1);
      double hpc = hp1*hpc2*spc[i];
      fpc[i] = hpc*spc[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector ind(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    if (gc - last_obs[i] > 0) {
      ind[i] = 1;
    } else {
      ind[i] = 0;
    }
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  NumericMatrix d1spc(sam_size, allpara_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1sp(i,j) = -sp2*log(sp1)*sp[i]*d_mp1(i,j);
          d1fp(i,j) = (pow(g1,-1)+log(sp1)-sp2*log(sp1))*fp[i]*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = -spc2*log(spc1)*spc[i]*d_mp1(i,j);
        } else {
          int jnd = j-mp1_size;
          d1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          d1fp(i,j) = (last_obs[i]*hp[i]-g1)*pow(g2,-1)*fp[i]*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = spc1*fpc[i]*d_mp2(i,jnd);
        }
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = exp(m2[i]);
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1fp(i,j) = (fp[i]*standard[i]/g2)*d_mp1(i,j);
          d1sp(i,j) = (fn[i]/g2)*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = (fnc[i]/g2)*d_mp1(i,j);
        } else {
          int jnd = j-mp1_size;
          d1fp(i,j) = (fp[i]/g2)*(pow(standard[i], 2)-1)*d_mp2(i,jnd);
          d1sp(i,j) = ((fn[i]*standard[i])/g2)*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = ((fnc[i]*standardc[i])/g2)*d_mp2(i,jnd);
        }
      }
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double spc1 = gc/g2;
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          d1fp(i,j) = fp[i]*(pow(g1,-1)+log(sp1)-2*last_obs[i]*pow(g1,-1)*log(sp1)*hp[i])*d_mp1(i,j);
          d1sp(i,j) = -last_obs[i]*pow(g1,-1)*log(sp1)*fp[i]*d_mp1(i,j);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = -gc*pow(g1,-1)*log(spc1)*fpc[i]*d_mp1(i,j);
        } else {
          int jnd = j-mp1_size;
          d1fp(i,j) = pow(g2,-1)*fp[i]*(-g1+2*last_obs[i]*hp[i])*d_mp1(i,jnd);
          d1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          d1hp(i,j) = (d1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*d1sp(i,j);
          d1spc(i,j) = spc1*fpc[i]*d_mp1(i,jnd);
        }
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      double th = general_hazard[i]+hp[i]*ind[i];
      double ts = (sp[i]-spc[i])*ind[i]+spc[i];
      double fe = fp[i];
      double tdh = pow(th,-1)*(d1hp(i,j)*ind[i]);
      double tds = pow(ts,-1)*((d1sp(i,j)-d1spc(i,j))*ind[i]+d1spc(i,j));
      double tdf = tdh+tds;
      double dfe = pow(fe,-1)*d1fp(i,j);
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status3[i]*dfe+
            cen_status4[i]*tdf);
    }
  }

  double sum_bi = sum(bi);
  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}

inline double surv::ite_curetime_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*y[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  NumericVector spc(sam_size), fpc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
      fpc[i] = hp1*pow(spc1, g1-1)*spc[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double gc = exp(m3[i]);
      fp[i] = fn[i]/(g2*last_obs[i]);
      fpc[i] = fnc[i]/(g2*gc);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = pow(1+spc2, -1);
      double hpc2 = pow(spc1, g1-1);
      double hpc = hp1*hpc2*spc[i];
      fpc[i] = spc[i]*hpc;
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double approx_sigma = pow(sam_size, -1/appro_weight);
  NumericVector ind_approx(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    ind_approx[i] = pow(1+exp(-(gc-last_obs[i])*pow(approx_sigma,-1)),-1);
  }

  double inner_cure = 0;
  for (int i = 0; i < cure_size; i++) {
    inner_cure += pow(y[i], 2);
  }
  double penalty = 0.5*smooth_para*inner_cure;

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double th = general_hazard[i]+hp[i]*ind_approx[i];
    double ts = (sp[i]-spc[i])*ind_approx[i]+spc[i];
    double tf = th*ts;
    double fe = fp[i];
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4-penalty);
  }

  double sum_bi = sum(bi);
  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline NumericVector surv::ite_curetime_score(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = y.size();

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*y[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_cr(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    for (int j = 0; j < cure_size; j++) {
      d_cr(i,j) = gc*cure_covariate(i,j);
    }
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  NumericVector spc(sam_size), fpc(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = exp(-spc2);
      fpc[i] = hp1*pow(spc1, g1-1)*spc[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      double gc = exp(m3[i]);
      standard[i] = (log(last_obs[i])-g1)/g2;
      standardc[i] = (log(gc)-g1)/g2;
    }
    sp = 1-pnorm(standard);
    spc = 1-pnorm(standardc);
    NumericVector fn = dnorm(standard);
    NumericVector fnc = dnorm(standardc);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      double gc = exp(m3[i]);
      fp[i] = fn[i]/(g2*last_obs[i]);
      fpc[i] = fnc[i]/(g2*gc);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc = exp(m3[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
      double spc1 = gc/g2;
      double spc2 = pow(spc1, g1);
      spc[i] = pow(1+spc2, -1);
      double hpc2 = pow(spc1, g1-1);
      double hpc = hp1*hpc2*spc[i];
      fpc[i] = spc[i]*hpc;
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double approx_sigma = pow(sam_size, -1/appro_weight);
  NumericVector ind_approx(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double gc = exp(m3[i]);
    ind_approx[i] = pow(1+exp(-(gc-last_obs[i])*pow(approx_sigma,-1)), -1);
  }

  NumericMatrix d1i(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < cure_size; j++) {
      d1i(i,j) = pow(approx_sigma,-1)*ind_approx[i]*(1-ind_approx[i])*d_cr(i,j);
    }
  }

  NumericVector d1p(cure_size);
  for (int i = 0; i < cure_size; i++) {
    d1p[i] = smooth_para*y[i];
  }

  NumericMatrix d1spc(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < cure_size; j++) {
      d1spc(i,j) = -fpc[i]*d_cr(i,j);
    }
  }

  NumericMatrix individual_score(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < cure_size; j++) {
      double th = general_hazard[i]+hp[i]*ind_approx[i];
      double ts = (sp[i]-spc[i])*ind_approx[i]+spc[i];
      double tdh = pow(th,-1)*hp[i]*d1i(i,j);
      double tds = pow(ts,-1)*((0-d1spc(i,j))*ind_approx[i]+(sp[i]-spc[i])*d1i(i,j)+d1spc(i,j));
      double tdf = tdh+tds;
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status4[i]*tdf-d1p[j]);
    }
  }

  double sum_bi = sum(bi);
  NumericVector out(cure_size);
  for (int i = 0; i < cure_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}

inline double surv::mixture_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size), cr(cure_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else if ((i >= mp1_size) && (i < mp1_size+mp2_size)) {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    } else {
      int jnd = i-mp1_size-mp2_size;
      cr[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cr[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double inner_cure = 0;
  for (int i = 0; i < cure_size; i++) {
    inner_cure += pow(cr[i], 2);
  }
  double penalty = 0.5*smooth_para*inner_cure;

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double pipi = pow(1+exp(-m3[i]), -1);
    double se = pipi+(1-pipi)*sp[i];
    double fe = (1-pipi)*fp[i];
    double he = fe*pow(se, -1);
    double th = general_hazard[i]+he;
    double ts = se;
    double tf = th*ts;
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4-penalty);
  }

  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline NumericVector surv::mixture_score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int model_size = mp1_size+mp2_size;
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size), cr(cure_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else if ((i >= mp1_size) && (i < model_size)) {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    } else {
      int jnd = i-model_size;
      cr[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cr[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_mp1(sam_size, mp1_size);
  NumericMatrix d_mp2(sam_size, mp2_size);
  NumericMatrix d_cr(sam_size, cure_size);
  if ((dist.compare("Weibull") == 0) || (dist.compare("log-logistic") == 0)) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = g1*mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = g2*mp2_covariate(i,j);
      }
      for (int j = 0; j < cure_size; j++) {
        d_cr(i,j) = gc*(1-gc)*cure_covariate(i,j);
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = mp2_covariate(i,j);
      }
      for (int j = 0; j < cure_size; j++) {
        d_cr(i,j) = gc*(1-gc)*cure_covariate(i,j);
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix md1sp(sam_size, model_size);
  NumericMatrix md1hp(sam_size, model_size);
  NumericMatrix md1fp(sam_size, model_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      for (int j = 0; j < model_size; j++) {
        if (j < mp1_size) {
          md1sp(i,j) = -sp2*log(sp1)*sp[i]*d_mp1(i,j);
          md1fp(i,j) = (pow(g1,-1)+log(sp1)-sp2*log(sp1))*fp[i]*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          md1fp(i,j) = (last_obs[i]*hp[i]-g1)*pow(g2,-1)*fp[i]*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      for (int j = 0; j < model_size; j++) {
        if (j < mp1_size) {
          md1fp(i,j) = (fp[i]*standard[i]/g2)*d_mp1(i,j);
          md1sp(i,j) = (fn[i]/g2)*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1fp(i,j) = (fp[i]/g2)*(pow(standard[i], 2)-1)*d_mp2(i,jnd);
          md1sp(i,j) = ((fn[i]*standard[i])/g2)*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      for (int j = 0; j < model_size; j++) {
        if (j < mp1_size) {
          md1fp(i,j) = fp[i]*(pow(g1,-1)+log(sp1)-2*last_obs[i]*pow(g1,-1)*log(sp1)*hp[i])*d_mp1(i,j);
          md1sp(i,j) = -last_obs[i]*pow(g1,-1)*log(sp1)*fp[i]*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1fp(i,j) = pow(g2,-1)*fp[i]*(-g1+2*last_obs[i]*hp[i])*d_mp1(i,jnd);
          md1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    double gc_inv = 1+exp(-m3[i]);
    double gc = pow(gc_inv, -1);
    double se = gc+(1-gc)*sp[i];
    double fe = (1-gc)*fp[i];
    double se_inv = pow(se, -1);
    double he = fe*se_inv;
    for (int j = 0; j < allpara_size; j++) {
      if (j < model_size) {
        d1sp(i,j) = (1-gc)*md1sp(i,j);
        d1hp(i,j) = (1-gc)*md1fp(i,j)*se_inv-pow(1-gc, 2)*fp[i]*md1sp(i,j)*pow(se_inv, 2);
        d1fp(i,j) = he*d1sp(i,j)+d1hp(i,j)*se;
      } else {
        int jnd = j-model_size;
        d1sp(i,j) = (1-sp[i])*d_cr(i,jnd);
        d1hp(i,j) = -fp[i]*se_inv*(1-(1-gc)*(1-sp[i])*se_inv)*d_cr(i,jnd);
        d1fp(i,j) = he*d1sp(i,j)+d1hp(i,j)*se;
      }
    }
  }

  NumericVector d1p(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    if (i >= model_size) {
      int jnd = i-model_size;
      d1p[i] = smooth_para*cr[jnd];
    } else {
      d1p[i] = 0;
    }
  }

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      double se = gc+(1-gc)*sp[i];
      double fe = (1-gc)*fp[i];
      double dfe;
      if (fe == 0) {
        dfe = 0;
      } else {
        dfe = pow(fe,-1)*d1fp(i,j);
      };
      double he = fe*pow(se, -1);
      double th = general_hazard[i]+he;
      double ts = se;
      double tdh = pow(th,-1)*d1hp(i,j);
      double tds = pow(ts,-1)*d1sp(i,j);
      double tdf = tdh+tds;
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status3[i]*dfe+
            cen_status4[i]*tdf-d1p[j]);
    }
  }

  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}

inline double surv::nmixture_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size), cr(cure_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else if ((i >= mp1_size) && (i < mp1_size+mp2_size)) {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    } else {
      int jnd = i-mp1_size-mp2_size;
      cr[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cr[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double inner_cure = 0;
  for (int i = 0; i < cure_size; i++) {
    inner_cure += pow(cr[i], 2);
  }
  double penalty = 0.5*smooth_para*inner_cure;

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double pipi = pow(1+exp(-m3[i]), -1);
    double se = pow(pipi, 1-sp[i]);
    double he = -fp[i]*log(pipi);
    double fe = se*he;
    double th = general_hazard[i]+he;
    double ts = se;
    double tf = ts*th;
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4-penalty);
  }

  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline NumericVector surv::nmixture_score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int model_size = mp1_size+mp2_size;
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size), cr(cure_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else if ((i >= mp1_size) && (i < model_size)) {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    } else {
      int jnd = i-model_size;
      cr[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cr[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_mp1(sam_size, mp1_size);
  NumericMatrix d_mp2(sam_size, mp2_size);
  NumericMatrix d_cr(sam_size, cure_size);
  if ((dist.compare("Weibull") == 0) || (dist.compare("log-logistic") == 0)) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = g1*mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = g2*mp2_covariate(i,j);
      }
      for (int j = 0; j < cure_size; j++) {
        d_cr(i,j) = gc*(1-gc)*cure_covariate(i,j);
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = mp2_covariate(i,j);
      }
      for (int j = 0; j < cure_size; j++) {
        d_cr(i,j) = gc*(1-gc)*cure_covariate(i,j);
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix md1sp(sam_size, model_size);
  NumericMatrix md1hp(sam_size, model_size);
  NumericMatrix md1fp(sam_size, model_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      for (int j = 0; j < model_size; j++) {
        if (j < mp1_size) {
          md1sp(i,j) = -sp2*log(sp1)*sp[i]*d_mp1(i,j);
          md1fp(i,j) = (pow(g1,-1)+log(sp1)-sp2*log(sp1))*fp[i]*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          md1fp(i,j) = (last_obs[i]*hp[i]-g1)*pow(g2,-1)*fp[i]*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      for (int j = 0; j < model_size; j++) {
        if (j < mp1_size) {
          md1fp(i,j) = (fp[i]*standard[i]/g2)*d_mp1(i,j);
          md1sp(i,j) = (fn[i]/g2)*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1fp(i,j) = (fp[i]/g2)*(pow(standard[i], 2)-1)*d_mp2(i,jnd);
          md1sp(i,j) = ((fn[i]*standard[i])/g2)*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      for (int j = 0; j < model_size; j++) {
        if (j < mp1_size) {
          md1fp(i,j) = fp[i]*(pow(g1,-1)+log(sp1)-2*last_obs[i]*pow(g1,-1)*log(sp1)*hp[i])*d_mp1(i,j);
          md1sp(i,j) = -last_obs[i]*pow(g1,-1)*log(sp1)*fp[i]*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1fp(i,j) = pow(g2,-1)*fp[i]*(-g1+2*last_obs[i]*hp[i])*d_mp1(i,jnd);
          md1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    double gc_inv = 1+exp(-m3[i]);
    double gc = pow(gc_inv, -1);
    double he = -fp[i]*log(gc);
    double se = pow(gc, 1-sp[i]);
    for (int j = 0; j < allpara_size; j++) {
      if (j < model_size) {
        d1sp(i,j) = -se*log(gc)*md1sp(i,j);
        d1hp(i,j) = -log(gc)*md1fp(i,j);
        d1fp(i,j) = he*d1sp(i,j)+d1hp(i,j)*se;
      } else {
        int jnd = j-model_size;
        d1sp(i,j) = (1-sp[i])*gc_inv*se*d_cr(i,jnd);
        d1hp(i,j) = -fp[i]*gc_inv*d_cr(i,jnd);
        d1fp(i,j) = he*d1sp(i,j)+d1hp(i,j)*se;
      }
    }
  }

  NumericVector d1p(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    if (i >= model_size) {
      int jnd = i-model_size;
      d1p[i] = smooth_para*cr[jnd];
    } else {
      d1p[i] = 0;
    }
  }

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      double se = pow(gc, 1-sp[i]);
      double he = -fp[i]*log(gc);
      double fe = se*he;
      double dfe;
      if (fe == 0) {
        dfe = 0;
      } else {
        dfe = pow(fe,-1)*d1fp(i,j);
      }
      double th = general_hazard[i]+he;
      double ts = se;
      double tdh = pow(th,-1)*d1hp(i,j);
      double tds = pow(ts,-1)*d1sp(i,j);
      double tdf = tdh+tds;
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status3[i]*dfe+
            cen_status4[i]*tdf-d1p[j]);
    }
  }

  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}

inline double surv::ite_model_mixture_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cure[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double pipi = pow(1+exp(-m3[i]), -1);
    double se = pipi+(1-pipi)*sp[i];
    double fe = (1-pipi)*fp[i];
    double he = fe*pow(se, -1);
    double th = general_hazard[i]+he;
    double ts = se;
    double tf = th*ts;
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4);
  }

  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline NumericVector surv::ite_model_mixture_score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cure[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_mp1(sam_size, mp1_size);
  NumericMatrix d_mp2(sam_size, mp2_size);
  if ((dist.compare("Weibull") == 0) || (dist.compare("log-logistic") == 0)) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = g1*mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = g2*mp2_covariate(i,j);
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    for (int i = 0; i < sam_size; i++) {
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = mp2_covariate(i,j);
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix md1sp(sam_size, allpara_size);
  NumericMatrix md1hp(sam_size, allpara_size);
  NumericMatrix md1fp(sam_size, allpara_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          md1sp(i,j) = -sp2*log(sp1)*sp[i]*d_mp1(i,j);
          md1fp(i,j) = (pow(g1,-1)+log(sp1)-sp2*log(sp1))*fp[i]*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          md1fp(i,j) = (last_obs[i]*hp[i]-g1)*pow(g2,-1)*fp[i]*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = exp(m2[i]);
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          md1fp(i,j) = (fp[i]*standard[i]/g2)*d_mp1(i,j);
          md1sp(i,j) = (fn[i]/g2)*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1fp(i,j) = (fp[i]/g2)*(pow(standard[i], 2)-1)*d_mp2(i,jnd);
          md1sp(i,j) = ((fn[i]*standard[i])/g2)*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          md1fp(i,j) = fp[i]*(pow(g1,-1)+log(sp1)-2*last_obs[i]*pow(g1,-1)*log(sp1)*hp[i])*d_mp1(i,j);
          md1sp(i,j) = -last_obs[i]*pow(g1,-1)*log(sp1)*fp[i]*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1fp(i,j) = pow(g2,-1)*fp[i]*(-g1+2*last_obs[i]*hp[i])*d_mp1(i,jnd);
          md1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    double gc_inv = 1+exp(-m3[i]);
    double gc = pow(gc_inv, -1);
    double se = gc+(1-gc)*sp[i];
    double se_inv = pow(se, -1);
    double he = (1-gc)*fp[i]*se_inv;
    for (int j = 0; j < allpara_size; j++) {
      d1sp(i,j) = (1-gc)*md1sp(i,j);
      d1hp(i,j) = (1-gc)*md1fp(i,j)*se_inv-pow(1-gc, 2)*fp[i]*md1sp(i,j)*pow(se_inv, 2);
      d1fp(i,j) = he*d1sp(i,j)+d1hp(i,j)*se;
    }
  }

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      double se = gc+(1-gc)*sp[i];
      double fe = (1-gc)*fp[i];
      double dfe;
      if (fe == 0) {
        dfe = 0;
      } else {
        dfe = pow(fe,-1)*d1fp(i,j);
      }
      double he = fe*pow(se, -1);
      double th = general_hazard[i]+he;
      double ts = se;
      double tdh = pow(th,-1)*d1hp(i,j);
      double tds = pow(ts,-1)*d1sp(i,j);
      double tdf = tdh+tds;
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status3[i]*dfe+
            cen_status4[i]*tdf);
    }
  }

  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}

inline double surv::ite_model_nmixture_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cure[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double pipi = pow(1+exp(-m3[i]), -1);
    double se = pow(pipi, 1-sp[i]);
    double he = -fp[i]*log(pipi);
    double fe = se*he;
    double th = general_hazard[i]+he;
    double ts = se;
    double tf = th*ts;
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4);
  }

  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline NumericVector surv::ite_model_nmixture_score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector mp1(mp1_size), mp2(mp2_size);
  for(int i = 0; i < allpara_size; i++) {
    if (i < mp1_size) {
      mp1[i] = y[i];
    } else {
      int jnd = i-mp1_size;
      mp2[jnd] = y[i];
    }
  }

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*cure[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_mp1(sam_size, mp1_size);
  NumericMatrix d_mp2(sam_size, mp2_size);
  if ((dist.compare("Weibull") == 0) || (dist.compare("log-logistic") == 0)) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = g1*mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = g2*mp2_covariate(i,j);
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    for (int i = 0; i < sam_size; i++) {
      for (int j = 0; j < mp1_size; j++) {
        d_mp1(i,j) = mp1_covariate(i,j);
      }
      for (int j = 0; j < mp2_size; j++) {
        d_mp2(i,j) = mp2_covariate(i,j);
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix md1sp(sam_size, allpara_size);
  NumericMatrix md1hp(sam_size, allpara_size);
  NumericMatrix md1fp(sam_size, allpara_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          md1sp(i,j) = -sp2*log(sp1)*sp[i]*d_mp1(i,j);
          md1fp(i,j) = (pow(g1,-1)+log(sp1)-sp2*log(sp1))*fp[i]*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          md1fp(i,j) = (last_obs[i]*hp[i]-g1)*pow(g2,-1)*fp[i]*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = exp(m2[i]);
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          md1fp(i,j) = (fp[i]*standard[i]/g2)*d_mp1(i,j);
          md1sp(i,j) = (fn[i]/g2)*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1fp(i,j) = (fp[i]/g2)*(pow(standard[i], 2)-1)*d_mp2(i,jnd);
          md1sp(i,j) = ((fn[i]*standard[i])/g2)*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      for (int j = 0; j < allpara_size; j++) {
        if (j < mp1_size) {
          md1fp(i,j) = fp[i]*(pow(g1,-1)+log(sp1)-2*last_obs[i]*pow(g1,-1)*log(sp1)*hp[i])*d_mp1(i,j);
          md1sp(i,j) = -last_obs[i]*pow(g1,-1)*log(sp1)*fp[i]*d_mp1(i,j);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        } else {
          int jnd = j-mp1_size;
          md1fp(i,j) = pow(g2,-1)*fp[i]*(-g1+2*last_obs[i]*hp[i])*d_mp1(i,jnd);
          md1sp(i,j) = sp1*fp[i]*d_mp2(i,jnd);
          md1hp(i,j) = (md1fp(i,j)/sp[i])-(fp[i]/pow(sp[i], 2))*md1sp(i,j);
        }
      }
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    double gc_inv = 1+exp(-m3[i]);
    double gc = pow(gc_inv, -1);
    double se = pow(gc, 1-sp[i]);
    double he = -log(gc)*fp[i];
    for (int j = 0; j < allpara_size; j++) {
      d1sp(i,j) = -se*log(gc)*md1sp(i,j);
      d1hp(i,j) = -log(gc)*md1fp(i,j);
      d1fp(i,j) = he*d1sp(i,j)+d1hp(i,j)*se;
    }
  }

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      double se = pow(gc, 1-sp[i]);
      double he = -fp[i]*log(gc);
      double fe = se*he;
      double dfe;
      if (fe == 0) {
        dfe = 0;
      } else {
        dfe = pow(fe,-1)*d1fp(i,j);
      }
      double th = general_hazard[i]+he;
      double ts = se;
      double tdh = pow(th,-1)*d1hp(i,j);
      double tds = pow(ts,-1)*d1sp(i,j);
      double tdf = tdh+tds;
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status3[i]*dfe+
            cen_status4[i]*tdf);
    }
  }

  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}

inline double surv::ite_curerate_mixture_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*y[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double inner_cure = 0;
  for (int i = 0; i < cure_size; i++) {
    inner_cure += pow(y[i], 2);
  }
  double penalty = 0.5*smooth_para*inner_cure;

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double pipi = pow(1+exp(-m3[i]), -1);
    double se = pipi+(1-pipi)*sp[i];
    double fe = (1-pipi)*fp[i];
    double he = fe*pow(se, -1);
    double th = general_hazard[i]+he;
    double ts = se;
    double tf = ts*th;
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4-penalty);
  }

  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline double surv::ite_curerate_nmixture_llh(NumericVector y) {
  int sam_size = last_obs.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*y[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  double inner_cure = 0;
  for (int i = 0; i < cure_size; i++) {
    inner_cure += pow(y[i], 2);
  }
  double penalty = 0.5*smooth_para*inner_cure;

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericVector individual_loglikelihood(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double pipi = pow(1+exp(-m3[i]), -1);
    double se = pow(pipi, 1-sp[i]);
    double he = -fp[i]*log(pipi);
    double fe = se*he;
    double th = general_hazard[i]+he;
    double ts = se;
    double tf = th*ts;
    double loglike1 = log(pow(ts, cen_status1[i]));
    double loglike2 = log(pow(ts, cen_status2[i]));
    double loglike3 = log(pow(fe, cen_status3[i]));
    double loglike4 = log(pow(tf, cen_status4[i]));
    individual_loglikelihood[i] = bi[i]*(loglike1+loglike2+loglike3+loglike4-penalty);
  }

  double out = sum(individual_loglikelihood)/sum_bi;
  return out;
}

inline NumericVector surv::ite_curerate_mixture_score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*y[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_cr(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    double gc_inv = 1+exp(-m3[i]);
    double gc = pow(gc_inv, -1);
    for (int j = 0; j < cure_size; j++) {
      d_cr(i,j) = gc*(1-gc)*cure_covariate(i,j);
    }
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    double gc_inv = 1+exp(-m3[i]);
    double gc = pow(gc_inv, -1);
    double se = gc+(1-gc)*sp[i];
    double se_inv = pow(se, -1);
    double he = (1-gc)*fp[i]*se_inv;
    for (int j = 0; j < allpara_size; j++) {
      d1sp(i,j) = (1-sp[i])*d_cr(i,j);
      d1hp(i,j) = -fp[i]*se_inv*(1-(1-gc)*(1-sp[i])*se_inv)*d_cr(i,j);
      d1fp(i,j) = he*d1sp(i,j)+d1hp(i,j)*se;
    }
  }

  NumericVector d1p(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    d1p[i] = smooth_para*y[i];
  }

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      double se = gc+(1-gc)*sp[i];
      double fe = (1-gc)*fp[i];
      double dfe;
      if (fe == 0) {
        dfe = 0;
      } else {
        dfe = pow(fe,-1)*d1fp(i,j);
      }
      double he = fe*pow(se, -1);
      double th = general_hazard[i]+he;
      double ts = se;
      double tdh = pow(th,-1)*d1hp(i,j);
      double tds = pow(ts,-1)*d1sp(i,j);
      double tdf = tdh+tds;
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status3[i]*dfe+
            cen_status4[i]*tdf-d1p[j]);
    }
  }

  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}

inline NumericVector surv::ite_curerate_nmixture_score(NumericVector y) {
  int sam_size = last_obs.size();
  int allpara_size = y.size();
  int mp1_size = mp1_covariate.ncol();
  int mp2_size = mp2_covariate.ncol();
  int cure_size = cure_covariate.ncol();

  NumericVector m1(sam_size), m2(sam_size), m3(sam_size);
  for (int i = 0; i < sam_size; i++) {
    double p1 = 0, p2 = 0, pc = 0;
    for (int j = 0; j < mp1_size; j++) {
      p1 += mp1_covariate(i,j)*mp1[j];
    }
    for (int j = 0; j < mp2_size; j++) {
      p2 += mp2_covariate(i,j)*mp2[j];
    }
    for (int j = 0; j < cure_size; j++) {
      pc += cure_covariate(i,j)*y[j];
    }
    m1[i] = p1;
    m2[i] = p2;
    m3[i] = pc;
  }

  NumericMatrix d_cr(sam_size, cure_size);
  for (int i = 0; i < sam_size; i++) {
    double gc_inv = 1+exp(-m3[i]);
    double gc = pow(gc_inv, -1);
    for (int j = 0; j < cure_size; j++) {
      d_cr(i,j) = gc*(1-gc)*cure_covariate(i,j);
    }
  }

  NumericVector sp(sam_size), hp(sam_size), fp(sam_size);
  if (dist.compare("Weibull") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = exp(-sp2);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2;
      fp[i] = sp[i]*hp[i];
    }
  } else if (dist.compare("log-normal") == 0) {
    NumericVector standard(sam_size), standardc(sam_size);
    for (int i = 0; i < sam_size; i++) {
      double g1 = m1[i];
      double g2 = m2[i];
      standard[i] = (log(last_obs[i])-g1)/g2;
    }
    sp = 1-pnorm(standard);
    NumericVector fn = dnorm(standard);
    for (int i = 0; i < sam_size; i++) {
      double g2 = m2[i];
      fp[i] = fn[i]/(g2*last_obs[i]);
      hp[i] = fp[i]/sp[i];
    }
  } else if (dist.compare("log-logistic") == 0) {
    for (int i = 0; i < sam_size; i++) {
      double g1 = exp(m1[i]);
      double g2 = exp(m2[i]);
      double sp1 = last_obs[i]/g2;
      double sp2 = pow(sp1, g1);
      sp[i] = pow(1+sp2, -1);
      double hp1 = g1/g2;
      double hp2 = pow(sp1, g1-1);
      hp[i] = hp1*hp2*sp[i];
      fp[i] = sp[i]*hp[i];
    }
  } else {
    stop("Weibull, log-normal, or log-logistic");
  }

  NumericMatrix d1sp(sam_size, allpara_size);
  NumericMatrix d1hp(sam_size, allpara_size);
  NumericMatrix d1fp(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    double gc_inv = 1+exp(-m3[i]);
    double gc = pow(gc_inv, -1);
    double se = gc+(1-gc)*sp[i];
    double he = -log(gc)*fp[i];
    for (int j = 0; j < allpara_size; j++) {
      d1sp(i,j) = (1-sp[i])*gc_inv*se*d_cr(i,j);
      d1hp(i,j) = -fp[i]*gc_inv*d_cr(i,j);
      d1fp(i,j) = he*d1sp(i,j)+d1hp(i,j)*se;
    }
  }

  NumericVector d1p(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    d1p[i] = smooth_para*y[i];
  }

  double sum_bi = sum(bi);
  if (!bl) {
    NumericVector bi(sam_size);
    for (int i = 0; i < sam_size; i++) {
      bi[i] = 1;
    }
  }

  NumericMatrix individual_score(sam_size, allpara_size);
  for (int i = 0; i < sam_size; i++) {
    for (int j = 0; j < allpara_size; j++) {
      double gc_inv = 1+exp(-m3[i]);
      double gc = pow(gc_inv, -1);
      double se = pow(gc, 1-sp[i]);
      double he = -fp[i]*log(gc);
      double fe = se*he;
      double dfe;
      if (fe == 0) {
        dfe = 0;
      } else {
        dfe = pow(fe,-1)*d1fp(i,j);
      }
      double th = general_hazard[i]+he;
      double ts = se;
      double tdh = pow(th,-1)*d1hp(i,j);
      double tds = pow(ts,-1)*d1sp(i,j);
      double tdf = tdh+tds;
      individual_score(i,j) =
        bi[i]*(
            cen_status1[i]*tds+
            cen_status2[i]*tds+
            cen_status3[i]*dfe+
            cen_status4[i]*tdf-d1p[j]);
    }
  }

  NumericVector out(allpara_size);
  for (int i = 0; i < allpara_size; i++) {
    out[i] = sum(individual_score(_,i))/sum_bi;
  }

  return out;
}








inline void forwarder_input_model(void* context, NumericVector y) {
  return static_cast<surv*>(context)->input_model(y);
}

inline void forwarder_input_cure(void* context, NumericVector y) {
  return static_cast<surv*>(context)->input_cure(y);
}

inline double forwarder_model_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->model_llh(y);
}

inline NumericVector forwarder_model_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->model_score(y);
}

inline double forwarder_curetime_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->curetime_llh(y);
}

inline NumericVector forwarder_curetime_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->curetime_score(y);
}

inline double forwarder_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->llh(y);
}

inline double forwarder_llh_indicator(void* context, NumericVector y) {
  return static_cast<surv*>(context)->llh_indicator(y);
}

inline NumericVector forwarder_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->score(y);
}

inline double forwarder_ite_model_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_model_llh(y);
}

inline NumericVector forwarder_ite_model_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_model_score(y);
}

inline double forwarder_ite_curetime_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_curetime_llh(y);
}

inline NumericVector forwarder_ite_curetime_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_curetime_score(y);
}

inline NumericVector forwarder_grid_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->grid_llh(y);
}

inline NumericVector forwarder_grid_mixture_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->grid_mixture_llh(y);
}

inline NumericVector forwarder_grid_nmixture_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->grid_nmixture_llh(y);
}

inline double forwarder_mixture_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->mixture_llh(y);
}

inline NumericVector forwarder_mixture_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->mixture_score(y);
}

inline double forwarder_nmixture_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->nmixture_llh(y);
}

inline NumericVector forwarder_nmixture_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->nmixture_score(y);
}

inline double forwarder_ite_model_mixture_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_model_mixture_llh(y);
}

inline double forwarder_ite_model_nmixture_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_model_nmixture_llh(y);
}

inline double forwarder_ite_curerate_mixture_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_curerate_mixture_llh(y);
}

inline double forwarder_ite_curerate_nmixture_llh(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_curerate_nmixture_llh(y);
}

inline NumericVector forwarder_ite_model_mixture_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_model_mixture_score(y);
}

inline NumericVector forwarder_ite_model_nmixture_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_model_nmixture_score(y);
}

inline NumericVector forwarder_ite_curerate_mixture_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_curerate_mixture_score(y);
}

inline NumericVector forwarder_ite_curerate_nmixture_score(void* context, NumericVector y) {
  return static_cast<surv*>(context)->ite_curerate_nmixture_score(y);
}


#endif
