#ifndef PARAMETRIC_BOOTSTRAP_SAMPLE_H
#define PARAMETRIC_BOOTSTRAP_SAMPLE_H

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
                                 double c4_of_c234);

#endif
