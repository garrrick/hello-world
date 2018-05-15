#ifndef GS_H
#define GS_H

#include <RcppArmadillo.h>

using namespace Rcpp;

List GS(double GS_from,
        double GS_to,
        NumericVector (*gl)(void*, NumericVector),
        void* context,
        int stop_count);

#endif
