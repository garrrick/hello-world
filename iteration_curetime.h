#ifndef ITERATION_CURETIME_H
#define ITERATION_CURETIME_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

List iteration_curetime(arma::vec p0,
                        double (*like)(void*, NumericVector),
                        NumericVector (*score_fun)(void*, NumericVector),
                        arma::mat Xe,
                        arma::vec Te,
                        void* context,
                        double tole,
                        int stop_count);

#endif
