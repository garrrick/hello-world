#ifndef INITIAL_H
#define INITIAL_H

#include "fun.h"
#include "iteration.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

List initial(std::string mtype,
             NumericMatrix mi1,
             NumericMatrix mi2,
             NumericMatrix ci,
             arma::mat X,
             arma::vec T,
             NumericVector ce,
             void* class_address,
             double tole,
             int stop_count);

#endif
