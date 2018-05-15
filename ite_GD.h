#ifndef ITEGD_H
#define ITEGD_H

#include "fun.h"
#include "iteration.h"
#include "iteration_curetime.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

List ite_GD(std::string mtype,
            NumericMatrix y0,
            NumericMatrix y1,
            NumericMatrix y2,
            arma::mat X,
            arma::vec T,
            NumericVector ce,
            void* class_address,
            double tole,
            int stop_count);

#endif
