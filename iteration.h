#ifndef ITERATION_H
#define ITERATION_H

#include <Rcpp.h>

using namespace Rcpp;

List iteration(NumericVector y0,
               double (*like)(void*, NumericVector),
               NumericVector (*score_fun)(void*, NumericVector),
               void* context,
               double tole,
               int stop_count);

#endif
