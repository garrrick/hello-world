#ifndef CEN_EST_H
#define CEN_EST_H

#include "fun.h"
#include "iteration.h"
#include <Rcpp.h>

using namespace Rcpp;

List cen_est(std::string mtype,
             NumericMatrix mi1,
             NumericMatrix mi2,
             void* class_address,
             double tole,
             int stop_count);

#endif
