#ifndef GD_H
#define GD_H

#include "fun.h"
#include "iteration.h"
#include <Rcpp.h>

using namespace Rcpp;

List GD(std::string mtype,
        NumericMatrix y0,
        NumericMatrix y1,
        NumericMatrix y2,
        void* class_address,
        double tole,
        int stop_count);

#endif
