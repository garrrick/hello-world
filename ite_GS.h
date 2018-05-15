#ifndef ITE_GS_H
#define ITE_GS_H

#include "fun.h"
#include "iteration.h"
#include "GS.h"
#include <Rcpp.h>

using namespace Rcpp;

List ite_GS(std::string mtype,
            NumericMatrix par1,
            NumericMatrix par2,
            double cu_from,
            double cu_to,
            void* class_address,
            double tole,
            int stop_count);

#endif
