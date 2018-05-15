#include "GS.h"
#include <Rcpp.h>

using namespace Rcpp;

List GS(double GS_from,
        double GS_to,
        NumericVector (*gl)(void*, NumericVector),
        void* context,
        int stop_count) {
  if (stop_count > 10) {
    stop_count = 10;
  }
  List out;
  double cg_lo = GS_from;
  double cg_up = GS_to;
  double sep = pow(10,-1);
  int space = 1+ceil((cg_up-cg_lo)*pow(sep,-1));
  NumericVector cure_grid(space);
  for (int i = 0; i < space; i++) {
    cure_grid[i] = cg_lo+i*sep;
  }
  int count = 1;
  NumericVector cul = gl(context, cure_grid);
  double mcul = max(cul);
  int mindex = 0;
  while (cul[mindex] < mcul) {
    mindex++;
  }
  double x0 = cure_grid[mindex];
  double x1 = 0;
  while(TRUE){
    count++;
    cg_lo = x0-sep;
    cg_up = x0+sep;
    sep *= .1;
    space = 1+ceil((cg_up-cg_lo)*pow(sep,-1));
    NumericVector cure_grid(space);
    for (int i = 0; i < space; i++) {
      cure_grid[i] = cg_lo+i*sep;
    }
    cul = gl(context, cure_grid);
    mcul = max(cul);
    mindex = 0;
    while (cul[mindex] < mcul) {
      mindex++;
    }
    x1 = cure_grid[mindex];
    bool check_c = (count == stop_count);
    if (check_c) {
      out["est"] = x1;
      out["count"] = count;
      out["loglikelihood"] = mcul;
      break;
    } else {
      x0 = x1;
    }
  }
  return out;
}
