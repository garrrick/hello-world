#include "iteration.h"
#include <Rcpp.h>

using namespace Rcpp;

List iteration(NumericVector y0,
               double (*like)(void*, NumericVector),
               NumericVector (*score_fun)(void*, NumericVector),
               void* context,
               double tole,
               int stop_count) {
  List est_out;
  LogicalVector check_y0_na = !is_na(y0);
  LogicalVector check_y0_nan = !is_nan(y0);
  LogicalVector check_y0_inf = !is_infinite(y0);
  LogicalVector check_y0 = check_y0_na & check_y0_nan & check_y0_inf;
  NumericVector l_value;
  l_value[0] = like(context, y0);
  LogicalVector check_y0_l_na = !is_na(l_value);
  LogicalVector check_y0_l_nan = !is_nan(l_value);
  LogicalVector check_y0_l_inf = !is_infinite(l_value);
  LogicalVector check_y0_l = check_y0_l_na & check_y0_l_nan & check_y0_l_inf;
  NumericVector s_values = score_fun(context, y0);
  LogicalVector check_y0_s_na = !is_na(s_values);
  LogicalVector check_y0_s_nan = !is_nan(s_values);
  LogicalVector check_y0_s_inf = !is_infinite(s_values);
  LogicalVector check_y0_s = check_y0_s_na & check_y0_s_nan & check_y0_s_inf;
  bool b_check_y0 = is_true(all(check_y0)) & is_true(all(check_y0_l)) & is_true(all(check_y0_s));
  int count = 0;
  double step_size = 1;
  NumericVector step_record;
  bool StopFlag = TRUE;
  if (b_check_y0) {
    while (StopFlag) {
      count++;
      step_record.push_back(step_size);
      NumericVector feasible_direction = score_fun(context, y0);
      double feasible_direction_sq = 0, dn = 0;
      for (int i = 0; i < y0.size(); i++) {
        feasible_direction_sq += pow(feasible_direction[i], 2);
        dn += pow(y0[i], 2);
      }
      double feasible_direction_length = sqrt(feasible_direction_sq);
      double y0_length = sqrt(dn);
      double fd_norm;
      if (y0_length != 0) {
        fd_norm = feasible_direction_length/y0_length;
      } else {
        fd_norm = feasible_direction_length;
      }
      bool check_converge_y0 = (fd_norm <= tole);
      if (check_converge_y0) {
        NumericVector est = y0;
        bool conv = check_converge_y0;
        int final_count = count;
        est_out["est"] = est;
        est_out["converge"] = conv;
        est_out["final_count"] = final_count;
        double l_value = like(context, y0);
        est_out["likelihood"] = l_value;
        est_out["gradient"] = feasible_direction;
        std::string c_status = std::to_string(count) + "-th iteration";
        est_out["converge_status"] = c_status;
        est_out["step_size"] = step_record;
        est_out["feasible_direction_length"] = feasible_direction_length;
        est_out["feasible_direction_norm"] = fd_norm;
        StopFlag = FALSE;
      } else {
        NumericVector y1 = y0+step_size*feasible_direction;
        LogicalVector check_y1_na = !is_na(y1);
        LogicalVector check_y1_nan = !is_nan(y1);
        LogicalVector check_y1_inf = !is_infinite(y1);
        LogicalVector check_y1 = check_y1_na & check_y1_nan & check_y1_inf;
        NumericVector l_value;
        l_value[0] = like(context, y1);
        LogicalVector check_y1_l_na = !is_na(l_value);
        LogicalVector check_y1_l_nan = !is_nan(l_value);
        LogicalVector check_y1_l_inf = !is_infinite(l_value);
        LogicalVector check_y1_l = check_y1_l_na & check_y1_l_nan & check_y1_l_inf;
        NumericVector s_values = score_fun(context, y1);
        LogicalVector check_y1_s_na = !is_na(s_values);
        LogicalVector check_y1_s_nan = !is_nan(s_values);
        LogicalVector check_y1_s_inf = !is_infinite(s_values);
        LogicalVector check_y1_s = check_y1_s_na & check_y1_s_nan & check_y1_s_inf;
        bool b_check_y1 = is_true(all(check_y1)) & is_true(all(check_y1_l)) & is_true(all(check_y1_s));
        bool check_stop = (count == stop_count);
        if (check_stop) {
          NumericVector est = y0;
          bool conv = check_converge_y0;
          int final_count = count;
          est_out["est"] = est;
          est_out["converge"] = conv;
          est_out["final_count"] = final_count;
          double l_value = like(context, y0);
          est_out["likelihood"] = l_value;
          est_out["gradient"] = feasible_direction;
          std::string c_status = std::to_string(count) + "-th iteration";
          est_out["converge_status"] = c_status;
          est_out["step_size"] = step_record;
          est_out["feasible_direction_length"] = feasible_direction_length;
          est_out["feasible_direction_norm"] = fd_norm;
          StopFlag = FALSE;
        } else {
          if (b_check_y1) {
            double w1 = like(context, y1);
            double w0 = like(context, y0);
            bool cp = w1 > w0;
            if (cp) {
              step_size = 1;
              y0 = y1;
            } else {
              step_size /= 2;
            }
          } else {
            NumericVector est = y0;
            bool conv = b_check_y1;
            int final_count = 0;
            est_out["est"] = est;
            est_out["converge"] = conv;
            est_out["final_count"] = final_count;
            double l_value = like(context, y0);
            est_out["likelihood"] = l_value;
            est_out["gradient"] = feasible_direction;
            std::string c_status = "Not feasible solution at " + std::to_string(count) + "-th iteration";
            est_out["converge_status"] = c_status;
            est_out["step_size"] = step_record;
            est_out["feasible_direction_length"] = feasible_direction_length;
            est_out["feasible_direction_norm"] = fd_norm;
            StopFlag = FALSE;
          }
        }
      }
    }
  } else {
    NumericVector est = y0;
    bool conv = b_check_y0;
    int final_count = 0;
    est_out["est"] = est;
    est_out["converge"] = conv;
    est_out["final_count"] = final_count;
    double l_value = like(context, y0);
    est_out["likelihood"] = l_value;
    est_out["gradient"] = s_values;
    std::string c_status = "Not feasible solution";
    est_out["converge_status"] = c_status;
    double unavailable_sign = -1;
    NumericVector unavailable_vector;
    unavailable_vector[0] = -1;
    est_out["step_size"] = unavailable_vector;
    est_out["feasible_direction_length"] = unavailable_sign;
    est_out["feasible_direction_norm"] = unavailable_sign;
  }
  return est_out;
}
