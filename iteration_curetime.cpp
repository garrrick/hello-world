#include "iteration_curetime.h"
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
                        int stop_count) {
  List est_out;
  NumericVector y0 = wrap(p0);
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
  NumericVector gradient = score_fun(context, y0);
  LogicalVector check_y0_s_na = !is_na(gradient);
  LogicalVector check_y0_s_nan = !is_nan(gradient);
  LogicalVector check_y0_s_inf = !is_infinite(gradient);
  LogicalVector check_y0_s = check_y0_s_na & check_y0_s_nan & check_y0_s_inf;
  bool b_check_y0 = is_true(all(check_y0)) & is_true(all(check_y0_l)) & is_true(all(check_y0_s));
  int e_size = Te.n_elem;
  int para_size = p0.n_elem;
  uvec working_set, unworking_set;
  mat active_set, inactive_set, ainv, im, projection_matrix;
  vec bx, check_active, gra, feasible_direction, dx;
  int no_active_constraints = 0;
  int no_inactive_constraints = 0;
  int count = 0;
  bool StopFlag = TRUE;
  bool de_constraint = FALSE;
  NumericVector step_record, Armijo_record;
  StringVector phase_record;
  double step_size = 1;
  int smallest_step_index;
  if (b_check_y0) {
    while (StopFlag) {
      count++;
      if (de_constraint == FALSE) {
        bx = Xe*p0;
        check_active = abs(bx-log(Te));
        no_active_constraints = 0;
        no_inactive_constraints = 0;
        for (uword i = 0; i < e_size; i++) {
          if (check_active(i) <= pow(10,-10)) {
            no_active_constraints++;
            working_set.resize(no_active_constraints);
            working_set(no_active_constraints-1) = i;
          } else {
            no_inactive_constraints++;
            unworking_set.resize(no_inactive_constraints);
            unworking_set(no_inactive_constraints-1) = i;
          }
        }
        active_set = Xe.rows(working_set);
        inactive_set = Xe.rows(unworking_set);
        std::string phase_report = std::to_string(no_active_constraints) + " constraints";
        phase_record.push_back(phase_report);
      }
      gra = as<vec>(gradient);
      if (no_active_constraints == 0) {
        feasible_direction = gra;
      } else {
        ainv = pinv(active_set*active_set.t());
        im = eye<mat>(para_size, para_size);
        projection_matrix = im-active_set.t()*ainv*active_set;
        feasible_direction = projection_matrix*gra;
      }
      double feasible_direction_sq = 0, dn = 0;
      for (uword i = 0; i < para_size; i++) {
        feasible_direction_sq += pow(feasible_direction(i), 2);
        dn += pow(p0(i), 2);
      }
      double feasible_direction_length = sqrt(feasible_direction_sq);
      double p0_length = sqrt(dn);
      double fd_norm;
      if (p0_length != 0) {
        fd_norm = feasible_direction_length/p0_length;
      } else {
        fd_norm = feasible_direction_length;
      }
      bool check_converge_p0 = (fd_norm <= tole);
      if (check_converge_p0 == TRUE) {
        if (no_active_constraints == 0) {
          vec est = p0;
          bool conv = check_converge_p0;
          int final_count = count;
          est_out["est"] = est;
          est_out["converge"] = conv;
          std::string c_status = "no active constraint";
          est_out["converge_status"] = c_status;
          est_out["step_size"] = step_record;
          est_out["Armijo_count"] = Armijo_record;
          est_out["final_count"] = final_count;
          double l_value = like(context, y0);
          est_out["likelihood"] = l_value;
          est_out["gradient"] = gra;
          est_out["smallest_step_index"] = smallest_step_index;
          est_out["feasible_direction_length"] = feasible_direction_length;
          est_out["feasible_direction_norm"] = fd_norm;
          est_out["active_set"] = active_set;
          StopFlag = FALSE;
          std::string phase_report = "converge and " + c_status;
          phase_record.push_back(phase_report);
          est_out["phase"] = phase_record;
        } else {
          vec Lagrange_multipliers = -ainv*active_set*gra;
          double Lmm = min(Lagrange_multipliers);
          if (Lmm >= 0) {
            vec est = p0;
            bool conv = check_converge_p0;
            int final_count = count;
            est_out["est"] = est;
            est_out["converge"] = conv;
            uword nac = 0;
            for (uword i = 0; i < no_active_constraints; i++) {
              if (Lagrange_multipliers(i) > 0) {
                nac++;
              }
            }
            std::string c_status = std::to_string(nac) + " active constraints";
            est_out["converge_status"] = c_status;
            est_out["step_size"] = step_record;
            est_out["Armijo_count"] = Armijo_record;
            est_out["final_count"] = final_count;
            double l_value = like(context, y0);
            est_out["likelihood"] = l_value;
            est_out["gradient"] = gra;
            est_out["smallest_step_index"] = smallest_step_index;
            est_out["feasible_direction_length"] = feasible_direction_length;
            est_out["feasible_direction_norm"] = fd_norm;
            est_out["active_set"] = active_set;
            StopFlag = FALSE;
            std::string phase_report = "d = 0 and " + c_status;
            phase_record.push_back(phase_report);
            est_out["phase"] = phase_record;
          } else {
            uword erase_index = 0;
            for (uword i = 0; i < Lagrange_multipliers.n_elem; i++) {
              double Lmi = Lagrange_multipliers(i);
              if (Lmi == Lmm) {
                erase_index = i;
              }
            }
            uword Xe_index = working_set(erase_index);
            active_set.shed_row(erase_index);
            working_set.resize(no_active_constraints-1);
            unworking_set.resize(no_inactive_constraints+1);
            unworking_set(no_inactive_constraints) = Xe_index;
            inactive_set = Xe.rows(unworking_set);
            no_active_constraints = active_set.n_rows;
            no_inactive_constraints = inactive_set.n_rows;
            de_constraint = TRUE;
            std::string phase_report = "d = 0 and remove 1 constraint";
            phase_record.push_back(phase_report);
          }
        }
      } else {
        dx = inactive_set*feasible_direction;
        bx = inactive_set*p0;
        vec subTe = Te.elem(unworking_set);
        vec steps = (bx-log(subTe))/dx;
        NumericVector step_size_candidate;
        int iac_size = unworking_set.n_elem;
        for (uword i = 0; i < iac_size; i++) {
          if (dx(i) > 0) {
            double chosen_step = steps(i);
            step_size_candidate.push_back(chosen_step);
          }
        }
        int ssc_size = step_size_candidate.size();
        if (ssc_size > 0) {
          double step_size_min = min(step_size_candidate);
          uword ssc_index = which_min(step_size_candidate);
          smallest_step_index = unworking_set(ssc_index);
          if (step_size_min < 1) {
            step_size = step_size_min;
          } else {
            step_size = 1;
          }
        } else {
          step_size = 1;
        }
        int Armijo_count = 0;
        bool Armijo_StopFlag = TRUE;
        while (Armijo_StopFlag) {
          vec p1 = p0+step_size*feasible_direction;
          NumericVector y1 = wrap(p1);
          LogicalVector check_y1_na = !is_na(y1);
          LogicalVector check_y1_nan = !is_nan(y1);
          LogicalVector check_y1_inf = !is_infinite(y1);
          LogicalVector check_y1 = check_y1_na & check_y1_nan & check_y1_inf;
          l_value[0] = like(context, y1);
          LogicalVector check_y1_l_na = !is_na(l_value);
          LogicalVector check_y1_l_nan = !is_nan(l_value);
          LogicalVector check_y1_l_inf = !is_infinite(l_value);
          LogicalVector check_y1_l = check_y1_l_na & check_y1_l_nan & check_y1_l_inf;
          gradient = score_fun(context, y1);
          LogicalVector check_y1_s_na = !is_na(gradient);
          LogicalVector check_y1_s_nan = !is_nan(gradient);
          LogicalVector check_y1_s_inf = !is_infinite(gradient);
          LogicalVector check_y1_s = check_y1_s_na & check_y1_s_nan & check_y1_s_inf;
          bool b_check_y1 = is_true(all(check_y1)) & is_true(all(check_y1_l)) & is_true(all(check_y1_s));
          Armijo_count++;
          if (b_check_y1) {
            bool check_Armijo_stop = (Armijo_count == 30);
            if (check_Armijo_stop) {
              y0 = y1;
              p0 = p1;
              Armijo_StopFlag = FALSE;
              step_record.push_back(step_size);
              Armijo_record.push_back(Armijo_count);
            } else {
              double w1 = like(context, y1);
              double w0 = like(context, y0);
              bool cp = w1 > w0;
              if (cp) {
                y0 = y1;
                p0 = p1;
                Armijo_StopFlag = FALSE;
                step_record.push_back(step_size);
                Armijo_record.push_back(Armijo_count);
              } else {
                step_size /= 2;
              }
            }
          } else {
            vec est = p0;
            bool conv = check_converge_p0;
            int final_count = count;
            est_out["est"] = est;
            est_out["converge"] = conv;
            std::string c_status = "Not feasible solution at " + std::to_string(Armijo_count) + "-th steps in " + std::to_string(count) + "-th iteration";
            est_out["converge_status"] = c_status;
            est_out["step_size"] = step_record;
            est_out["Armijo_count"] = Armijo_record;
            est_out["final_count"] = final_count;
            double l_value = like(context, y0);
            est_out["likelihood"] = l_value;
            est_out["gradient"] = gra;
            est_out["smallest_step_index"] = smallest_step_index;
            est_out["feasible_direction_length"] = feasible_direction_length;
            est_out["feasible_direction_norm"] = fd_norm;
            est_out["active_set"] = active_set;
            Armijo_StopFlag = FALSE;
            StopFlag = FALSE;
            std::string phase_report = c_status;
            phase_record.push_back(phase_report);
            est_out["phase"] = phase_record;
          }
        }
        de_constraint = FALSE;
        std::string phase_report = "climbing with step size " + std::to_string(step_size);
        phase_record.push_back(phase_report);
      }
      bool check_stop = (count == stop_count);
      if (check_stop) {
        vec est = p0;
        bool conv = check_converge_p0;
        int final_count = count;
        est_out["est"] = est;
        est_out["converge"] = conv;
        std::string c_status = std::to_string(count) + "-th iteration";
        est_out["converge_status"] = c_status;
        est_out["step_size"] = step_record;
        est_out["Armijo_count"] = Armijo_record;
        est_out["final_count"] = final_count;
        double l_value = like(context, y0);
        est_out["likelihood"] = l_value;
        est_out["gradient"] = gra;
        est_out["smallest_step_index"] = smallest_step_index;
        est_out["feasible_direction_length"] = feasible_direction_length;
        est_out["feasible_direction_norm"] = fd_norm;
        est_out["active_set"] = active_set;
        StopFlag = FALSE;
        std::string phase_report = "stop at maximal " + c_status;
        phase_record.push_back(phase_report);
        est_out["phase"] = phase_record;
      }
    }
  } else {
    vec est = p0;
    int final_count = 0;
    est_out["est"] = est;
    bool conv = FALSE;
    est_out["converge"] = conv;
    std::string c_status = "Not feasible solution";
    est_out["converge_status"] = c_status;
    NumericVector ua_NumericVector;
    ua_NumericVector[0] = -1;
    est_out["step_size"] = ua_NumericVector;
    est_out["Armijo_count"] = ua_NumericVector;
    est_out["final_count"] = final_count;
    double l_value = like(context, y0);
    est_out["likelihood"] = l_value;
    est_out["gradient"] = gra;
    uvec ua_uvec;
    int ua_int = -1;
    est_out["smallest_step_index"] = ua_int;
    double ua_double = -1;
    est_out["feasible_direction_length"] = ua_double;
    est_out["feasible_direction_norm"] = ua_double;
    est_out["active_set"] = active_set;
    vec ua_vec;
    std::string phase_report = c_status;
    phase_record.push_back(phase_report);
    est_out["phase"] = phase_record;
  }
  return est_out;
}
