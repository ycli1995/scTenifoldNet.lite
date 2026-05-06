#include <RcppEigen.h>

#include "data_manipulation.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd matrix_rc_norm(Eigen::MatrixXd data, double scale_factor) {
  for (int i = 0; i < data.cols(); ++i) {
    double sum = data.col(i).sum();
    for (int j = 0; j < data.rows(); ++j) {
      data(j, i) = data(j, i) / sum * scale_factor;
    }
  }
  return data;
}
