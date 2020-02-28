#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppEigen;

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::setNbThreads(1);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}
