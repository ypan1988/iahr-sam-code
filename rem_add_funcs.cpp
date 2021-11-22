#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double obj_fun(arma::mat A, arma::vec U2) {
  return arma::as_scalar(U2.t() * A * U2);
}

// [[Rcpp::export]]
double remove_one(arma::mat A, arma::uword i, arma::vec u) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  arma::uvec i2(1);
  i2(0) = i;
  uidx.shed_row(i);
  // arma::mat A2(arma::size(A));

  double d = A(i, i);
  arma::vec b = A.submat(uidx, i2);

  // step 1: compute A*
  A = A.submat(uidx, uidx) - b * b.t() / d;

  // return obj_fun(A,u(uidx));
  return arma::as_scalar(u(uidx).t() * A * u(uidx));
}

// [[Rcpp::export]]
arma::mat remove_one_mat(arma::mat A, arma::uword i) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  arma::uvec i2(1);
  i2(0) = i;
  uidx.shed_row(i);
  // arma::mat A2(arma::size(A));

  double d = A(i, i);
  arma::vec b = A.submat(uidx, i2);

  // step 1: compute A*
  A = A.submat(uidx, uidx) - b * b.t() / d;

  // return obj_fun(A,u(uidx));
  return A;
}

// [[Rcpp::export]]
double add_one(const arma::mat &A, double sigma_jj, const arma::vec &f,
               const arma::vec &u) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1);
  // step 1: compute A*
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  for (arma::uword j = 0; j < A2.n_rows - 1; j++) {
    A2(j, A2.n_cols - 1) = 0;
    A2(A2.n_rows - 1, j) = 0;
  }
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;

  // // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
       ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));

  // // // step 4: compute K3_inv
  A2 = A2 -
       ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));

  return obj_fun(A2, u);
}

// [[Rcpp::export]]
arma::mat add_one_mat(arma::mat A, double sigma_jj, arma::vec f) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1);
  // step 1: compute A*
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  for (arma::uword j = 0; j < A2.n_rows - 1; j++) {
    A2(j, A2.n_cols - 1) = 0;
    A2(A2.n_rows - 1, j) = 0;
  }
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;

  // // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
       ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));

  // // // step 4: compute K3_inv
  A2 = A2 -
       ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));

  return A2;
}
