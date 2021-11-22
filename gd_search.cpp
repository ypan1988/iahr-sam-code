#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double obj_fun(const arma::mat &A, const arma::vec &U2) {
  return arma::as_scalar(U2.t() * A * U2);
}

// [[Rcpp::export]]
double remove_one(const arma::mat &A, arma::uword i, const arma::vec &u) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  uidx.shed_row(i);

  double d = A(i, i);
  arma::vec b = A.submat(uidx, arma::uvec({i}));

  return obj_fun(A.submat(uidx, uidx) - b * b.t() / d, u(uidx));
}

// [[Rcpp::export]]
arma::mat remove_one_mat(const arma::mat &A, arma::uword i) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  uidx.shed_row(i);

  double d = A(i, i);
  arma::vec b = A.submat(uidx, arma::uvec({i}));

  return A.submat(uidx, uidx) - b * b.t() / d;
}

// [[Rcpp::export]]
double add_one(const arma::mat &A, double sigma_jj, const arma::vec &f,
               const arma::vec &u) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1, arma::fill::zeros);
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  // for(arma::uword j =0; j < A2.n_rows-1; j++){
  //   A2(j,A2.n_cols-1) = 0;
  //   A2(A2.n_rows-1,j) = 0;
  // }
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;

  // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
       ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));

  // step 4: compute K3_inv
  A2 = A2 -
       ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));

  return obj_fun(A2, u);
}

// [[Rcpp::export]]
arma::mat add_one_mat(const arma::mat &A, double sigma_jj, const arma::vec &f) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1, arma::fill::zeros);
  // step 1: compute A*
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  for (arma::uword j = 0; j < A2.n_rows - 1; j++) {
    A2(j, A2.n_cols - 1) = 0;
    A2(A2.n_rows - 1, j) = 0;
  }
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;

  // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
       ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));

  // step 4: compute K3_inv
  A2 = A2 -
       ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));

  return A2;
}

arma::uvec std_setdiff(arma::uvec &x, arma::uvec &y) {
  std::vector<int> a = arma::conv_to<std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to<std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;

  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));

  return arma::conv_to<arma::uvec>::from(out);
}

// [[Rcpp::export]]
double ChooseSwap(arma::uvec idx_in, arma::mat A, arma::mat sig, arma::vec u,
                  arma::uvec &out2, arma::mat &out3) {
  // generate the complete index
  arma::vec idx = arma::linspace(0, sig.n_rows - 1, sig.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);

  // get the index not included in complete index
  arma::uvec idx_out = std_setdiff(uidx, idx_in);

  // find one index from idx_in to remove
  // which results in largest val of remove_one()
  arma::vec u_idx_in = u.elem(idx_in);
  int idx_rm = 0;
  double maxval = remove_one(A, 0, u_idx_in);
  for (std::size_t i = 1; i < idx_in.n_elem; ++i) {
    double val = remove_one(A, i, u_idx_in);
    if (val > maxval) idx_rm = i, maxval = val;
  }

  // compute the new A without the index of idx_rm
  arma::mat rm1A = remove_one_mat(A, idx_rm);

  // remove index idx_rm from idx_in
  idx_in.shed_row(idx_rm);

  // find one index from idx_out to add (swap)
  // which results in largest val of add_one()
  arma::uword idx_swap = 0;
  maxval = add_one(rm1A, sig(idx_out(0), idx_out(0)),
                   sig.submat(idx_in, arma::uvec({idx_out(0)})),
                   u.elem(arma::join_cols(idx_in, arma::uvec({idx_out(0)}))));

  for (arma::uword i = 1; i < idx_out.n_elem; ++i) {
    arma::uword ii = idx_out(i);
    double val =
        add_one(rm1A, sig(ii, ii), sig.submat(idx_in, arma::uvec({ii})),
                u.elem(arma::join_cols(idx_in, arma::uvec({ii}))));
    if (val > maxval) idx_swap = ii, maxval = val;
  }

  // compute the new A with the index of idx_swap
  arma::mat newA = add_one_mat(rm1A, sig(idx_swap, idx_swap),
                               sig.submat(idx_in, arma::uvec({idx_swap})));
  idx_in = arma::join_cols(idx_in, arma::uvec({idx_swap}));

  out2 = idx_in;
  out3 = newA;
  return obj_fun(newA, u.elem(idx_in));
}

// [[Rcpp::export]]
arma::uvec Grad(arma::uvec idx_in, arma::mat A, arma::mat sig, arma::vec u,
                double tol = 1e-9, bool trace = true) {
  double new_val = obj_fun(A, u.elem(idx_in));
  double diff = 1.0;
  int i = 0;
  while (diff > tol) {
    double val = new_val;
    i = i + 1;

    arma::uvec out2;
    arma::mat out3;
    new_val = ChooseSwap(idx_in, A, sig, u, out2, out3);
    diff = new_val - val;
    if (diff > 0) {
      A = out3;
      idx_in = out2;
    }
    if (trace) Rcpp::Rcout << "\rIter: " << i << " " << diff << std::endl;
  }
  return (idx_in);
}
