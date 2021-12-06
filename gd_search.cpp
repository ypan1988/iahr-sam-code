#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

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
  arma::vec val_out_vec(idx_in.n_elem, arma::fill::zeros);
  for (std::size_t i = 0; i < idx_in.n_elem; ++i) {
    val_out_vec(i) = remove_one(A, i, u_idx_in);
  }
  int idx_rm = val_out_vec.index_max();

  // compute the new A without the index of idx_rm
  arma::mat rm1A = remove_one_mat(A, idx_rm);

  // remove index idx_rm from idx_in
  idx_in.shed_row(idx_rm);

  // find one index from idx_out to add (swap)
  // which results in largest val of add_one()
  arma::vec val_in_vec(idx_out.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < idx_out.n_elem; ++i) {
    arma::uword ii = idx_out(i);
    val_in_vec(i) =
          add_one(rm1A, sig(ii, ii), sig.submat(idx_in, arma::uvec({ii})),
                  u.elem(arma::join_cols(idx_in, arma::uvec({ii}))));
  }
  arma::uword idx_swap = idx_out(val_in_vec.index_max());

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

// [[Rcpp::export]]
double ChooseSwapRobust(arma::uword nlist, arma::uvec idx_in, const arma::mat &A_list,
                        const arma::mat &sig_list, const arma::vec &u_list, const arma::vec &weights,
                        arma::uvec &out2, arma::mat &out3) {
  const arma::uword u_nrows = u_list.n_rows / nlist;
  const arma::uword A_nrows = A_list.n_rows / nlist;

  // generate the complete index
  arma::vec idx = arma::linspace(0, u_nrows - 1, u_nrows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);

  // get the index not included in complete index
  arma::uvec idx_out = std_setdiff(uidx, idx_in);

  // find one index from idx_in to remove
  // which results in largest val of remove_one()
  arma::mat val_out_mat(idx_in.n_elem, nlist, arma::fill::zeros);
#pragma omp parallel for
  for (std::size_t j = 0; j < nlist; ++j) {
    arma::mat A = A_list.rows(j*A_nrows, (j+1)*A_nrows-1);
    arma::vec u = u_list.rows(j*u_nrows, (j+1)*u_nrows-1);
    arma::vec u_idx_in = u.elem(idx_in);
    for (std::size_t i = 0; i < idx_in.n_elem; ++i) {
      val_out_mat(i,j) = remove_one(A, i, u_idx_in);
    }
  }

  arma::vec val_out_vec = val_out_mat * weights;
  int idx_rm = val_out_vec.index_max();

  // compute the new A without the index of idx_rm
  arma::mat rm1A_list;
  for (std::size_t j = 0; j < nlist; ++j) {
    arma::mat A = A_list.rows(j*A_nrows, (j+1)*A_nrows-1);
    rm1A_list = arma::join_cols(rm1A_list, remove_one_mat(A, idx_rm));
  }

  // remove index idx_rm from idx_in
  idx_in.shed_row(idx_rm);

  // find one index from idx_out to add (swap)
  // which results in largest val of add_one()
  arma::mat val_in_mat(idx_out.n_elem,nlist,arma::fill::zeros);
#pragma omp parallel for
  for (arma::uword j = 0; j < nlist; ++j) {
    arma::mat sig = sig_list.rows(j*u_nrows, (j+1)*u_nrows-1);
    arma::vec u = u_list.subvec(j*u_nrows, (j+1)*u_nrows-1);
    arma::mat rm1A = rm1A_list.rows(j*(A_nrows-1), (j+1)*(A_nrows-1)-1);
    for (arma::uword i = 0; i < idx_out.n_elem; ++i) {
      arma::uword ii = idx_out(i);
      val_in_mat(i,j) =
        add_one(rm1A, sig(ii, ii), sig.submat(idx_in, arma::uvec({ii})),
                u.elem(arma::join_cols(idx_in, arma::uvec({ii}))));
    }
  }

  arma::vec val_in_vec = val_in_mat * weights;
  arma::uword indexmax = val_in_vec.index_max();
  arma::uword idx_swap = idx_out(indexmax);

  // compute the new A with the index of idx_swap
  arma::mat newA_list;
  for (std::size_t j = 0; j < nlist; ++j) {
    arma::mat sig = sig_list.rows(j*u_nrows, (j+1)*u_nrows-1);
    arma::mat rm1A = rm1A_list.rows(j*(A_nrows-1), (j+1)*(A_nrows-1)-1);
    newA_list = arma::join_cols(newA_list, add_one_mat(rm1A, sig(idx_swap, idx_swap),
                                                       sig.submat(idx_in, arma::uvec({idx_swap}))));
  }

  idx_in = arma::join_cols(idx_in, arma::uvec({idx_swap}));

  out2 = idx_in;
  out3 = newA_list;

  return val_in_vec(indexmax);
}

// [[Rcpp::export]]
arma::uvec GradRobust(arma::uword nlist, arma::uvec idx_in, arma::mat A_list,
                      arma::mat sig_list, arma::vec u_list, arma::vec weights,
                      double tol = 1e-9, bool trace = true) {
  const arma::uword u_nrows = u_list.n_rows / nlist;
  const arma::uword A_nrows = A_list.n_rows / nlist;

  arma::vec new_val_vec(nlist, arma::fill::zeros);
  for (arma::uword j = 0; j < nlist; ++j) {
    arma::mat A = A_list.rows(j*A_nrows, (j+1)*A_nrows-1);
    arma::vec u = u_list.subvec(j*u_nrows, (j+1)*u_nrows-1);
    new_val_vec(j) = obj_fun(A, u.elem(idx_in));
  }
  double new_val = arma::dot(new_val_vec, weights);

  double diff = 1.0;
  int i = 0;
  while (diff > tol) {
    double val = new_val;
    i = i + 1;

    arma::uvec out2;
    arma::mat out3;
    new_val = ChooseSwapRobust(nlist, idx_in, A_list, sig_list, u_list, weights,
                               out2, out3);

    diff = new_val - val;
    if (diff > 0) {
      A_list = out3;
      idx_in = out2;
    }
    if (trace) Rcpp::Rcout << "\rIter: " << i << " " << diff << std::endl;
  }
  return (idx_in);
}
