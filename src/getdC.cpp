
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export(getdC)]]
Rcpp::List getdC(const arma::mat& X, const arma::colvec& v, const double& s2, const double& tau2, const arma::colvec& y) {
  int ncX = X.n_cols;
  mat I_ncx = eye(ncX,ncX);
  mat diagv = diagmat(v);
  mat Vinv = trans(X) * diagv * X/s2 + I_ncx/tau2;
  mat Vinvchol = chol(Vinv);
  vec dgVinvchol = diagvec(Vinvchol);
  mat V = Vinv.i(); // need chol2inv()
  double Vinvldet = sum(log(dgVinvchol));
  vec bhat = V * trans(X) * diagv * y/s2;
  vec d = (-0.5 * Vinvldet) + (0.5 * trans(bhat) * Vinv * bhat);
  Rcpp::List L = Rcpp::List::create(d,bhat,Vinvldet,V,Vinvchol,Vinv);
  return L;
}

