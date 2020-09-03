//' @export
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// getd<-function(X,v,s2,tau2,y){
//   Vinv<-t(X)%*%diag(v)%*%X/s2+diag(ncol(X))/tau2
//   Vinv.chol<-chol(Vinv)
//   V<-chol2inv(Vinv.chol)#solve(Vinv)
//   Vinv.ldet<-sum(log(diag(Vinv.chol)))#determinant(Vinv)$mod
//   bhat<-V%*%t(X)%*%diag(v)%*%y/s2 #regression eqn
//   d <- -.5*Vinv.ldet + .5*t(bhat)%*%Vinv%*%bhat
//   return(list(d=d,bhat=bhat,Vinv.ldet=Vinv.ldet,V=V,Vinv.chol=Vinv.chol,Vinv=Vinv))
// }

// [[Rcpp::export]]
Rcpp::List getdC(const arma::mat& X, const arma::colvec& v, const double& s2, const double& tau2, const arma::colvec& y) {
  int ncX = X.n_cols;
  arma::mat I_ncx = eye(ncX,ncX);
  arma::mat diagv = diagmat(v);
  arma::mat Vinv = trans(X) * diagv * X/s2 + I_ncx/tau2;
  arma::mat Vinvchol = chol(Vinv);
  arma::vec dgVinvchol = diagvec(Vinvchol);
  arma::mat V = Vinv.i(); // need chol2inv()
  double Vinvldet = sum(log(dgVinvchol));
  arma::vec bhat = V * trans(X) * diagv * y/s2;
  arma::vec d = (-0.5 * Vinvldet) + (0.5 * trans(bhat) * Vinv * bhat);
  Rcpp::List L = Rcpp::List::create(d,bhat,Vinvldet,V,Vinvchol,Vinv);
  return L;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//


