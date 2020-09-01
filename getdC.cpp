#include <Rcpp.h>
# include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

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
List getdC(NumericMatrix X, NumericMatrix v, double s2, double tau2, NumericVector y, Function chol, Function chol2inv, Function diag) {
  int ncX = X.n_cols;
  NumericMatrix Vinv = trans(X) %*% diag(v) %*% X/s2 + diag(ncX)/tau2;
  NumericMatrix Vinv.chol = chol(Vinv);
  NumericMatrix V = chol2inv(Vinv.chol);
  NumericMatrix Vinv.ldet = sum(log(diag(Vinv.chol)));
  NumericMatrix bhat = V %*% t(X) %*% diag(v) %*% y/s2;
  double d = -.5*Vinv.ldet + .5*t(bhat)%*%Vinv%*%bhat;
  List L = List::create(d,bhat,Vinv.ldet,V,Vinv.chol,Vinv);
  return L;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


