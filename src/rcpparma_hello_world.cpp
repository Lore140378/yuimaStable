#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
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

// [[Rcpp::export]]
NumericVector mysum(int lengthx, arma::rowvec x, arma::rowvec a, arma::colvec b, double pi) {
  NumericVector mydens(lengthx);
  //for (int i = 0;i<lengthx;i++){
  //  mydens(i) = arma::as_scalar(a * cos(b*x(i)));
  //}
  
  mydens = a * cos(b*x);
  return mydens/pi;
  // as.numeric(t(nodes^(1/beta0-1)/beta0*weights)%*%cos(nodes^(1/beta0)%*%t(x))/pi)
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector mytail(int lengthx, arma::vec x, arma::vec exp1, arma::vec V, arma::vec logw) {
  NumericVector dens(lengthx);
  arma::colvec unitd(exp1.n_elem,1);
  arma::colvec dummyx(exp1.n_elem);
  arma::colvec dummy(exp1.n_elem);
  arma::vec dummy2(exp1.n_elem);
  arma::vec dummy3(exp1.n_elem);
  for (int i = 0;i<lengthx;i++){
    arma::colvec dummyx = unitd.ones()*arma::as_scalar(x[i]);
    dummy = exp(exp1%log(dummyx));//pow(dummyx,exp1);
    dummy2= dummy%V;
    dummy3=exp(-dummy2+logw);
    dens(i) = arma::as_scalar((dummy2.t()*dummy3));
  }
  return dens;
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List CompMatAsym(int pos, arma::mat Deriv_A, arma::vec estim_jump, arma::mat logDeriv_C){
  int nrA = Deriv_A.n_rows;
  int nrC = logDeriv_C.n_rows;
  arma::mat dumA(nrA,nrA);
  arma::mat dumC(nrC+1,nrC+1);
  arma::mat o(1,1);
  for (int i = 0;i<pos;i++){
    arma::vec da = Deriv_A.col(i);
    dumA = dumA+da*da.t()/arma::as_scalar(estim_jump(i)*estim_jump(i));
    arma::colvec dc = logDeriv_C.col(i);
    arma::mat matc = dc*dc.t();
    arma::mat Gmatc0 = arma::join_rows(dc,matc);
    arma::mat Gmatc1 = arma::join_rows(o.ones(),dc.t());
    arma::mat Gmatc = arma::join_cols(Gmatc1,Gmatc0);
    dumC = dumC+Gmatc;
  }
  return Rcpp::List::create(
    Rcpp::Named("avec") = dumA,
    Rcpp::Named("amat") = dumC
  );
}
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List CompMatAsym1(int pos, arma::mat Deriv_A, arma::vec estim_jump, arma::mat logDeriv_C, arma::vec W){
  int nrA = Deriv_A.n_rows;
  int nrC = logDeriv_C.n_rows;
  double dumW;
  arma::vec da(nrA);
  arma::mat dumA(nrA,nrA);
  arma::mat dumAinc(nrA,nrA);
  arma::mat dumA2(nrA,nrA);
  arma::colvec dc(nrC);
  arma::mat dumC(nrC+1,nrC+1);
  arma::mat dumC2(nrC+1,nrC+1);
  arma::mat o(1,1);
  for (int i = 0;i<pos;i++){
    da = Deriv_A.col(i);
    dumW = arma::as_scalar(W[i]);
    dumAinc = da*da.t()/arma::as_scalar(estim_jump(i)*estim_jump(i))*dumW;
    dumA = dumA + dumAinc;
    dumA2 = dumA2 + dumAinc*dumW;
    dc = logDeriv_C.col(i);
    arma::mat matc = dc*dc.t();
    arma::mat Gmatc0 = arma::join_rows(dc,matc);
    arma::mat Gmatc1 = arma::join_rows(o.ones(),dc.t());
    arma::mat Gmatc = arma::join_cols(Gmatc1,Gmatc0);
    arma::mat dumCinc = Gmatc*dumW;
    dumC = dumC+dumCinc;
    dumC2 = dumC2+dumCinc*dumW;
  }
  return Rcpp::List::create(
    Rcpp::Named("avec") = dumA,
    Rcpp::Named("amat") = dumC,
    Rcpp::Named("avec2") = dumA2,
    Rcpp::Named("amat2") = dumC2
  );
}

/*
 * for(i in pos){
 dumA <- Deriv_A[,i]%*%t(Deriv_A[,i])/estim_jump[i]^2*W[i]
 dumA2 <- dumA*W[i] 
 dumC <- rbind(cbind(1,t(logDeriv_C[,i])),
 cbind(logDeriv_C[,i],logDeriv_C[,i]%*%t(logDeriv_C[,i])))*W[i]
 
#           dumC <- logDeriv_C[,i]%*%t(logDeriv_C[,i])*W[i] 
 dumC2 <- dumC*W[i]
 first_comp<- first_comp+dumA
 first_comp2 <- first_comp2+dumA2
 second_comp <-  second_comp+dumC
 second_comp2 <- second_comp2+dumC2
 }
 * 
 */

/*
 for(i in pos){
 dumA <- Deriv_A[,i]%*%t(Deriv_A[,i])/estim_jump[i]^2
 
 dumC <- rbind(cbind(1,t(logDeriv_C[,i])),
 cbind(logDeriv_C[,i],logDeriv_C[,i]%*%t(logDeriv_C[,i])))
 first_comp<- first_comp+dumA
 second_comp <-  second_comp+dumC
 }
 /*
 */
