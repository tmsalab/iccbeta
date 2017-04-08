#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title Compute ICCbeta intraclass correlation
//' @description Computes the proportion of variance attributed to group heterogeneity in slopes as described by Aguinis & Culpepper (2015).
//' The function uses lmer model information. 
//' @usage icc_beta(X,l2id,T,vy)
//' @param X The design \code{matrix} of fixed effects from a lmer model.  
//' @param l2id A \code{vector} that identifies group membership. The vector must be coded as a sequence from 1 to J, the number of groups.
//' @param T A \code{matrix} of the estimated variance-covariance matrix of a lmer model fit.
//' @return vy The variance of the dependent variable.
//' @author Steven A Culpepper

// [[Rcpp::export]]
Rcpp::List icc_beta(const arma::mat& X,const arma::vec& l2id,const arma::mat& T,double vy){
  unsigned int N = l2id.n_elem;
  unsigned int p = X.n_cols;
  unsigned int J=max(l2id);
  arma::mat SUMS = arma::zeros<arma::mat>(J,p);
  arma::mat means = arma::zeros<arma::mat>(J,p);
  arma::vec Nj = arma::zeros<arma::vec>(J);
  arma::vec ONE = arma::ones<arma::vec>(N);
  
  for(unsigned int i=0;i<N;i++){
    //assume l2id goes from 0 to J-1
    Nj(l2id(i)-1) += ONE(i);
    SUMS.row(l2id(i)-1) += X.row(i);
  }
  
  for(unsigned int j=0;j<J;j++){
    if(Nj(j)>0.0){
      means.row(j) = SUMS.row(j)/Nj(j);
    }
  }

  arma::mat XcpXc = arma::zeros<arma::mat>(p,p);
  for(unsigned int i=0;i<N;i++){
    XcpXc += ( X.row(i) - means.row(l2id(i)-1) ).t() * ( X.row(i) - means.row(l2id(i)-1) );
  }

  double rho_beta = trace(T * XcpXc/(sum(Nj)-1.0))/vy;  
  
  return Rcpp::List::create(Rcpp::Named("J",J),
                            Rcpp::Named("means",means),
                            Rcpp::Named("XcpXc",XcpXc),
                            Rcpp::Named("Nj",Nj),
                            Rcpp::Named("rho_beta",rho_beta)
                            );
}
