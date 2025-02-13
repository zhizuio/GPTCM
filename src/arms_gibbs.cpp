// Gibbs sampling for multivariate ARMS

#include <cmath>

#include "arms.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// The following evaluation functions might have to be in C not C++ and be put into the extern arms C code
double log_dens_betas(double par, 
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) {

    double logpost = 0.;

    return logpost;
}


double myfunc(double x, void *mydata) {
  int a = 0;
  // TBC...
}

//' Multivariate ARMS via Gibbs sampler
//'
//' @param n number of variates to draw
//' @param initialPoints this can be a matrix in multivariate case
//'
// [[Rcpp::export]]
arma::mat arms_gibbs(
  int n,
  arma::vec initialPoints,
  arma::vec minRange,
  arma::vec maxRange,
  bool metropolis, 
  
  arma::vec& betaJ, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) 
{

  unsigned int p = betaJ.n_elem;

  int ninit = 1;
  double anArray[] = {5., 16.};
  double *xl = anArray;
  double *xr = anArray;
  //double myfunc = .1;
  void *mydata = nullptr;
  int dometrop = 1;
	double *xprev = anArray;
  double *xsamp = anArray;

  arma::mat samp = arma::zeros<arma::mat>(n, p);
  for (unsigned int i = 0; i < n; ++i)
  {
    for (unsigned int j = 0; j < p; ++j)
    {
      // # here 'initialPoints' can be a vector/matrix for initializing meshgrid values
      double minD = minRange[0]; double maxD = maxRange[0];
      
      double tmp = ARMS::arms_simple (
        ninit, xl,  xr,
        myfunc, mydata,
        dometrop, xprev, xsamp);

      samp(i, j) = tmp;
      betaJ[j] = tmp;
    }
  }

  return samp;
}