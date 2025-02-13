/* Gibbs sampling for multivariate ARS*/

#include "ars.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Gibbs sampling for multivariate ARS
//'
//' @param n number of variates to draw
//' @param initialPoints this can be a matrix in multivariate case
//'
// [[Rcpp::export]]
arma::mat ars_gibbs(
  int n,
  arma::vec initialPoints,
  arma::vec minRange,
  arma::vec maxRange,
  
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) 
{

  unsigned int p = xis.n_elem;

  // initial abscissae
  /*arma::vec initT = initialPoints;

  if (initialPoints.n_elem == 1) 
  {
    initT.insert_rows(1, initialPoints[0]*arma::ones<arma::vec>(p-1));
  }*/
  
  // boundaries of the support of h
  /*
  arma::vec minD = minRange;
  arma::vec maxD = maxRange;
  
  if (minRange.n_elem == 1)
  {
    //arma::vec minD(p, arma::fill::value(minRange));
    minD.insert_rows(1, minRange[0]*arma::ones<arma::vec>(p-1));
  }
  if (maxRange.n_elem == 1)
  {
    //arma::vec maxD(p, arma::fill::value(maxRange));
    maxD.insert_rows(1, maxRange[0]*arma::ones<arma::vec>(p-1));
  } 
  
  std::cout << "...Debug minD=" << minD.t() << "\n" <<
  "...maxD=" << maxD.t() << "\n";
  */

  arma::mat samp = arma::zeros<arma::mat>(n, p);
  for (unsigned int i = 0; i < n; ++i)
  {
    for (unsigned int j = 0; j < p; ++j)
    {
      // # here 'initialPoints' can be a vector/matrix for initializing meshgrid values
      double minD = minRange[0]; double maxD = maxRange[0];
      arma::vec tmp = ars(1, initialPoints, minD, maxD,
        j,
        xis, 
        vA, 
        vB, 
        datX0, 
        datProportion, 
        datEvent, 
        weibullS);
        /*arma::vec tmp = ars(1, initialPoints, maxD[j], maxD[j],
          j,
          xis, 
          vA, 
          vB, 
          datX0, 
          datProportion, 
          datEvent, 
          weibullS);*/

          samp(i, j) = tmp[0];
          xis[j] = tmp[0];
    }
  }

  return samp;
}