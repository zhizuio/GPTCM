// Gibbs sampling for variance parameters

#include "simple_gibbs.h"
#include <stdio.h>


// update \xi's variance vSq
double sampleV(
  int n, 
  double vA, 
  double vB, 
  arma::vec xi
)
{
  xi.shed_row(0);
  double vA_post = wA + 0.5 * arma::accu(xi != 0.);
  double vB_post = wB + 0.5 * arma::accu(xi % xi);

  double vSq = 1. / R::rgamma(vA_post, 1. / vB_post);

  return vSq;
}

// update \zetas' variance wSq
double sampleW(
  int n, 
  double wA, 
  double wB, 
  arma::mat zetas
)
{
  double wA_post = wA + 0.5 * arma::accu(zetas != 0.);
  double wB_post = wB + 0.5 * arma::accu(zetas % zetas);

  double wSq = 1. / R::rgamma(wA_post, 1. / wB_post);

  return wSq;
}

// update \zeta0's variance w0Sq
double sampleW(
  int n, 
  double wA, 
  double wB, 
  arma::vec zeta0
)
{
  double wA_post = wA + 0.5 * arma::accu(zetas != 0.);
  double wB_post = wB + 0.5 * arma::accu(zetas % zetas);

  double wSq = 1. / R::rgamma(wA_post, 1. / wB_post);

  return wSq;
  
}

// update \betas' variance tauSq
double sampleTau(
  int n, 
  double tauA, 
  double tauB, 
  arma::mat betas
)
{
  double tauA_post = tauA + 0.5 * arma::accu(betas != 0.);
  double tauB_post = tauB + 0.5 * arma::accu(betas % betas);

  double tauSq = 1. / R::rgamma(tauA_post, 1. / tauB_post);

  return tauSq;
}