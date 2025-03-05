/* header file for updating variances using classical Gibbs sampler */

#ifndef SIMPLE_GIBBS_H
#define SIMPLE_GIBBS_H

#include <cmath>
#include <RcppArmadillo.h>


double sampleV(
  int n, 
  double vA, 
  double vB, 
  arma::vec xi
); 

double sampleW(
  int n, 
  double wA, 
  double wB, 
  arma::mat zetas
); 

double sampleW(
  int n, 
  double wA, 
  double wB, 
  arma::vec zeta0
); 

double sampleTau(
  int n, 
  double tauA, 
  double tauB, 
  arma::mat betas
); 

#endif
