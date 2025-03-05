/* header file for univariate and multivariate arms for all parameters */

#ifndef ARMS_GIBBS_H
#define ARMS_GIBBS_H

#include <cmath>

#include "arms.h"
#include "eval_func.h"

#include <RcppArmadillo.h>


arma::mat arms_gibbs_xi(
  int n,
  int nsamp, 
  int ninit,
  double minRange,
  double maxRange,

  int metropolis, 
  bool simple,
  double convex,
  int npoint,
  
  arma::vec currentPars, 
  double vA, 
  double vB, 
  arma::mat datX0, 
  arma::mat datProportion, 
  arma::ivec datEvent, 
  arma::mat weibullS
); 

arma::mat arms_gibbs_beta
( 
  int n,
  int nsamp, 
  int ninit,
  double minRange,
  double maxRange,

  int metropolis, 
  bool simple,
  double convex,
  int npoint,
  
  arma::mat currentPars, 
  double tauSq,
  double kappa,
  arma::cube datX, 
  arma::vec datTheta,
  arma::mat datMu,
  arma::mat datProportion, 
  arma::ivec datEvent, 
  arma::vec datTime,
  arma::mat weibullS
); 

arma::mat arms_gibbs_zeta
( 
  int n,
  int nsamp, 
  int ninit,
  double minRange,
  double maxRange,

  int metropolis, 
  bool simple,
  double convex,
  int npoint,
  
  arma::mat currentPars, 
  double w0Sq, 
  double wSq, 
  double kappa,
  bool dirichlet,
  arma::cube datX, 
  arma::vec datTheta,
  arma::mat datProportionConst,
  arma::ivec datEvent, 
  arma::mat weibullS, 
  arma::mat weibullLambda
);

arma::vec arms_phi
( 
  int n,
  int nsamp, 
  int ninit,
  double minRange,
  double maxRange,

  int metropolis, 
  bool simple,
  double convex,
  int npoint,
  
  double currentPars, 
  double Delta,
  arma::mat datProportion,
  arma::mat datProportionConst
); 

arma::vec arms_kappa
( 
  int n,
  int nsamp, 
  int ninit,
  double minRange,
  double maxRange,

  int metropolis, 
  bool simple,
  double convex,
  int npoint,
  
  double currentPars, 
  double kappaA,
  double kappaB,
  bool invGamma,
  arma::vec datTheta,
  arma::mat datMu,
  arma::mat datProportion, 
  arma::ivec datEvent, 
  arma::vec datTime
);

#endif
