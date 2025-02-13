/* header file for ars function */

#ifndef ARS_H
#define ARS_H

#include <cmath>
#include <vector>
#include <utility>

#include <RcppArmadillo.h>

using namespace std;

// (point value, gradient value)
//typedef pair<double, double> pointgrad;

//typedef std::function<pointgrad(double x)> log_dens_callback;

typedef struct {
  int k;
  vector<double> T;
  vector<double> H;
  vector<double> Hprime;
  vector<double> Z;
  vector<double> Q;
} algo_state;

typedef std::vector<double> stdvec;

arma::vec log_dens_xi0(double par, 
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS);

// a 2-vector with un-normalized log density and its first derivative
arma::vec log_dens_xi(
  double par, 
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS);

arma::vec ars_internal(
  int n,
  std::vector<double> initAbscissae,
  double minD,
  double maxD,
  algo_state* state,// = NULL,
  long maxiter,// = 0, 

  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS);

arma::vec ars(
  int n,
  arma::vec initialPoints,
  double minRange,
  double maxRange,
  
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS);

Rcpp::List ars_debug(
  int n,
  arma::vec initialPoints,
  double minRange,
  double maxRange,
  long maxiter,
  
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS);


#endif