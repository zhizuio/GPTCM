/* header file for evaluation functions (i.e. log densities) */

#ifndef EVAL_FUNC_H
#define EVAL_FUNC_H

  #include <stdio.h>
  #include <RcppArmadillo.h>


  inline double upperbound1 = 70.;
  inline double upperbound2 = 69.;

  typedef std::vector<double> stdvec;

    
  typedef struct common_data
  {
    // members
    arma::vec currentPars;
    unsigned int jj;
    double *xl;
    double *xr;
    double vA;
    double vB;
    arma::mat datX;
    arma::mat datProportion;
    arma::uvec datEvent;
    arma::mat weibullS;
  }dataS;

  void create_mydata
  (
    arma::vec currentPars,
    unsigned int jj,
    double *xl,
    double *xr,
    double vA, 
    double vB,
    arma::mat datX,
    arma::mat datProportion,
    arma::uvec datEvent,
    arma::mat weibullS,
    dataS *abc_data
  );

  double myfunc
  (
    double par, 
    void *abc_data
  );

  double log_dens_betas
  (
    double par, 
    void *abc_data
  );

  double log_dens_zetas
  (
    double par, 
    void *abc_data
  );
  


#endif