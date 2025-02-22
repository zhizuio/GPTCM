/* header file for evaluation functions (i.e. log densities) */

#ifndef EVAL_FUNC_H
#define EVAL_FUNC_H

  #include <stdio.h>
  #include <RcppArmadillo.h>


  inline double upperbound = 700.;

  typedef std::vector<double> stdvec;

    
  typedef struct common_data
  {
    // members
    double *currentPars;
    unsigned int jj;

    int l;
    int p;
    int L;
    int N;
    double tauSq;
    double kappa;
    double *datTheta; 
    double *datMu; 
    double *datTime; 

    double *xl;
    double *xr;
    double vA;
    double vB;
    double *datX;
    double *datProportion;
    int *datEvent;
    double *weibullS;
  }dataS;


  double log_dens_xis
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