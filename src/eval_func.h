/* header file for evaluation functions (i.e. log densities) */

#ifndef EVAL_FUNC_H
#define EVAL_FUNC_H

  #include <stdio.h>
  #include <RcppArmadillo.h>


  inline double upperbound = 700.;
  inline double upperbound3 = 170.;
  inline double lowerbound = 1.0e-10;

  typedef std::vector<double> stdvec;

    
  typedef struct common_data
  {
    // members
    double *currentPars;
    
    int jj;
    int l;
    int p;
    int L;
    int N;

    double vA;
    double vB;
    double tauSq;
    double w0Sq;
    double wSq;
    double phi;
    double Delta;
    bool dirichlet;

    double *xl;
    double *xr;

    double kappa;
    double kappaA;
    double kappaB;
    bool invGamma;
    double *datTheta; 
    double *datMu; 
    double *datX;
    double *datProportion;
    double *datProportionConst;
    int *datEvent;
    double *datTime;
    double *weibullS;
    double *weibullLambda;
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
  
  double log_dens_phi
  (
    double par, 
    void *abc_data
  );

  double log_dens_kappa
  (
    double par, 
    void *abc_data
  );

  double pdfTruncNorm
  (
    double x,
    double m, 
    double sd, 
    double lower, 
    double upper
  );

#endif