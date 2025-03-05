/* header file for the main MCMC loop */

#ifndef DRIVE_H
#define DRIVE_H

  #include <stdio.h>
  #include <RcppArmadillo.h>
    
  typedef struct hyperpar_data
  {
    // members
    double vSq;
    double vA;
    double vB;
    double tauSq;
    double tauA;
    double tauB;
    double wSq;
    double wA;
    double wB;
    double w0Sq;
    double w0A;
    double w0B;

    double kappaA;
    double kappaB;
  }hyperparS;

  typedef struct range_data
  {
    // members
    double xiMin;
    double xiMax;
    double zetaMin;
    double zetaMax;
    double kappaMin;
    double kappaMax;
    double betaMin;
    double betaMax;
  }rangeS;

#endif