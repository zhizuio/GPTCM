// Gibbs sampling for multivariate ARMS

#include <cmath>

#include "arms.h"
#include "eval_func.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


//' Multivariate ARMS via Gibbs sampler
//'
//' @param n Number of samples to draw
//' @param every How many samples to draw for generating each sample; only the last 'n' draws will be kept; but it seems not work for every>1, not know why
//' @param ninit Number of initials as meshgrid values for envelop search
//' @param convex Adjustment for convexity (non-negative value, default 1.0)
//' @param npoint Maximum number of envelope points
//'
// [[Rcpp::export]]
arma::mat arms_gibbs_xi(
  int n,
  int every, 
  int ninit,
  arma::vec minRange,
  arma::vec maxRange,

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
  arma::mat weibullS) 
{
  /* IF every>1, all samples are the same; not know why
  if (n > every)
  {
    std::printf("Arguments 'n' should not be larger than argument 'every'.\n");
    exit (0);
  }
  */

  //int ninit = initialPoints.n_elem; // # here 'initialPoints' can be a vector/matrix (in univariate/multivariate cases) for initializing meshgrid values

  // number of parameters to be updated
  int p = currentPars.n_elem; // = 1;
  int L = datProportion.n_cols;
  int N = datProportion.n_rows;

  double minD;
  double maxD;
  minD = minRange[0]; // [j]
  maxD = maxRange[0]; // [j]
  double *xl; xl = &minD;
  double *xr; xr = &maxD;

  int dometrop = metropolis;

  //dataS *mydata = create_mydata(currentPars, j, vA, vB, datX0, datProportion, datEvent, weibullS);
  dataS *mydata = (dataS *)malloc(sizeof (dataS));
  //create_mydata(currentPars, j, xl, xr, vA, vB, datX0, datProportion, datEvent, weibullS, mydata);
  mydata->currentPars = currentPars.memptr();
  mydata->p = p;
  mydata->L = L;
  mydata->N = N;
  mydata->xl = xl;
  mydata->xr = xr;
  mydata->vA = vA;
  mydata->vB = vB;
  mydata->datX = datX0.memptr();
  mydata->datProportion = datProportion.memptr();
  mydata->datEvent = datEvent.memptr();
  mydata->weibullS = weibullS.memptr();
  
  // convert arma::mat to C double pointer or array
  // TBA...

  arma::mat samp = arma::zeros<arma::mat>(p, n + 1);
  samp.col(0) = currentPars;

  for (unsigned int i = 0; i < n; ++i)
  {
    for (unsigned int j = 0; j < p; ++j)
    {
      mydata->jj = j;
      //double initi = samp(j, i);  //samp(j, i)
      //double *xprev; xprev = &initi;
      double xprev = samp(j, i);
      double *xsamp = (double*)malloc(every * sizeof(double));
      double qcent[1], xcent[1];
      int neval, ncent = 0;

      int err;
      if (simple)
      {
        err = ARMS::arms_simple (
          ninit, xl,  xr,
          log_dens_xis, mydata,
          dometrop, &xprev, xsamp);
      } else {
        arma::vec xinit0 = arma::linspace( minD+1.0e-10, maxD-1.0e-10, ninit );
        double xinit[ninit];
        for (unsigned int i = 0; i < ninit; ++i)
          xinit[i] = xinit0[i];

        err = ARMS::arms (
          xinit, ninit, xl,  xr,
          log_dens_xis, mydata,
          &convex, npoint,
          dometrop, &xprev, xsamp,
          every, qcent, xcent, ncent, &neval);
          /*
          std::cout << "...debug arms(): ninit=" << ninit <<
          "; xinit[0]="  << xinit[0] << 
          "; xinit[4]="  << xinit[4] << 
          "; xinit[ninit-1]="  << xinit[ninit-1] << 
          "; xl=" << *xl << 
          "; xr=" << *xr << 
          "; convex=" << convex << 
          "; npoint=" << npoint << 
          "; dometrop=" << dometrop << 
          "; xprev=" << *xprev << 
          "; xsamp=" << *xsamp << 
          "; every=" << every << 
          "; qcent=" << *qcent << 
          "; xcent=" << *xcent << 
          "; ncent=" << ncent << 
          "; neval=" << neval << 
          "\n"; */
      }

      // check ARMS validity
      if (err > 0)
       std::printf("In ARMS::arms_(): error code in ARMS = %d.\n", err);
      if (isnan(xsamp[every-1]))
        std::printf("In ARMS::arms_(): NaN generated, possibly due to overflow in (log-)density (e.g. with densities involving exp(exp(...))).\n");
      if (xsamp[n-1] < minD || xsamp[every-1] > maxD)
        std::printf("In ARMS::arms_(): %d-th sample out of range [%f, %f] (fused domain). Got %f.\n", every, *xl, *xr, xsamp[every-1]);

      //xprev = xsamp[every - 1];
      //for (unsigned int i = 0; i < n; ++i) samp(j, i + 1) = xsamp[every - n + i];
      currentPars[j] = xsamp[every - 1];
      samp(j, i + 1) = xsamp[every - 1];

      free(xsamp);
    }
  }
  free(mydata);
  // remove the inital values in the first column of samp
  samp.shed_col(0);

  //std::cout << "...debug sample=" << samp.col(0).t() << "\n";

  return samp;
}



//' Multivariate ARMS via Gibbs sampler for beta
//'
//' @param n Number of samples to draw
//' @param every How many samples to draw for generating each sample; only the last 'n' draws will be kept; but it seems not work for every>1, not know why
//' @param ninit Number of initials as meshgrid values for envelop search
//' @param convex Adjustment for convexity (non-negative value, default 1.0)
//' @param npoint Maximum number of envelope points
//'
// [[Rcpp::export]]
arma::mat arms_gibbs_beta(
  int n,
  int every, 
  int ninit,
  arma::vec minRange,
  arma::vec maxRange,

  int metropolis, 
  bool simple,
  double convex,
  int npoint,
  
  arma::mat currentPars, 
  double vA, 
  double vB, 
  double tauSq,
  double kappa,
  arma::cube datX, 
  arma::vec datTheta,
  arma::mat datMu,
  arma::mat datProportion, 
  arma::ivec datEvent, 
  arma::vec datTime,
  arma::mat weibullS) 
{
  // dimensions
  int p = currentPars.n_rows; 
  int L = datProportion.n_cols;
  int N = datProportion.n_rows;
  double minD;
  double maxD;
  minD = minRange[0]; // [j]
  maxD = maxRange[0]; // [j]
  double *xl; xl = &minD;
  double *xr; xr = &maxD;

  int dometrop = metropolis;

  dataS *mydata = (dataS *)malloc(sizeof (dataS));
  mydata->currentPars = currentPars.memptr();
  mydata->p = p;
  mydata->L = L;
  mydata->N = N;
  mydata->xl = xl;
  mydata->xr = xr;
  mydata->vA = vA;
  mydata->vB = vB;
  mydata->tauSq = tauSq,
  mydata->kappa = kappa,
  mydata->datTheta = datTheta.memptr();
  mydata->datMu = datMu.memptr();
  mydata->datProportion = datProportion.memptr();
  mydata->datEvent = datEvent.memptr();
  mydata->datTime = datTime.memptr();
  mydata->weibullS = weibullS.memptr();
  
  for (int l = 0; l < L; ++l)
  {
    // Gibbs sampling
    for (int j = 0; j < p; ++j)
    {
      mydata->jj = j;
      mydata->l = l;
      mydata->datX = (datX.slice(l)).memptr();
      
      //double initi = currentPars(j, l);  //samp(j, i)
      //double *xprev; xprev = &initi;
      double xprev = currentPars(j, l);
      double *xsamp = (double*)malloc(every * sizeof(double));
      double qcent[1], xcent[1];
      int neval, ncent = 0;

      int err;
      if (simple)
      {
        err = ARMS::arms_simple (
          ninit, xl,  xr,
          log_dens_betas, mydata,
          dometrop, &xprev, xsamp);
      } else {
        arma::vec xinit0 = arma::linspace( minD+1.0e-10, maxD-1.0e-10, ninit );
        double xinit[ninit];
        for (int i = 0; i < ninit; ++i)
          xinit[i] = xinit0[i];

        err = ARMS::arms (
          xinit, ninit, xl,  xr,
          log_dens_betas, mydata,
          &convex, npoint,
          dometrop, &xprev, xsamp,
          every, qcent, xcent, ncent, &neval);
      }

      // check ARMS validity
      if (err > 0)
       std::printf("In ARMS::arms_(): error code in ARMS = %d.\n", err);
      if (isnan(xsamp[every-1]))
        std::printf("In ARMS::arms_(): NaN generated, possibly due to overflow in (log-)density (e.g. with densities involving exp(exp(...))).\n");
      if (xsamp[n-1] < minD || xsamp[every-1] > maxD)
        std::printf("In ARMS::arms_(): %d-th sample out of range [%f, %f] (fused domain). Got %f.\n", every, *xl, *xr, xsamp[every-1]);

      currentPars(j, l) = xsamp[every - 1];

      // if put 'create_mydata' out of for-loop, the following updates can for elements of pointer *mydata
      datMu.col(l) = arma::exp( datX.slice(l) *  currentPars.col(l));
      arma::vec lambdas = datMu.col(l) / std::tgamma(1. + 1./kappa);
      weibullS.col(l) = arma::exp( -arma::pow(datTime / lambdas, kappa) );

      mydata->currentPars = currentPars.memptr();
      mydata->datMu = datMu.memptr(); // update this due to its change with updated coefficients
      mydata->weibullS = weibullS.memptr(); // update this due to its change with updated coefficients

      free(xsamp);
    }
  }
  free(mydata);

  // Assembling output
  /*Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("betas") = currentPars,
    Rcpp::Named("mus") = datMu,
    Rcpp::Named("weibullS") = weibullS
  );
  */
  return currentPars;
}
