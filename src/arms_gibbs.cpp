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
  arma::uvec datEvent, 
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
  unsigned int p = currentPars.n_elem; // = 1;

  int dometrop = metropolis;
  
  // convert arma::mat to C double pointer or array
  // TBA...

  arma::mat samp = arma::zeros<arma::mat>(p, n + 1);
  samp.col(0) = currentPars;

  double minD;
  double maxD;
  for (unsigned int i = 0; i < n; ++i)
  {
    for (unsigned int j = 0; j < p; ++j)
    {
      minD = minRange[0]; // [j]
      maxD = maxRange[0]; // [j]
      double *xl; xl = &minD;
      double *xr; xr = &maxD;
      double initi = samp(j, i);  //samp(j, i)
      double *xprev; xprev = &initi;
      double *xsamp = (double*)malloc(every * sizeof(double));
      double qcent[1], xcent[1];
      int neval, ncent = 0;

      //dataS *mydata = create_mydata(currentPars, j, vA, vB, datX0, datProportion, datEvent, weibullS);
      dataS *mydata = (dataS *)malloc(sizeof (dataS));
      create_mydata(currentPars, j, xl, xr, vA, vB, datX0, datProportion, datEvent, weibullS, mydata);

      //*xprev = samp(j, i);
      int err;
      if (simple)
      {
        err = ARMS::arms_simple (
          ninit, xl,  xr,
          myfunc, mydata,
          dometrop, xprev, xsamp);
      } else {
        arma::vec xinit0 = arma::linspace( minD+1.0e-10, maxD-1.0e-10, ninit );
        double xinit[ninit];
        for (unsigned int i = 0; i < ninit; ++i)
          xinit[i] = xinit0[i];

        err = ARMS::arms (
          xinit, ninit, xl,  xr,
          myfunc, mydata,
          &convex, npoint,
          dometrop, xprev, xsamp,
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

      free(mydata);
    }
  }
  // remove the inital values in the first column of samp
  samp.shed_col(0);

  //std::cout << "...debug sample=" << samp.col(0).t() << "\n";

  return samp;
}