// Gibbs sampling for multivariate ARMS

#include <cmath>

#include "arms.h"
#include "eval_func.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


//' Multivariate ARMS via Gibbs sampler
//'
//' @param n Number of samples to draw
//' @param nsamp How many samples to draw for generating each sample; only the last draw will be kept
//' @param ninit Number of initials as meshgrid values for envelop search
//' @param convex Adjustment for convexity (non-negative value, default 1.0)
//' @param npoint Maximum number of envelope points
//'
// [[Rcpp::export]]
arma::mat arms_gibbs_xi(
  int n,
  int nsamp, 
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
  // number of parameters to be updated
  int p = currentPars.n_elem; // = 1;
  int L = datProportion.n_cols;
  int N = datProportion.n_rows;

  int dometrop = metropolis;

  double minD;
  double maxD;
  minD = minRange[0]; // [j]
  maxD = maxRange[0]; // [j]
  double *xl; xl = &minD;
  double *xr; xr = &maxD;

  double xinit[ninit];
  if (!simple)
  {
    arma::vec xinit0 = arma::linspace( minD+1.0e-10, maxD-1.0e-10, ninit );
    for (int i = 0; i < ninit; ++i)
      xinit[i] = xinit0[i];
  }

  //dataS *mydata = create_mydata(currentPars, j, vA, vB, datX0, datProportion, datEvent, weibullS);
  dataS *mydata = (dataS *)malloc(sizeof (dataS));
  //create_mydata(currentPars, j, xl, xr, vA, vB, datX0, datProportion, datEvent, weibullS, mydata);
  mydata->currentPars = currentPars.memptr();
  mydata->p = p;
  mydata->L = L;
  mydata->N = N;
  // mydata->xl = xl;
  // mydata->xr = xr;
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

  for (int i = 0; i < n; ++i)
  {
    // Gibbs sampling
    for (int j = 0; j < p; ++j)
    {
      mydata->jj = j;
      //double initi = samp(j, i);  //samp(j, i)
      //double *xprev; xprev = &initi;
      double xprev = samp(j, i);
      double *xsamp = (double*)malloc(nsamp * sizeof(double));
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
        err = ARMS::arms (
          xinit, ninit, xl,  xr,
          log_dens_xis, mydata,
          &convex, npoint,
          dometrop, &xprev, xsamp,
          nsamp, qcent, xcent, ncent, &neval);
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
          "; nsamp=" << nsamp << 
          "; qcent=" << *qcent << 
          "; xcent=" << *xcent << 
          "; ncent=" << ncent << 
          "; neval=" << neval << 
          "\n"; */
      }

      // check ARMS validity
      if (err > 0)
       std::printf("In ARMS::arms_(): error code in ARMS = %d.\n", err);
      if (isnan(xsamp[nsamp-1]))
        std::printf("In ARMS::arms_(): NaN generated, possibly due to overflow in (log-)density (e.g. with densities involving exp(exp(...))).\n");
      if (xsamp[nsamp-1] < minD || xsamp[nsamp-1] > maxD)
        std::printf("In ARMS::arms_(): %d-th sample out of range [%f, %f] (fused domain). Got %f.\n", nsamp, *xl, *xr, xsamp[nsamp-1]);

      //xprev = xsamp[nsamp - 1];
      //for (int i = 0; i < n; ++i) samp(j, i + 1) = xsamp[nsamp - n + i];
      currentPars[j] = xsamp[nsamp - 1];
      samp(j, i + 1) = xsamp[nsamp - 1];

      mydata->currentPars = currentPars.memptr();

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
//' @param nsamp How many samples to draw for generating each sample; only the last draw will be kept
//' @param ninit Number of initials as meshgrid values for envelop search
//' @param convex Adjustment for convexity (non-negative value, default 1.0)
//' @param npoint Maximum number of envelope points
//'
// [[Rcpp::export]]
arma::mat arms_gibbs_beta( 
  /* make a subfunction arms_gibbs for only vector betas that can be used for (varying-length) variable selected vector*/
  int n,
  int nsamp, 
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
  int N = datX.n_rows;
  int p = datX.n_cols; 
  int L = datX.n_slices;

  int dometrop = metropolis;

  // objects for arms()
  double minD;
  double maxD;
  minD = minRange[0]; // [j]
  maxD = maxRange[0]; // [j]
  double *xl; xl = &minD;
  double *xr; xr = &maxD;

  double xinit[ninit];
  if (!simple)
  {
    arma::vec xinit0 = arma::linspace( minD+1.0e-10, maxD-1.0e-10, ninit );
    for (int i = 0; i < ninit; ++i)
      xinit[i] = xinit0[i];
  }

  dataS *mydata = (dataS *)malloc(sizeof (dataS));
  mydata->currentPars = currentPars.memptr();
  mydata->p = p;
  mydata->L = L;
  mydata->N = N;
  // mydata->xl = xl;
  // mydata->xr = xr;
  // mydata->vA = vA;
  // mydata->vB = vB;
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
      mydata->datX = datX.slice(l).memptr();
      
      //double initi = currentPars(j, l);  //samp(j, i)
      //double *xprev; xprev = &initi;
      double xprev = currentPars(j, l);
      double *xsamp = (double*)malloc(nsamp * sizeof(double));
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
        err = ARMS::arms (
          xinit, ninit, xl,  xr,
          log_dens_betas, mydata,
          &convex, npoint,
          dometrop, &xprev, xsamp,
          nsamp, qcent, xcent, ncent, &neval);
      }

      // check ARMS validity
      if (err > 0)
       std::printf("In ARMS::arms_(): error code in ARMS = %d.\n", err);
      if (isnan(xsamp[nsamp-1]))
        std::printf("In ARMS::arms_(): NaN generated, possibly due to overflow in (log-)density (e.g. with densities involving exp(exp(...))).\n");
      if (xsamp[nsamp-1] < minD || xsamp[nsamp-1] > maxD)
        std::printf("In ARMS::arms_(): %d-th sample out of range [%f, %f] (fused domain). Got %f.\n", nsamp, *xl, *xr, xsamp[nsamp-1]);

      // std::cout << "...debug xsamp[1:nsamp]=";
      // for( int i=0; i<nsamp; ++i) {
      //   std::cout << ", " << xsamp[i] ;
      // }
      // std::cout <<  "/n";

      currentPars(j, l) = xsamp[nsamp - 1];

      // if put 'create_mydata' out of for-loop, the following updates can for elements of pointer *mydata
      arma::vec logMu_l = datX.slice(l) * currentPars.col(l);
      logMu_l.elem(arma::find(logMu_l > upperbound)).fill(upperbound);
      datMu.col(l) = arma::exp( logMu_l );
      // arma::vec lambdas = datMu.col(l) / std::tgamma(1. + 1./kappa);
      // weibullS.col(l) = arma::exp( -arma::pow(datTime / lambdas, kappa) );
      arma::vec lambdas = arma::pow( datTime / (datMu.col(l) / std::tgamma(1. + 1./kappa)), kappa);
      lambdas.elem(arma::find(lambdas > upperbound)).fill(upperbound);
      weibullS.col(l) = arma::exp( -lambdas );

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
