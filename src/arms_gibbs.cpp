// Gibbs sampling for multivariate ARMS

#include <cmath>

#include "arms.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


inline double upperbound1 = 70.;
inline double upperbound2 = 69.;

typedef std::vector<double> stdvec;

/*
// redefine the struct dataS, since it cannot be imported from extern "C"{} in 'arms.cpp'
typedef struct common_data
{
   //
   double *currentPars;
   unsigned int jj;
   double vA;
   double vB;
   double *datX0;
   double *datProportion;
   double *datEvent;
   double *weibullS;
}dataS;
*/

typedef struct common_data
{
   // members
   arma::vec currentPars;
   unsigned int jj;
   double *xl;
   double *xr;
   double vA;
   double vB;
   arma::mat datX0;
   arma::mat datProportion;
   arma::uvec datEvent;
   arma::mat weibullS;
}dataS;

void create_mydata(
  arma::vec currentPars,
  unsigned int jj,
  double *xl,
  double *xr,
  double vA, 
  double vB,
  arma::mat datX0,
  arma::mat datProportion,
  arma::uvec datEvent,
  arma::mat weibullS,
  dataS *abc_data) 
{
  // create objects for conditional (log-)density
  abc_data->currentPars = currentPars;
  abc_data->jj = jj;
  abc_data->xl = xl;
  abc_data->xr = xr;
  abc_data->vA = vA;
  abc_data->vB = vB;
  abc_data->datX0 = datX0;
  abc_data->datProportion = datProportion;
  abc_data->datEvent = datEvent;
  abc_data->weibullS = weibullS;

}

// The following evaluation functions might have to be in C not C++ and be put into the extern arms C code
double log_dens_xis(double par, 
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  arma::mat datX0, 
  arma::mat datProportion, 
  arma::uvec datEvent, 
  arma::mat& weibullS) {

    double logpost = 0.;

    return logpost;
}

/*
double myfunc(double x, void *mydata) {
  int a = 0;
  // TBC...
}

double myfunc(double par, 
  void *mydata) {
  double h;

  h = std::log( 4.0 * std::sqrt(2.0) * std::sqrt(par) * std::exp(-0.5 * par * (3.0 - 1.0) * (3.0 - 1.0)) / std::sqrt(M_PI) );

  return h;
}
*/

// If 'myfunc()' has to be defined in C code, all vec/mat elements in 'struct common_data{}' and 'void create_mydata()' 
//   need to be defined as 'double *', since array in C will copy big data and occupy too much memory
double myfunc(
  double par, 
  void *abc_data) 
{
  double h = 0.0; 

  // dataS * mydata_parm = *(dataS *)mydata; // error: cannot initialize a variable of type 'dataS *' with an rvalue of type 'void *'
  // *mydata_parm = &abc_data; // error: no viable overloaded '='

  // Allocation of a zero-initialized memory block of (num*size) bytes
  // 'calloc' returns a void* (generic pointer) on success of memory allocation; type-cast it due to illegal in C++ but illegal in C
  dataS *mydata_parm = (dataS *)calloc(sizeof(dataS), sizeof(dataS));
  *mydata_parm = *(dataS *)abc_data;

  //if (par >= *(mydata_parm->xl) && par <= *(mydata_parm->xr)) // this is judged in ARMS::initial()
  //{
    stdvec xis0 = arma::conv_to<stdvec>::from(mydata_parm->currentPars);
    xis0.erase(xis0.begin());
  
    arma::vec& xis = mydata_parm->currentPars;
    unsigned int jj = mydata_parm->jj;
  
    xis[jj] = par;
    double vSq = 10.;
    if (jj > 0) 
    {
      int ans = std::count(xis0.begin(), xis0.end(), 0.);
      double vA_tmp = mydata_parm->vA + 0.5 * (double)ans;
      double vB_tmp = mydata_parm->vB + 0.5 * 
        std::inner_product( xis0.begin(), xis0.end(), xis0.begin(), 0. );
      vSq = 1. / R::rgamma(vA_tmp, 1. / vB_tmp);
    }
  
    // compute the log density
    double logprior = - par * par / vSq / 2.;
  
    arma::vec eta = mydata_parm->datX0 * xis;
    eta.elem(arma::find(eta > upperbound1)).fill(upperbound1);
    arma::vec thetas = arma::exp( eta );
  
    double logpost_first = arma::accu( eta.elem(arma::find(mydata_parm->datEvent)) );
    arma::vec logpost_second = arma::zeros<arma::vec>(mydata_parm->datProportion.n_rows);
    for(unsigned int ll=0; ll<mydata_parm->datProportion.n_cols; ++ll) 
    {
      logpost_second += mydata_parm->datProportion.col(ll) % mydata_parm->weibullS.col(ll);
    }
  
    double logpost_second_sum = - arma::accu(thetas % (1. - logpost_second));
  
    h = logpost_first + logpost_second_sum + logprior;
  //} else {
  //  h = -1.0e100;
  //}
  
  return h;
}


//' Multivariate ARMS via Gibbs sampler
//'
//' @param n Number of samples to draw
//' @param every How many samples to draw for generating each sample; only the last 'n' draws will be kept; but it seems not work for every>1, not know why
//' @param ninit Number of initials as meshgrid values for envelop search
//' @param convex Adjustment for convexity (non-negative value, default 1.0)
//' @param npoint Maximum number of envelope points
//'
// [[Rcpp::export]]
arma::mat arms_gibbs(
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