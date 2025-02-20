// Gibbs sampling for multivariate ARMS

#include <cmath>

#include "arms.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


inline double upperbound1 = 70.;
inline double upperbound2 = 69.;

typedef std::vector<double> stdvec;

typedef struct common_data
{
   /* members */
   arma::vec currentPars;
   unsigned int jj;
   double vA;
   double vB;
   arma::mat datX0;
   arma::mat datProportion;
   arma::uvec datEvent;
   arma::mat weibullS;
}dataS;
/*
dataS * create_mydata (
  arma::vec currentPars,
  unsigned int jj,
  double vA, 
  double vB,
  arma::mat datX0,
  arma::mat datProportion,
  arma::uvec datEvent,
  arma::mat weibullS)
{
  dataS *tmp = (dataS *)malloc(sizeof (dataS));
  tmp->currentPars = currentPars;
  tmp->jj = jj;
  tmp->vA = vA;
  tmp->vB = vB;
  tmp->datX0 = datX0;
  tmp->datProportion = datProportion;
  tmp->datEvent = datEvent;
  tmp->weibullS = weibullS;

  return tmp;
}
*/
void create_mydata(
  arma::vec currentPars,
  unsigned int jj,
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
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
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

double myfunc(
  double par, 
  void *abc_data) 
{
  double h = 0.1; 

  //dataS * mydata_parm = *(dataS *)mydata;
  //dataS *mydata_parm = (dataS *)malloc( sizeof(dataS));
  dataS *mydata_parm = (dataS *)calloc(sizeof(dataS), sizeof(dataS));
  //*mydata_parm = &abc_data;

  stdvec xis0 = arma::conv_to<stdvec>::from(mydata_parm->currentPars);
  xis0.erase(xis0.begin());

  arma::vec& xis = mydata_parm->currentPars;
  unsigned int jj = mydata_parm->jj;

  xis[jj] = par;
  double vSq;
  if (jj == 0) 
  {
    vSq = 10.;
  } else {
    int ans = std::count(xis0.begin(), xis0.end(), 0.);
    double vA_tmp = mydata_parm->vA + 0.5 * (double)ans;
    double vB_tmp = mydata_parm->vB + 0.5 * 
      std::inner_product( xis0.begin(), xis0.end(), xis0.begin(), 0. );
    vSq = 1. / R::rgamma(vA_tmp, 1. / vB_tmp);
  }

  // compute the log density
  double logprior = - par * par / vSq / 2.;
  double logprior2 = - par / vSq;

  arma::vec eta = mydata_parm->datX0 * xis;
  eta.elem(arma::find(eta > upperbound1)).fill(upperbound1);
  arma::vec thetas = arma::exp( eta );

  double logpost_first = arma::accu( eta.elem(arma::find(mydata_parm->datEvent)) );
  double logpost_first2;
  if(jj != 0)
  {
    arma::uvec singleIdx_jj = { jj };
    logpost_first2 = arma::accu( mydata_parm->datX0.submat(arma::find(mydata_parm->datEvent), singleIdx_jj) );
  } else {
    logpost_first2 = std::count(mydata_parm->datEvent.begin(), mydata_parm->datEvent.end(), 1) + 0.;
  }

  arma::vec logpost_second = arma::zeros<arma::vec>(mydata_parm->datProportion.n_rows);
  for(unsigned int ll=0; ll<mydata_parm->datProportion.n_cols; ++ll) 
  {
    logpost_second += mydata_parm->datProportion.col(ll) % mydata_parm->weibullS.col(ll);
  }

  double logpost_second_sum = - arma::accu(thetas % (1. - logpost_second));

  h = logpost_first + logpost_second_sum + logprior;

  std::cout << "...debug: myfunc_h = " << h << "\n";
  
  return h;
}


//' Multivariate ARMS via Gibbs sampler
//'
//' @param n number of variates to draw
//' @param initialPoints this can be a matrix in multivariate case
//'
// [[Rcpp::export]]
arma::mat arms_gibbs(
  int n,
  arma::vec initialPoints,
  arma::vec minRange,
  arma::vec maxRange,

  int metropolis, 
  
  arma::vec currentPars, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat weibullS) 
{

  unsigned int p = 1;//currentPars.n_elem;

  int ninit = initialPoints.n_elem; // # here 'initialPoints' can be a vector/matrix (in univariate/multivariate cases) for initializing meshgrid values
  int dometrop = metropolis;
  // double *xl; 
  // double *xr; 
  // double *xprev; 

  //double myfunc = .1;

  // structure holding data for (log-)density
  //void *mydata = nullptr;
  //dataS *mydata0;
  //dataS *mydata0 = (dataS *)malloc(sizeof (dataS));
  //void *mydata = &mydata0;

//std::cout << "*mydata0.vA=" << mydata0->vA << "; *mydata0.datEvent" << mydata0->datEvent.t() << "\n"; 

  arma::mat samp = arma::zeros<arma::mat>(p, n);
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
      double initi = initialPoints[0]; // [j]
      double *xprev; xprev = &initi;
      double *xsamp = (double*)malloc(n * sizeof(double));

      //dataS *mydata = create_mydata(currentPars, j, vA, vB, datX0, datProportion, datEvent, weibullS);
      dataS *mydata = (dataS *)malloc(sizeof (dataS));
      create_mydata(currentPars, j, vA, vB, datX0, datProportion, datEvent, weibullS, mydata);

      double tmp = 0.0;
      int err = ARMS::arms_simple (
        ninit, xl,  xr,
        myfunc, mydata,
        dometrop, xprev, xsamp);

      samp(j, i) = tmp;
      currentPars[j] = tmp;

      free(mydata);
    }
  }

  // free memory after usage
  //free(mydata0);

  return samp;
}