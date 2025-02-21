/* Evaluation functions (i.e. log densities) for ARS and ARMS*/

#include "eval_func.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]



void create_mydata(
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
  dataS *abc_data) 
{
  // create objects for conditional (log-)density
  abc_data->currentPars = currentPars;
  abc_data->jj = jj;
  abc_data->xl = xl;
  abc_data->xr = xr;
  abc_data->vA = vA;
  abc_data->vB = vB;
  abc_data->datX = datX;
  abc_data->datProportion = datProportion;
  abc_data->datEvent = datEvent;
  abc_data->weibullS = weibullS;

}

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
  
    arma::vec eta = mydata_parm->datX * xis;
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
