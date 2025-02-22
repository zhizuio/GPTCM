/* Evaluation functions (i.e. log densities) for ARS and ARMS*/

#include "eval_func.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// log-density for coefficient xis
double log_dens_xis(
  double par, 
  void *abc_data) 
{
  // If myfunc() has to be defined in C code, all vec/mat elements in 'struct common_data{}' and 'void create_mydata()' 
  //   need to be defined as 'double *', since array in C will copy big data and occupy too much memory

  double h = 0.0; 

  // dataS * mydata_parm = *(dataS *)mydata; // error: cannot initialize a variable of type 'dataS *' with an rvalue of type 'void *'
  // *mydata_parm = &abc_data; // error: no viable overloaded '='

  // Allocation of a zero-initialized memory block of (num*size) bytes
  // 'calloc' returns a void* (generic pointer) on success of memory allocation; type-cast it due to illegal in C++ but illegal in C
  dataS *mydata_parm = (dataS *)calloc(sizeof(dataS), sizeof(dataS));
  *mydata_parm = *(dataS *)abc_data;

  //if (par >= *(mydata_parm->xl) && par <= *(mydata_parm->xr)) // this is judged in ARMS::initial()
  //{
    arma::vec xis(mydata_parm->currentPars, mydata_parm->p, true);
    //stdvec xis0 = arma::conv_to<stdvec>::from(mydata_parm->currentPars);
    stdvec xis0 = arma::conv_to<stdvec>::from(xis);
    xis0.erase(xis0.begin());

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

    //arma::mat X(mydata_parm->datX, n, p, false);  // use auxiliary memory
    //std::cout << "...X:\n" << X << "\n";
  
    //arma::vec eta = mydata_parm->datX * xis;
    arma::vec eta = arma::mat(mydata_parm->datX, mydata_parm->N, mydata_parm->p, false) * xis;
    eta.elem(arma::find(eta > upperbound)).fill(upperbound);
    arma::vec thetas = arma::exp( eta );
  
    //double logpost_first = arma::accu( eta.elem(arma::find(mydata_parm->datEvent)) );
    double logpost_first = arma::accu( eta.elem(arma::find(
      arma::ivec(mydata_parm->datEvent, mydata_parm->N, false)
    )) );
    arma::vec logpost_second = arma::zeros<arma::vec>(mydata_parm->N);
    arma::mat datProportion(mydata_parm->datProportion, mydata_parm->N, mydata_parm->L, false);
    arma::mat weibullS(mydata_parm->weibullS, mydata_parm->N, mydata_parm->L, false);
    for(unsigned int ll=0; ll<mydata_parm->L; ++ll) 
    {
      //logpost_second += mydata_parm->datProportion.col(ll) % mydata_parm->weibullS.col(ll);
      logpost_second += datProportion.col(ll) % weibullS.col(ll);
    }
  
    double logpost_second_sum = - arma::accu(thetas % (1. - logpost_second));
  
    h = logpost_first + logpost_second_sum + logprior;
  //} else {
  //  h = -1.0e100;
  //}
  
  return h;
}


// log-density for coefficient betas
double log_dens_betas(
  double par, 
  void *abc_data) 
{
  double h = 0.0; 

  // Allocation of a zero-initialized memory block of (num*size) bytes
  // 'calloc' returns a void* (generic pointer) on success of memory allocation; type-cast it due to illegal in C++ but illegal in C
  dataS *mydata_parm = (dataS *)calloc(sizeof(dataS), sizeof(dataS));
  *mydata_parm = *(dataS *)abc_data;

  //arma::vec& pars = mydata_parm->currentPars;
  arma::mat pars(mydata_parm->currentPars, mydata_parm->p, mydata_parm->L, true);
  pars(mydata_parm->jj, mydata_parm->l) = par;

  arma::mat mu_tmp(mydata_parm->datMu, mydata_parm->N, mydata_parm->L, true);
  arma::vec logMu_l = arma::mat(mydata_parm->datX, mydata_parm->N, mydata_parm->p, false) * 
    pars.col(mydata_parm->l);
  logMu_l.elem(arma::find(logMu_l > upperbound)).fill(upperbound);
  mu_tmp.col(mydata_parm->l) = arma::exp(logMu_l);

  //arma::mat weibullS_tmp = mydata_parm->weibullS;
  //arma::mat weibullS(mydata_parm->weibullS, mydata_parm->N, mydata_parm->L, false);
  //arma::mat weibullS_tmp = weibullS;
  arma::mat weibullS_tmp(mydata_parm->weibullS, mydata_parm->N, mydata_parm->L, true);
  //arma::mat lambdas = mu_tmp / std::tgamma(1. + 1./mydata_parm->kappa);
  //weibullS_tmp.col(mydata_parm->l) = arma::exp( -arma::pow(mydata_parm->datTime / lambdas.col(mydata_parm->l), mydata_parm->kappa) );
  arma::mat weibull_lambdas = mu_tmp / std::tgamma(1. + 1./mydata_parm->kappa);
  weibullS_tmp.col(mydata_parm->l) = arma::exp( - arma::pow(
    arma::vec(mydata_parm->datTime, mydata_parm->N, false) / weibull_lambdas.col(mydata_parm->l), 
    mydata_parm->kappa) );

  // compute log density
  double logprior = - par * par / mydata_parm->tauSq / 2.;

  arma::vec logpost_first = arma::zeros<arma::vec>(mydata_parm->N);
  arma::mat datProportion(mydata_parm->datProportion, mydata_parm->N, mydata_parm->L, false);
  for(int ll=0; ll<(mydata_parm->L); ++ll) 
  {
    //logpost_first += mydata_parm->datProportion.col(ll) % 
    //  arma::pow(mydata_parm->weibullS.col(ll), -mydata_parm->kappa) % weibullS_tmp.col(ll);
    logpost_first += datProportion.col(ll) % 
      arma::pow(weibull_lambdas.col(ll), - mydata_parm->kappa) % weibullS_tmp.col(ll);
  }

  //double logpost_first_sum = arma::accu( arma::log( logpost_first.elem(arma::find(mydata_parm->datEvent)) ) );
  double logpost_first_sum = arma::accu( arma::log( logpost_first.elem(arma::find(
    arma::ivec(mydata_parm->datEvent, mydata_parm->N, false)
  )) ) );

  //double logpost_second_sum = arma::accu(mydata_parm->datTheta % 
  //  mydata_parm->datProportion.col(mydata_parm->l) % weibullS_tmp.col(mydata_parm->l));
  double logpost_second_sum = arma::accu(arma::vec(mydata_parm->datTheta, mydata_parm->N, false) % 
    datProportion.col(mydata_parm->l) % weibullS_tmp.col(mydata_parm->l));

  h = logpost_first_sum + logpost_second_sum + logprior;
  
  return h;
}
