// Main function for the MCMC loop

#include "arms_gibbs.h"
#include "simple_gibbs.h"

// [[Rcpp::depends(RcppArmadillo)]]


//' Main function for the MCMC loop
//'
//' @param nIter Number of MCMC iterations
//' @param burnin Length of MCMC burn-in period
//' @param n Number of samples to draw
//' @param nsamp How many samples to draw for generating each sample; only the last draw will be kept
//' @param ninit Number of initials as meshgrid values for envelop search
//' @param convex Adjustment for convexity (non-negative value, default 1.0)
//' @param npoint Maximum number of envelope points
//' @param dirichlet Not yet implemented
//'
// [[Rcpp::export]]
Rcpp::List mcmc(
  int nIter,
  int burnin,

  int n,
  int nsamp, 
  int ninit,
  int metropolis, 
  bool simple,
  double convex,
  int npoint,
  bool dirichlet,
  bool proportion_model,
  
  const Rcpp::List initList, 
  const Rcpp::List rangeList,
  const Rcpp::List hyperparList, 

  arma::ivec datEvent, 
  arma::vec datTime,
  arma::cube datX, 
  arma::cube datX0, 
  arma::mat datProportionConst) 
{
  // It might be more efficient for all input data with data type &
  // const Rcpp::List& initList, 
  // const Rcpp::List& rangeList,
  // const Rcpp::List& hyperpar, 

  // const Rcpp::List& survObj, 
  // const arma::cube& datX, 
  // const arma::cube& datX0, 
  // const arma::mat& datProportionConst


  // dimensions
  int N = datX.n_rows;
  int p = datX.n_cols; 
  int L = datX.n_slices;

  // hyperparameters
  hyperparS *hyperpar = (hyperparS *)malloc(sizeof (hyperparS));
  hyperpar->vSq = Rcpp::as<double>(hyperpar["vSq"]);;
  hyperpar->vA = Rcpp::as<double>(hyperpar["vA"]);;
  hyperpar->vB = Rcpp::as<double>(hyperpar["vB"]);;
  hyperpar->tauSq = Rcpp::as<double>(hyperpar["tauSq"]);;
  hyperpar->tauA = Rcpp::as<double>(hyperpar["tauA"]);;
  hyperpar->tauB = Rcpp::as<double>(hyperpar["tauB"]);;
  hyperpar->wSq = Rcpp::as<double>(hyperpar["wSq"]);;
  hyperpar->wA = Rcpp::as<double>(hyperpar["wA"]);;
  hyperpar->wB = Rcpp::as<double>(hyperpar["wB"]);;
  hyperpar->w0Sq = Rcpp::as<double>(hyperpar["w0Sq"]);;
  hyperpar->w0A = Rcpp::as<double>(hyperpar["w0A"]);;
  hyperpar->w0B = Rcpp::as<double>(hyperpar["w0B"]);;
  hyperpar->kappaA = Rcpp::as<double>(hyperpar["kappaA"]);;
  hyperpar->kappaB = Rcpp::as<double>(hyperpar["kappaB"]);;

  // ranges of parameters
  rangeS *range = (rangeS *)malloc(sizeof (rangeS));
  range->xiMin = Rcpp::as<double>(hyperpar["xiMin"]);;
  range->xiMax = Rcpp::as<double>(hyperpar["xiMax"]);;
  range->zetaMin = Rcpp::as<double>(hyperpar["zetaMin"]);;
  range->zetaMax = Rcpp::as<double>(hyperpar["zetaMax"]);;
  range->kappaMin = Rcpp::as<double>(hyperpar["kappaMin"]);;
  range->kappaMax = Rcpp::as<double>(hyperpar["kappaMax"]);;
  range->betaMin = Rcpp::as<double>(hyperpar["betaMin"]);;
  range->betaMax = Rcpp::as<double>(hyperpar["betaMax"]);;

  // initial values of parameters
  arma::vec xi = Rcpp::as<arma::vec>(initList["xi"]);;
  arma::mat zetas = Rcpp::as<arma::mat>(initList["zetas"]);;
  arma::mat betas = Rcpp::as<arma::mat>(initList["betas"]);;
  double kappas = Rcpp::as<double>(initList["kappas"]);;

  // initializing mcmc results
  arma::vec vSq_mcmc = arma::zeros<arma::vec>(1+m);
  arma::mat xi_mcmc = arma::zeros<arma::mat>(1+m, xi.n_elem);
  arma::vec wSq_mcmc = arma::zeros<arma::vec>(1+m);
  arma::mat zeta_mcmc = arma::zeros<arma::mat>(1+m, p*L);
  arma::vec tauSq_mcmc = arma::zeros<arma::vec>(1+m);
  arma::mat beta_mcmc = arma::zeros<arma::mat>(1+m, p*L);
  arma::vec kappa_mcmc = arma::zeros<arma::vec>(1+m);

  // initializing relevant quantities; can be declared like arma::mat&
  arma::vec datTheta; 
  arma::vec datMu; 
  arma::mat datProportion;
  arma::mat weibullS;
  arma::mat weibullLambda;
  
  // main MCMC loop
  for (int m=0; m<nIter; ++m)
  {

    // update \xi's variance vSq
    hyperpar->vSq = sampleV(1, hyperpar->vA, hyperpar->vB, xi);
    vSq_mcmc[1+m] = hyperpar->vSq;

    // update \xi's in cure fraction
    xi = arms_gibbs_xi
    (
      n,
      nsamp, 
      ninit,
      range->xiMin,
      range->xiMax,
    
      metropolis, 
      simple,
      convex,
      npoint,
      
      xi, 
      hyperpar->vA, 
      hyperpar->vB, 
      datX0, 
      datProportion, 
      datEvent, 
      weibullS
    );
    xi_mcmc.row(1+m) = xi;

    // update cure rate based on the new xi
    thetas = arma::exp(datX0 * xi);

    // update parameters in the proportion model
    if(proportion_model)
    {
      if(dirichlet)
      {
        // update \zetas' variance wSq
        hyperpar->wSq = sampleW(1, hyperpar->wA, hyperpar->wB, zetas.rows(1, p));
        wSq_mcmc[1+m] = hyperpar->wSq;
        if(w0IGamma)
        {
          hyperpar->w0Sq = sampleW(1, hyperpar->w0A, hyperpar->w0B, zetas.row(0));
        }

        zetas = arms_gibbs_zeta
        ( 
          int n,
          int nsamp, 
          int ninit,
          range->zetaMin,
          range->zetaMax,
        
          int metropolis, 
          bool simple,
          double convex,
          int npoint,
          
          zetas, 
          hyperpar->w0Sq, 
          hyperpar->wSq, 
          kappas,
          dirichlet,
          datX, 
          datTheta,
          datProportionConst,
          datEvent, 
          weibullS, 
          weibullLambda
        );
        zeta_mcmc.row(1+m) = arma::vectorise(zetas);

        // update Dirichlet's concentrations and proportions based on the new zetas
        arma::mat alphas = arma::zeros<arma::mat>(N, L);
        for(int ll=0; ll<L; ++ll) 
        {
          alphas.col(ll) = arma::exp( currentPars(0, ll) + datX.slice(ll) * currentPars.submat(1, ll, p, ll) );
        }
        alphas.elem(arma::find(alphas > upperbound3)).fill(upperbound3);
        alphas.elem(arma::find(alphas < lowerbound)).fill(lowerbound);
        datProportion = alphas / arma::repmat(arma::sum(alphas, 1), 1, L);

      } else {
        std::printf("Warning: In arms_gibbs_zeta(), Dirichlet modeling with logit/alr-link is not implement!\n");
        break;
      }
    }

    // update Weibull's shape parameter kappa
    kappas = arms_kappa
    ( 
      n,
      nsamp, 
      ninit,
      range->kappaMin,
      range->kappaMax,
    
      metropolis, 
      simple,
      convex,
      npoint,
      
      kappas, 
      hyperpar->kappaA,
      hyperpar->kappaB,
      kappaIGamma,
      datTheta,
      datMu,
      datProportion, 
      datEvent, 
      datTime
    ); 
    kappa_mcmc[1+m] <- kappas

    // update Weibull's quantities based on the new kappas
    for(int l=0; l<L; ++l) 
    {
      weibullLambda.col(l) = datMu.col(l) / std::tgamma(1.0+1.0/kappas);
      weibullS.col(l) = arma::exp(- arma::pow( datTime/weibullLambda.col(l), kappas));
    }

    // update \betas' variance tauSq
    hyperpar->tauSq = sampleTau(1, hyperpar->tauA, hyperpar->tauB, zetas.rows(1, p));
    tauSq[1+m] = hyperpar->tauSq;

    // update \betas in non-cure fraction
    betas = arms_gibbs_beta
    (
      n,
      nsamp, 
      ninit,
      range->betaMin,
      range->betaMax,
    
      metropolis, 
      simple,
      convex,
      npoint,
      
      betas, 
      hyperpar->tauSq,
      kappas,
      datX, 
      datTheta,
      datMu,
      datProportion, 
      datEvent, 
      datTime,
      weibullS
    );
    beta_mcmc.row(1+m) = arma::vectorise(betas);

    // update Weibull's quantities based on the new betas
    for(int l=0; l<L; ++l) 
    {
      arma::vec logMu_l = datX.slice(l) * betas.col(l);
      logMu_l.elem(arma::find(logMu_l > upperbound)).fill(upperbound);
      datMu.col(l) = arma::exp( logMu_l );
      weibullLambda.col(l) = datMu.col(l) / std::tgamma(1.0+1.0/kappas);
      weibullS.col(l) = arma::exp(- arma::pow( datTime/weibullLambda.col(l), kappas));
    }

  }
  
  // wrap all outputs
  Rcpp::List results;
  results["samples"] = variates;
  results["T"] = Rcpp::NumericVector(state.T.begin(), state.T.end());
  results["k"] = state.k;
  results["Z"] = state.Z;
  results["H"] = state.H;
  results["Hprime"] = state.Hprime;
  results["Q"] = state.Q;

  return results;
}
