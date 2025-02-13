// Adaptive rejection sampling algorithm for R based on Gilks & Wild (1992).
// Allows efficient sampling from univariate densities that are log concave.
//
// Adapted from https://github.com/kuperov/adsample/blob/master/src/adsample.cpp (Alex Cooper, 2017)

#include <cmath>
#include <cstdlib>  // for abs()
#include <vector>
#include <utility>

using namespace std;

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

inline double upperbound1 = 70.;
inline double upperbound2 = 69.;

// (point value, gradient value)
//typedef pair<double, double> pointgrad;

// callback to evaluate the (possibly expensive) log density function,
// which returns the ordinate and the gradient evaluated at the abscissa x.
// Note that while the gradient need not be normalized (it can differ from
// the true gradient up to a constant multiple), the gradient and ordinate
// must share the *same* constant.
//typedef std::function<pointgrad(double x)> log_dens_callback;

// container for the state of the algorithm, for returning state back to
// R for debugging
typedef struct {
  int k;
  vector<double> T;
  vector<double> H;
  vector<double> Hprime;
  vector<double> Z;
  vector<double> Q;
} algo_state;

typedef std::vector<double> stdvec;

arma::vec log_dens_xi0(double par, 
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) {
  // a 2-vector with un-normalized log density and its first derivative
  arma::vec h = arma::zeros<arma::vec>(2); 

  h[0] = std::log( 4.0 * std::sqrt(2.0) * std::sqrt(par) * std::exp(-0.5 * par * (3.0 - 1.0) * (3.0 - 1.0)) / std::sqrt(M_PI) );

  h[1] = 0.5 * (-par * (3.0 - 1.0) * (3.0 - 1.0) + 1.0) / par;

  return h;
}

// a 2-vector with un-normalized log density and its first derivative
arma::vec log_dens_xi(
  double par, 
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) 
{
  arma::vec h = arma::zeros<arma::vec>(2); 

  //arma::vec xis0 = xis;
  //xis0.shed_row(0);
  stdvec xis0 = arma::conv_to<stdvec>::from(xis);
  xis0.erase(xis0.begin());

  xis[jj] = par;
  double vSq;
  if (jj == 0) 
  {
    vSq = 10.;
  } else {
    int ans = std::count(xis0.begin(), xis0.end(), 0.);
    double vA_tmp = vA + 0.5 * (double)ans;
    //arma::uvec tmp = arma::find(xis0);
    //double vA_tmp = vA + 0.5 * (double)tmp.n_elem;//arma::as_scalar(tmp.n_elem);
    //double vB_tmp = vB + 0.5 * arma::accu(xis0 % xis0);
    double vB_tmp = vB + 0.5 * 
      std::inner_product( xis0.begin(), xis0.end(), xis0.begin(), 0. );
    vSq = 1. / R::rgamma(vA_tmp, 1. / vB_tmp);
  }

  // compute the log density
  double logprior = - par * par / vSq / 2.;
  double logprior2 = - par / vSq;

  arma::vec eta = datX0 * xis;
  eta.elem(arma::find(eta > upperbound1)).fill(upperbound1);
  arma::vec thetas = arma::exp( eta );
  
  //double logpost_first = arma::accu( arma::log(thetas.elem(arma::find(datEvent == 1))) );
  double logpost_first = arma::accu( eta.elem(arma::find(datEvent)) );
  double logpost_first2;
  if(jj != 0)
  {
    arma::uvec singleIdx_jj = { jj };
    logpost_first2 = arma::accu( datX0.submat(arma::find(datEvent), singleIdx_jj) );
  } else {
    //logpost_first2 = arma::accu(datEvent);
    logpost_first2 = std::count(datEvent.begin(), datEvent.end(), 1) + 0.;
  }

  arma::vec logpost_second = arma::zeros<arma::vec>(datProportion.n_rows);
  for(unsigned int ll=0; ll<datProportion.n_cols; ++ll) 
  {
    logpost_second += datProportion.col(ll) % weibullS.col(ll);
  }

  double logpost_second_sum = - arma::accu(thetas % (1. - logpost_second));

  //h[0] = logpost_first - arma::accu(thetas % (1. - logpost_second)) + logprior;
  h[0] = logpost_first + logpost_second_sum + logprior;

  // compute the first derivative of the log density
  h[1] = logpost_first2 + logpost_second_sum + logprior2;

  return h;
}

// Draw n variates from the unnormalized density h using adaptive sampling.
//
// IMPORTANT: the density MUST be univariate and log concave in this parameter.
//
// params:
//    n - number of variates to draw
//    state - optional pointer to a map for storing algorithm state. Leave as
//            NULL to discard state at end of algorithm
//    maxiter - maximum number of iterations, zero for no limit
arma::vec ars_internal(
  int n,
  std::vector<double> initAbscissae,
  double minD,
  double maxD,
  algo_state* state,// = NULL,
  long maxiter,// = 0, 

  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) 
{
  // vector<double> variates(n);
  arma::vec variates = arma::zeros<arma::vec>(n); 

  int k = initAbscissae.size();

  vector<double> T(initAbscissae); // [x_1, x_2, ..., x_k]
  vector<double> H(k);             // [h(x_1), h(x_2), ..., h(x_k)]
  vector<double> Hprime(k);        // [h'(x_1), h'(x_2), ..., h'(x_k)]
  for (int i=0; i<k; i++) {
    //auto pgrad = h(T[i]);
    //auto pgrad = log_dens(T[i]);

    auto pgrad = log_dens_xi(
      T[i], 
      jj,
      xis,
      vA,
      vB,
      datX0,
      datProportion,
      datEvent,
      weibullS
    );
    //std::cout << "...debug  pgrad=" << pgrad << "\n";

    H[i] = pgrad[0];
    Hprime[i] = pgrad[1];
  }

  // Compute the value of z_i. This is equation (1) in the paper,
  // and unlike the other variables in this program the subscript used
  // in the paper lines up exactly with the i parameter used here.
  auto z = [&T, &H, &Hprime](int i) {
    return (H[i] - H[i-1] - T[i]*Hprime[i] + T[i-1]*Hprime[i-1])/
      (Hprime[i-1] - Hprime[i]);
  };

  // partition of piecewise envelope function
  vector<double> Z(k + 1);
  Z[0] = minD;
  Z[k] = maxD;
  for (int i=1; i < k; i++) {
    Z[i] = z(i);
  }

  // Compute the area under the unnormalized envelope curve Q*s_k(x) between
  // z_i and z_{i+1}, ie the (i+1)th piecewise segment of Q*s_k(k).
  //
  // Note this is not the normalized s_k(x) given in equation (3), but instead
  //     Q*s_k(x) = exp(u(x)) = exp(h(x_1) + (x - x_1)*h'(x)),
  // ie Q has the value of the integral over D following the / in equation (3).
  //
  auto env_area_pt = [&T, &H, &Hprime, &Z](int i) {
    //return std::exp(H[i] - T[i]*Hprime[i])* (std::exp(Z[i+1]*Hprime[i]) - std::exp(Z[i]*Hprime[i]))/Hprime[i];
    
    double e1 = H[i] - T[i] * Hprime[i];
    double e21 = Z[i+1] * Hprime[i];
    double e22 = Z[i] * Hprime[i];
    
    double e1_trunc = e1 > upperbound1 ? upperbound1 : e1; // truncation by considering e1 > e2 or vice versa
    //double e21_trunc = e21 > upperbound1 ? upperbound1 : e21; 
    double e21_trunc = e21 < upperbound1 ? e21 : (e21>e22 ? upperbound1 : upperbound2); 
    double e22_trunc = e22 < upperbound1 ? e22 : (e22>e21 ? upperbound1 : upperbound2); 

    return std::exp(e1_trunc) * (std::exp(e21_trunc) - std::exp(e22_trunc)) / Hprime[i];
    
  };

  // Qtot and Q are area under s_k(x) partitioned by Zs
  vector<double> Q(k);
  double Qtot = 0.;
  for (int i=0; i<k; i++) {
    double area = env_area_pt(i);
    Q[i] = area;
    Qtot += area;
  }

  int nsampled = 0;
  for (long iter = 0; nsampled < n && (maxiter == 0 || iter >= maxiter);
       iter++) {
    // step 1: draw from s_k(x)
    //double w = as<double>(runif(1));
    double w = R::runif(0., 1.);
    // renormalize this draw to the measure of envelope s_k(x)
    double wq = Qtot * w;
    int i = 0;  // index for piecewise chunk of s_k(x)
    double q_partial = 0; // area to the left of ith chunk
    while (q_partial + Q[i] < wq)
      q_partial += Q[i++];
    
      // xstar is the x ordinate corresponding to w
    //double xstar = std::log((wq - q_partial) * Hprime[i] * std::exp(T[i]*Hprime[i] - H[i]) + std::exp(Z[i] * Hprime[i])) / Hprime[i];
    
      double e1 = T[i] * Hprime[i] - H[i];
      double e2 = Z[i] * Hprime[i];
      double e1_trunc = e1 < upperbound1 ? e1 : upperbound1; // truncation by considering e1 > e2 or vice versa
      double e2_trunc = e2 < upperbound1 ? e2 : upperbound1; 
  
      //double tmp = std::log((wq - q_partial) * Hprime[i] * std::exp(e1_trunc) + std::exp(e2_trunc)) / Hprime[i];
      // the following condition is needed due to the numerical truncation sometimes
      //double xstar = (tmp > minD && tmp < maxD) ? tmp : ( (tmp > maxD) ? maxD : minD );
      //double xstar = tmp > maxD ? maxD : ( tmp < minD ? minD : tmp );
    double xstar = std::log((wq - q_partial) * Hprime[i] * std::exp(e1_trunc) + std::exp(e2_trunc)) / Hprime[i];

    
/*
    std::cout << "...debug  xstar=" << xstar << 
    "; wq=" << wq << 
    "; q_partial=" << q_partial << 
    "; i=" << i << 
    "; Hprime[i]=" << Hprime[i] << 
    "; T[i]=" << T[i] << 
    "; Hprime[i]=" << Hprime[i] << 
    ";  H[i]=" <<  H[i] << 
    "; Z[i]=" << Z[i] << 
    "\n";
    if(std::isnan(xstar)) break; */

    // step 2: squeezing test (sec 2.2.2) using cheap upper and lower hulls
    //w = as<double>(Rcpp::runif(1));
    w = R::runif(0., 1.);
    double u_xstar = H[i] + (xstar - T[i])*Hprime[i];  // equation (2)
    double l_xstar = ((i > 0) && (i < k-1)) ?
    ((T[i+1] - xstar)*H[i] + (xstar - T[i])*H[i+1])/(T[i+1] - T[i])  // eq (4)
      :
      -std::numeric_limits<double>::infinity();   // outside support of l(x)
    if (w <= std::exp(l_xstar - u_xstar)) {
      variates[nsampled++] = xstar;
      continue;
    }

    // squeezing test has failed, so we now pay the cost of evaluating h(x)
    // and h'(x), and will add an additional abscissa to T below
    //auto h_xstar_pair = h(xstar);
    //auto h_xstar_pair = log_dens(xstar);

    auto h_xstar_pair = log_dens_xi(
      xstar, 
      jj,
      xis,
      vA,
      vB,
      datX0,
      datProportion,
      datEvent,
      weibullS
    );
    double h_xstar = h_xstar_pair[0];
    double h_xstar_prime = h_xstar_pair[1];
    // rejection test (sec 2.2.2 contd)
    if (w <= std::exp(h_xstar - u_xstar)) {
      variates[nsampled++] = xstar;
      // no continue here since we still need to save h(xstar) and h'(xstar)
    }

    // step 3: save h(xstar) and h'(xstar) and recompute the others
    k += 1;
    if (xstar < T[i]) {
      // xstar inserted to LEFT of x_j
      int j = i;

      T.insert(T.begin()+j, xstar);
      H.insert(H.begin()+j, h_xstar);
      Hprime.insert(Hprime.begin()+j, h_xstar_prime);
      Z.insert(Z.begin()+(j+1), z(j+1));
      if (j != 0) { // LH boundary is never updated
        Z[j] = z(j);
      }
      Qtot -= Q[i];
      Q.insert(Q.begin()+j, env_area_pt(j));
      // we redraw the boundary of region to right
      Q[j + 1] = env_area_pt(j+1);
      Qtot += Q[j + 1] + Q[j];
      if (j != 0) {
        Qtot -= Q[j-1];
        // we redraw the boundary of region to left
        Q[j-1] = env_area_pt(j-1);
        Qtot += Q[j-1];
      }
    } else {
      // xstar inserted to RIGHT of x_j
      int j = i + 1;
      /* printf("RHS j = %d\n", j); */
      T.insert(T.begin()+j, xstar);
      H.insert(H.begin()+j, h_xstar);
      Hprime.insert(Hprime.begin()+j, h_xstar_prime);
      Z.insert(Z.begin()+j, z(j));
      if (j + 1 < k) {  // RH boundary is never updated
        Z[j+1] = z(j+1);
      }
      // Z[j] = z(j, T, H, Hprime);
      Qtot -= Q[i];
      Q.insert(Q.begin()+j, env_area_pt(j));
      // we redraw the boundary of region to left
      Q[j-1] = env_area_pt(j-1);
      Qtot += Q[j] + Q[j-1];
      if (j < k-1) {
        Qtot -= Q[j+1];
        // we redraw the boundary of region to right
        Q[j+1] = env_area_pt(j+1);
        Qtot += Q[j+1];
      }
    }
  }

  if (state != NULL) {
    state->k = k;
    state->T = T;
    state->H = H;
    state->Hprime = Hprime;
    state->Z = Z;
    state->Q = Q;
  }

  return variates;
}

//' Rcpp wrapper for adsamp algorithm
//'
//' @param n number of variates to draw
//' @param log_dens log density function, a function with a single parameter
//' @param initialPoints at least 2 points in support to seed algorithm
//' @param minRange lower bound of support
//' @param maxRange upper bound of support
//'
// [[Rcpp::export]]
arma::vec ars(
  int n,
  arma::vec initialPoints,
  double minRange,
  double maxRange,
  
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) 
{

  // lambda for getting log density and its slope
  /*
  log_dens_callback h = [log_dens](double x) {
    NumericVector pointslope = log_dens(x);
    double y = pointslope[0];
    double yprime = pointslope[1];
    return std::pair<double, double>(y, yprime);
  };
  */
  // initial abscissae
  std::vector<double> initT(initialPoints.begin(), initialPoints.end());
  // boundaries of the support of h
  //double minD = minRange;
  //double maxD = maxRange;

  long maxiter = 0;
  //auto samp = ars_internal(n, h, initT, minD, maxD, NULL, 0);
  auto samp = ars_internal(n, initT, minRange, maxRange, NULL, maxiter, 
                            jj,
                            xis, 
                            vA, 
                            vB, 
                            datX0, 
                            datProportion, 
                            datEvent, 
                            weibullS);
  return Rcpp::NumericVector(samp.begin(), samp.end());
}

//' Duplicate wrapper for adsamp algorithm for debugging. This version
//' returns the algorithm state as well as the variates.
//'
//' @param n number of variates to draw
//' @param log_dens log density function, a function with a single parameter
//' @param initialPoints at least 2 points in support to seed algorithm
//' @param minRange lower bound of support
//' @param maxRange upper bound of support
//' @param maxiter maximum number of iterations; zero if no limit
//'
// [[Rcpp::export]]
Rcpp::List ars_debug(
  int n,
  arma::vec initialPoints,
  double minRange,
  double maxRange,
  long maxiter,
  
  unsigned int jj,
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) 
{

  // lambda for getting log density and its slope
  /*
  log_dens_callback h = [log_dens](double x) {
    NumericVector pointslope = log_dens(x);
    double y = pointslope[0];
    double yprime = pointslope[1];
    return std::pair<double, double>(y, yprime);
  };
  */
  // initial abscissae
  std::vector<double> initT(initialPoints.begin(), initialPoints.end());
  // boundaries of the support of h
  //double minD = minRange;
  //double maxD = maxRange;

  algo_state state;
  //auto samp = ars_internal(n, h, initT, minD, maxD, &state, maxiter);
  auto samp = ars_internal(n, initT, minRange, maxRange, &state, maxiter, 
                            jj,
                            xis, 
                            vA, 
                            vB, 
                            datX0, 
                            datProportion, 
                            datEvent, 
                            weibullS);
  auto variates = Rcpp::NumericVector(samp.begin(), samp.end());

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

//' Multivariate ARS via Gibbs sampler
//'
//' @param n number of variates to draw
//' @param initialPoints this can be a matrix in multivariate case
//'
// [[Rcpp::export]]
arma::mat arsGibbs(
  int n,
  arma::vec initialPoints,
  arma::vec minRange,
  arma::vec maxRange,
  
  arma::vec& xis, 
  double vA, 
  double vB, 
  const arma::mat datX0, 
  const arma::mat datProportion, 
  const arma::uvec datEvent, 
  arma::mat& weibullS) 
{

  unsigned int p = xis.n_elem;

  // initial abscissae
  /*arma::vec initT = initialPoints;

  if (initialPoints.n_elem == 1) 
  {
    initT.insert_rows(1, initialPoints[0]*arma::ones<arma::vec>(p-1));
  }*/
  
  // boundaries of the support of h
  /*
  arma::vec minD = minRange;
  arma::vec maxD = maxRange;
  
  if (minRange.n_elem == 1)
  {
    //arma::vec minD(p, arma::fill::value(minRange));
    minD.insert_rows(1, minRange[0]*arma::ones<arma::vec>(p-1));
  }
  if (maxRange.n_elem == 1)
  {
    //arma::vec maxD(p, arma::fill::value(maxRange));
    maxD.insert_rows(1, maxRange[0]*arma::ones<arma::vec>(p-1));
  } 
  
  std::cout << "...Debug minD=" << minD.t() << "\n" <<
  "...maxD=" << maxD.t() << "\n";
  */

  arma::mat samp = arma::zeros<arma::mat>(n, p);
  for (unsigned int i = 0; i < n; ++i)
  {
    for (unsigned int j = 0; j < p; ++j)
    {
      // # here 'initialPoints' can be a vector/matrix for initializing meshgrid values
      //double minD = minRange[0]; double maxD = maxRange[0];
      arma::vec tmp = ars(1, initialPoints, minRange[0], maxRange[0],
        j,
        xis, 
        vA, 
        vB, 
        datX0, 
        datProportion, 
        datEvent, 
        weibullS);
        /*arma::vec tmp = ars(1, initialPoints, maxD[j], maxD[j],
          j,
          xis, 
          vA, 
          vB, 
          datX0, 
          datProportion, 
          datEvent, 
          weibullS);*/

          samp(i, j) = tmp[0];
          xis[j] = tmp[0];
    }
  }

  return samp;
}