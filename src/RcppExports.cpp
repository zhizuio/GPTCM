// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// arms_gibbs_xi
arma::mat arms_gibbs_xi(int n, int every, int ninit, arma::vec minRange, arma::vec maxRange, int metropolis, bool simple, double convex, int npoint, arma::vec currentPars, double vA, double vB, arma::mat datX0, arma::mat datProportion, arma::ivec datEvent, arma::mat weibullS);
RcppExport SEXP _GPTCM_arms_gibbs_xi(SEXP nSEXP, SEXP everySEXP, SEXP ninitSEXP, SEXP minRangeSEXP, SEXP maxRangeSEXP, SEXP metropolisSEXP, SEXP simpleSEXP, SEXP convexSEXP, SEXP npointSEXP, SEXP currentParsSEXP, SEXP vASEXP, SEXP vBSEXP, SEXP datX0SEXP, SEXP datProportionSEXP, SEXP datEventSEXP, SEXP weibullSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type every(everySEXP);
    Rcpp::traits::input_parameter< int >::type ninit(ninitSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type minRange(minRangeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type maxRange(maxRangeSEXP);
    Rcpp::traits::input_parameter< int >::type metropolis(metropolisSEXP);
    Rcpp::traits::input_parameter< bool >::type simple(simpleSEXP);
    Rcpp::traits::input_parameter< double >::type convex(convexSEXP);
    Rcpp::traits::input_parameter< int >::type npoint(npointSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type currentPars(currentParsSEXP);
    Rcpp::traits::input_parameter< double >::type vA(vASEXP);
    Rcpp::traits::input_parameter< double >::type vB(vBSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type datX0(datX0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type datProportion(datProportionSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type datEvent(datEventSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type weibullS(weibullSSEXP);
    rcpp_result_gen = Rcpp::wrap(arms_gibbs_xi(n, every, ninit, minRange, maxRange, metropolis, simple, convex, npoint, currentPars, vA, vB, datX0, datProportion, datEvent, weibullS));
    return rcpp_result_gen;
END_RCPP
}
// ars
arma::vec ars(int n, arma::vec initialPoints, double minRange, double maxRange, unsigned int jj, arma::vec& xis, double vA, double vB, const arma::mat datX0, const arma::mat datProportion, const arma::uvec datEvent, arma::mat& weibullS);
RcppExport SEXP _GPTCM_ars(SEXP nSEXP, SEXP initialPointsSEXP, SEXP minRangeSEXP, SEXP maxRangeSEXP, SEXP jjSEXP, SEXP xisSEXP, SEXP vASEXP, SEXP vBSEXP, SEXP datX0SEXP, SEXP datProportionSEXP, SEXP datEventSEXP, SEXP weibullSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initialPoints(initialPointsSEXP);
    Rcpp::traits::input_parameter< double >::type minRange(minRangeSEXP);
    Rcpp::traits::input_parameter< double >::type maxRange(maxRangeSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type xis(xisSEXP);
    Rcpp::traits::input_parameter< double >::type vA(vASEXP);
    Rcpp::traits::input_parameter< double >::type vB(vBSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type datX0(datX0SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type datProportion(datProportionSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type datEvent(datEventSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type weibullS(weibullSSEXP);
    rcpp_result_gen = Rcpp::wrap(ars(n, initialPoints, minRange, maxRange, jj, xis, vA, vB, datX0, datProportion, datEvent, weibullS));
    return rcpp_result_gen;
END_RCPP
}
// ars_debug
Rcpp::List ars_debug(int n, arma::vec initialPoints, double minRange, double maxRange, long maxiter, unsigned int jj, arma::vec& xis, double vA, double vB, const arma::mat datX0, const arma::mat datProportion, const arma::uvec datEvent, arma::mat& weibullS);
RcppExport SEXP _GPTCM_ars_debug(SEXP nSEXP, SEXP initialPointsSEXP, SEXP minRangeSEXP, SEXP maxRangeSEXP, SEXP maxiterSEXP, SEXP jjSEXP, SEXP xisSEXP, SEXP vASEXP, SEXP vBSEXP, SEXP datX0SEXP, SEXP datProportionSEXP, SEXP datEventSEXP, SEXP weibullSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initialPoints(initialPointsSEXP);
    Rcpp::traits::input_parameter< double >::type minRange(minRangeSEXP);
    Rcpp::traits::input_parameter< double >::type maxRange(maxRangeSEXP);
    Rcpp::traits::input_parameter< long >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type xis(xisSEXP);
    Rcpp::traits::input_parameter< double >::type vA(vASEXP);
    Rcpp::traits::input_parameter< double >::type vB(vBSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type datX0(datX0SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type datProportion(datProportionSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type datEvent(datEventSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type weibullS(weibullSSEXP);
    rcpp_result_gen = Rcpp::wrap(ars_debug(n, initialPoints, minRange, maxRange, maxiter, jj, xis, vA, vB, datX0, datProportion, datEvent, weibullS));
    return rcpp_result_gen;
END_RCPP
}
// ars_gibbs
arma::mat ars_gibbs(int n, arma::vec initialPoints, arma::vec minRange, arma::vec maxRange, arma::vec& par, double vA, double vB, const arma::mat datX0, const arma::mat datProportion, const arma::uvec datEvent, arma::mat& weibullS);
RcppExport SEXP _GPTCM_ars_gibbs(SEXP nSEXP, SEXP initialPointsSEXP, SEXP minRangeSEXP, SEXP maxRangeSEXP, SEXP parSEXP, SEXP vASEXP, SEXP vBSEXP, SEXP datX0SEXP, SEXP datProportionSEXP, SEXP datEventSEXP, SEXP weibullSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initialPoints(initialPointsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type minRange(minRangeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type maxRange(maxRangeSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type par(parSEXP);
    Rcpp::traits::input_parameter< double >::type vA(vASEXP);
    Rcpp::traits::input_parameter< double >::type vB(vBSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type datX0(datX0SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type datProportion(datProportionSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type datEvent(datEventSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type weibullS(weibullSSEXP);
    rcpp_result_gen = Rcpp::wrap(ars_gibbs(n, initialPoints, minRange, maxRange, par, vA, vB, datX0, datProportion, datEvent, weibullS));
    return rcpp_result_gen;
END_RCPP
}
// sampleGamma
Rcpp::List sampleGamma(arma::umat Gammas, bool gamma_sampler_bandit, double a_pi, double b_pi);
RcppExport SEXP _GPTCM_sampleGamma(SEXP GammasSEXP, SEXP gamma_sampler_banditSEXP, SEXP a_piSEXP, SEXP b_piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type Gammas(GammasSEXP);
    Rcpp::traits::input_parameter< bool >::type gamma_sampler_bandit(gamma_sampler_banditSEXP);
    Rcpp::traits::input_parameter< double >::type a_pi(a_piSEXP);
    Rcpp::traits::input_parameter< double >::type b_pi(b_piSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleGamma(Gammas, gamma_sampler_bandit, a_pi, b_pi));
    return rcpp_result_gen;
END_RCPP
}
