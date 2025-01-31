#ifndef sampleGamma_H
#define sampleGamma_H

#include <RcppArmadillo.h>

// functions
Rcpp::List sampleGamma(arma::umat ,
                      bool ,
                      double ,
                       double );
arma::umat gammaMC3Proposal(arma::umat& ,
                            arma::uvec& ,
                            const int );
arma::umat gammaBanditProposal( 
    arma::umat& , 
    arma::uvec& , 
    int ); // steppedGamma , updateIdx , outcomeIdx
                           
arma::uvec randWeightedIndexSampleWithoutReplacement
(
    unsigned int ,
    const arma::vec& ,
    unsigned int
);
arma::uvec randWeightedIndexSampleWithoutReplacement
(
    unsigned int ,
    unsigned int
);
arma::uword randWeightedIndexSampleWithoutReplacement
(
    unsigned int ,
    const arma::vec&
);
double logPDFWeightedIndexSampleWithoutReplacement(const arma::vec& ,
                                                   const arma::uvec& );
double logspace_add(double ,
                    double );
                 
double logPDFBernoulli( int ,
                       double );
double logPDFBernoulli(const arma::uvec&  ,
                       double );

double logPGamma(const arma::umat& ,
                 const arma::vec& ,
                 arma::mat* ,
                 const double ,
                 const double );
double logPGamma(const arma::umat& ,
                 const arma::vec& );
double logPGamma(const arma::umat& ,
                 arma::mat* mrfG,
                 const double ,
                 const double );

double logLikelihood(arma::uvec&  );

//arma::umat Gammas;
//arma::mat Pi;
//double a_pi, b_pi;
//double a_mrf, b_mrf;
//double logP_gamma, log_likelihood;

bool gamma_type_mrf;

// Bandit-sampling related quantities
int n_updates_bandit;
arma::vec banditZeta;
arma::mat banditAlpha;
arma::mat banditBeta;
arma::vec mismatch;
arma::vec normalised_mismatch;
arma::vec normalised_mismatch_backwards;

double banditLimit;
double banditIncrement;


#endif
