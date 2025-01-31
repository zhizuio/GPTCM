#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include "sampleGamma.h"

// [[Rcpp::export]]
Rcpp::List sampleGamma(arma::umat Gammas,
                      bool gamma_sampler_bandit,
                      double a_pi, 
                      double b_pi)
{
    //double logP_gamma = logPGamma( Gammas );
    //double log_likelihood = logLikelihood( Gammas );

    arma::umat proposedGamma;
    arma::uvec updateIdx;

    int p = Gammas.n_rows;
    int L = Gammas.n_cols;

    // decide on one component
    int componentUpdateIdx = std::floor( R::runif( 0, L ) );

    // update VS probability corresponding to the decided component
    double pi = R::rbeta( a_pi + arma::sum(Gammas.col(componentUpdateIdx)), 
                    b_pi + p - arma::sum(Gammas.col(componentUpdateIdx)) );
    
    // Update the proposed Gamma
    switch ( gamma_sampler_bandit )
    {
        case 1 :
            proposedGamma = gammaBanditProposal( Gammas, updateIdx, componentUpdateIdx );
            // logProposalRatio = ??? // to be updated
            break;
            
        case 0 :
            proposedGamma = gammaMC3Proposal( Gammas, updateIdx, componentUpdateIdx );
            break;
            
        default:
            break;
    }
    
    // note only one outcome is updated
    // update log probabilities

    double logProposalGammaRatio = 0.;
    double logProposalBetaRatio = 0.;
    double logLikelihoodRatio = 0.;

    // compute logProposalGammaRatio, i.e. proposedGammaPrior - logP_gamma
    for(unsigned int i=0; i<updateIdx.n_elem; ++i)
    {
        logProposalGammaRatio += logPDFBernoulli( proposedGamma(i,componentUpdateIdx), pi ) - logPDFBernoulli( Gammas(i,componentUpdateIdx), pi );
    }

    // compute logProposalBetaRatio given proposedGamma, i.e. proposedBetaPrior - logP_beta
    // here resample beta_componentUpdateIdx and compute 
    // sampleBetaKGivenGammas = to be added
    // ...
    logProposalBetaRatio = 0.;

    // compute logLikelihoodRatio, i.e. proposedLikelihood - log_likelihood
    logLikelihoodRatio = 0.;
    
    double logAccProb = logProposalGammaRatio + logProposalBetaRatio + logLikelihoodRatio;
    
    if( R::runif(0,1) < logAccProb )
    {
        Gammas = proposedGamma;
        
        //logP_gamma = proposedGammaPrior;
        //log_likelihood = proposedLikelihood;
        
        // ++gamma_acc_count;
        //gamma_acc_count += 1.; // / updatedOutcomesIdx.n_elem;
    }
    
    // after A/R, update bandit Related variables
    if( gamma_sampler_bandit )
    {
        for(arma::uvec::iterator iter = updateIdx.begin(); iter != updateIdx.end(); ++iter)
        {
            // FINITE UPDATE
            if( banditAlpha(*iter,componentUpdateIdx) + banditBeta(*iter,componentUpdateIdx) <= banditLimit )
            {
                banditAlpha(*iter,componentUpdateIdx) += banditIncrement * Gammas(*iter,componentUpdateIdx);
                banditBeta(*iter,componentUpdateIdx) += banditIncrement * (1-Gammas(*iter,componentUpdateIdx));
            }
        }
    }
    
    return Rcpp::List::create(
                              Rcpp::Named("Gammas") = Gammas,
                              Rcpp::Named("proposedGamma") = proposedGamma,
                              Rcpp::Named("logAccProb") = logAccProb,
                              Rcpp::Named("logProposalGammaRatio") = logProposalGammaRatio,
                              Rcpp::Named("logProposalBetaRatio") = logProposalBetaRatio,
                              Rcpp::Named("logLikelihoodRatio") = logLikelihoodRatio
                                  ); 
    //return Gammas;
}
/*
arma::vec samplePi(arma::umat& Gammas)
{
    // re-initialize the variable selection probabilities
    arma::vec Pi = arma::zeros<arma::vec>( Gammas.n_cols ); 

    for(unsigned int j=0; j<Gammas.n_cols; ++j)
    {
        Pi(j) = R::rbeta( a_pi + arma::sum(Gammas.col(j)), b_pi + p arma::sum(Gammas.col(j)));
    }

    return Pi;
}
*/
arma::umat gammaMC3Proposal(arma::umat& Gammas,
                            arma::uvec& updateIdx,
                            const int componentUpdateIdx_ )
{
    arma::umat mutantGamma = Gammas;
    int p = Gammas.n_rows;
    int n_updates_MC3 = std::ceil( p / 5 ); //arbitrary number, should I use something different?
    updateIdx = arma::uvec( n_updates_MC3 );
    
    for( unsigned int i=0; i<n_updates_MC3; ++i)
    {
        updateIdx(i) = std::floor( R::runif( 0, p ) );    // note that I might be updating multiple times the same coeff
    }
    
    for( auto i : updateIdx)
    {
        mutantGamma(i,componentUpdateIdx_) = ( R::runif(0,1) < 0.5 )? Gammas(i,componentUpdateIdx_) : 1-Gammas(i,componentUpdateIdx_); // could simply be ( 0.5 ? 1 : 0) ;
    }
    
    return mutantGamma ; // pass this to the outside, it's the (symmetric) logProposalRatio
}


// sampler for proposed updates on Gammas
arma::umat gammaBanditProposal(arma::umat& Gammas,
                           arma::uvec& updateIdx,
                           int componentUpdateIdx_ )
{
    
    int p = Gammas.n_elem;
    arma::umat mutantGamma = Gammas;
    double logProposalRatio{0.};
    
    n_updates_bandit = 4; // this needs to be low as its O(n_updates!)

    int nVSPredictors = n_updates_bandit; // to be changed
    
    // Sample Zs (only for relevant component)
    for(unsigned int i=0; i<nVSPredictors; ++i)
    {
        banditZeta(i) = R::rbeta(banditAlpha(i,componentUpdateIdx_),banditAlpha(i,componentUpdateIdx_));
    }
    
    // Create mismatch (only for relevant outcome)
    for(unsigned int i=0; i<nVSPredictors; ++i)
    {
        mismatch(i) = (mutantGamma(i,componentUpdateIdx_)==0)?(banditZeta(i)):(1.-banditZeta(i));   //mismatch
    }
    
    // Normalise
    // mismatch = arma::log(mismatch); //logscale ??? TODO
    // normalised_mismatch = mismatch - Utils::logspace_add(mismatch);
    
    normalised_mismatch = mismatch / arma::as_scalar(arma::sum(mismatch));
    
    if( R::runif(0,1) < 0.5 )   // one deterministic update
    {
        // Decide which to update
        updateIdx = arma::zeros<arma::uvec>(1);
        updateIdx(0) = randWeightedIndexSampleWithoutReplacement(nVSPredictors,normalised_mismatch); // sample the one
        
        // Update
        mutantGamma(updateIdx(0),componentUpdateIdx_) = 1 - Gammas(updateIdx(0),componentUpdateIdx_); // deterministic, just switch

        // Compute logProposalRatio probabilities
        normalised_mismatch_backwards = mismatch;
        normalised_mismatch_backwards(updateIdx(0)) = 1. - normalised_mismatch_backwards(updateIdx(0)) ;
        
        // normalised_mismatch_backwards = normalised_mismatch_backwards - Utils::logspace_add(normalised_mismatch_backwards);
        normalised_mismatch_backwards = normalised_mismatch_backwards / arma::as_scalar(arma::sum(normalised_mismatch_backwards));
        
        logProposalRatio = ( std::log( normalised_mismatch_backwards(updateIdx(0)) ) ) -
        ( std::log( normalised_mismatch(updateIdx(0)) ) );
        
    }else{
        
        /*
         n_updates_bandit random (bern) updates
         Note that we make use of column indexing here for armadillo matrices
         */

        logProposalRatio = 0.;
        // Decide which to update
        updateIdx = arma::zeros<arma::uvec>(n_updates_bandit);
        updateIdx = randWeightedIndexSampleWithoutReplacement(nVSPredictors,normalised_mismatch,n_updates_bandit); // sample n_updates_bandit indexes
        
        normalised_mismatch_backwards = mismatch; // copy for backward proposal
        
        // Update
        for(unsigned int i=0; i<n_updates_bandit; ++i)
        {
            mutantGamma(updateIdx(i),componentUpdateIdx_) = R::rbinom( 1, banditZeta(updateIdx(i))); // random update
            
            normalised_mismatch_backwards(updateIdx(i)) = 1.- normalised_mismatch_backwards(updateIdx(i));
            
            logProposalRatio += logPDFBernoulli(Gammas(updateIdx(i),componentUpdateIdx_),banditZeta(updateIdx(i))) -
            logPDFBernoulli(mutantGamma(updateIdx(i),componentUpdateIdx_),banditZeta(updateIdx(i)));
        }

        // note that above I might be resampling a value equal to the current one, thus not updating da facto ... TODO
        
        // Compute logProposalRatio probabilities
        // normalised_mismatch_backwards = normalised_mismatch_backwards - Utils::logspace_add(normalised_mismatch_backwards);
        normalised_mismatch_backwards = normalised_mismatch_backwards / arma::as_scalar(arma::sum(normalised_mismatch_backwards));
        
        logProposalRatio += logPDFWeightedIndexSampleWithoutReplacement(normalised_mismatch_backwards,updateIdx) -
        logPDFWeightedIndexSampleWithoutReplacement(normalised_mismatch,updateIdx);
    }
    
    //return logProposalRatio; // pass this to the outside
    return mutantGamma;
            
}
                                  
// subfunctions used for bandit proposal

arma::uvec randWeightedIndexSampleWithoutReplacement
(
    unsigned int populationSize,    // size of set sampling from
    const arma::vec& weights,       // (log) probability for each element
    unsigned int sampleSize         // size of each sample
) // sample is a zero-offset indices to selected items, output is the subsampled population.
{
    // note I can do everything in the log scale as the ordering won't change!
    arma::vec tmp = Rcpp::rexp( populationSize, 1. );
    arma::vec score = tmp - weights;
    arma::uvec result = arma::sort_index( score,"ascend" );

    return result.subvec(0,sampleSize-1);
}

// Overload with equal weights
arma::uvec randWeightedIndexSampleWithoutReplacement
(
    unsigned int populationSize,    // size of set sampling from
    unsigned int sampleSize         // size of each sample
) // sample is a zero-offset indices to selected items, output is the subsampled population.
{
    // note I can do everything in the log scale as the ordering won't change!
    arma::vec score = Rcpp::rexp( populationSize, 1. );
    arma::uvec result = arma::sort_index( score,"ascend" );
    
    return result.subvec(0,sampleSize-1);
}

// overload with sampleSize equal to one
arma::uword randWeightedIndexSampleWithoutReplacement
(
    unsigned int populationSize,    // size of set sampling from
    const arma::vec& weights     // probability for each element
) // sample is a zero-offset indices to selected items, output is the subsampled population.
{
    // note I can do everything in the log scale as the ordering won't change!

    double u = R::runif(0,1);
    double tmp = weights(0);
    unsigned int t = 0;

    while(u > tmp)
    {
        // tmp = Utils::logspace_add(tmp,logWeights(++t));
        tmp += weights(++t);
    }

    return t;
}

// logPDF rand Weighted Indexes (need to implement the one for the original starting vector?)
double logPDFWeightedIndexSampleWithoutReplacement(const arma::vec& weights,
                                                   const arma::uvec& indexes)
{
    // arma::vec logP_permutation = arma::zeros<arma::vec>((int)std::tgamma(indexes.n_elem+1));  //too big of a vector
    double logP_permutation = 0.; double tmp;

    std::vector<unsigned int> v = arma::conv_to<std::vector<unsigned int>>::from(arma::sort(indexes));
    // vector should be sorted at the beginning.

    arma::uvec current_permutation;
    arma::vec current_weights;

    do {
        current_permutation = arma::conv_to<arma::uvec>::from(v);
        current_weights = weights;
        tmp = 0.;

        while( current_permutation.n_elem > 0 )
        {
           tmp += log(current_weights(current_permutation(0)));
           current_permutation.shed_row(0);
           current_weights = current_weights/arma::sum(current_weights(current_permutation));   // this will gets array weights that do not sum to 1 in total, but will only use relevant elements
        }

        logP_permutation = logspace_add( logP_permutation,tmp );

    } while (std::next_permutation(v.begin(), v.end()));

    return logP_permutation;
}

double logspace_add(double a,
                    double b)
{

    if(a <= std::numeric_limits<float>::lowest())
          return b;
    if(b <= std::numeric_limits<float>::lowest())
          return a;
    return std::max(a, b) + std::log( (double)(1. + std::exp( (double)-std::abs((double)(a - b)) )));
}

double logPDFBernoulli(int x, double pi)
{
    if( x > 1. ||  x < 0. )
        return -std::numeric_limits<double>::infinity();
    else
        return x*std::log(pi) + (1-x)*std::log(1.-pi);
}

double logPDFBernoulli(const arma::uvec& x, double pi)
{
	double n = (double)x.n_elem;
	double k = arma::sum(x);

	return k*std::log(pi) + (n-k)*std::log(1.-pi);
}

// these below are general interfaces
double logPGamma(const arma::umat& Gammas , 
                    const arma::vec& pi ,
                    arma::mat* mrfG,
                    const double a_mrf ,
                    const double b_mrf )
{
    double logP_gamma;

    switch ( gamma_type_mrf )
    {
        case 1 :
            logP_gamma = logPGamma( Gammas , mrfG, a_mrf , b_mrf );
            break;

        case 0 :
            logP_gamma = logPGamma( Gammas , pi );
            break;

        default:
            break;
    }
    return logP_gamma;
}

// this is the simple hierarchical prior
double logPGamma( const arma::umat& Gammas , 
                    const arma::vec& pi )
{
    // if( gamma_type_mrf )
    //     break;
    
    double logP = 0.;
    for(unsigned int j=0; j<Gammas.n_cols; ++j)
    {
        logP += logPDFBernoulli( Gammas.col(j) , pi(j) );
        // logP += logPDFBinomial( arma::sum( Gammas.row(j) ) , nOutcomes , pi(j) ); // do we care about the binomial coeff? I don't think so..
    }
    return logP;
}

// this is the MRF prior
double logPGamma( const arma::umat& Gammas, 
                    arma::mat* mrfG,
                    const double a, 
                    const double b )
{
    // if( !gamma_type_mrf )
    //     break;
    
    //arma::mat externalMRFG = mrfG->cols( arma::linspace<arma::uvec>(0,2,3) );
    
    double logP = 0.;
    // calculate the linear and quadratic parts in MRF by using all edges of G
    arma::vec gammaVec = arma::conv_to< arma::vec >::from(arma::vectorise(Gammas));
    double quad_mrf = 0.;
    double linear_mrf = 0.;
    //int count_linear_mrf = 0; // If the MRF graph matrix has diagonals 0, count_linear_mrf is always 0.
    for( unsigned i=0; i < mrfG->n_rows; ++i )
    {
        if( (*mrfG)(i,0) != (*mrfG)(i,1) ){
            quad_mrf += 2.0 * gammaVec( (*mrfG)(i,0) ) * gammaVec( (*mrfG)(i,1) ) * (*mrfG)(i,2);
        }else{
                if( gammaVec( (*mrfG)(i,0) ) == 1 ){
                    linear_mrf += (*mrfG)(i,2); // should this be 'linear_mrf += e * (externalMRFG)(i,2)'?
                    //count_linear_mrf ++;
                }
        }
    }
    //logP = arma::as_scalar( linear_mrf + d * (arma::accu( Gammas ) - count_linear_mrf) + e * 2.0 * quad_mrf );
    // Should logP be the following?
    logP = arma::as_scalar( a * arma::accu( Gammas ) + b * (linear_mrf + quad_mrf) );
    
    return logP;
}

double logLikelihood(arma::uvec& Gammas )
{
    double logP = 0.;
    // TBC...
}

