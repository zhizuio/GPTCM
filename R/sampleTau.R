#' @title Gibbs sampling
#'
#' @description
#' Gibbs sampling for the variance tauSq of nonzero betas
#'
#' @name sampleTau
#'
#' @importFrom stats rgamma
#' 
#' @param n.sample
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
sampleTau <- function(n.sample, externalTauA, externalTauB, externalBetas) {
  tauA.post <- externalTauA + 0.5 * sum(externalBetas != 0)
  tauB.post <- externalTauB + 0.5 * sum(externalBetas^2)
  tauSq <- 1 / rgamma(n.sample, shape = tauA.post, rate = tauB.post)

  return(tauSq)
}
