#' @title Gibbs sampling
#'
#' @description
#' Gibbs sampling for the variance v of xi
#'
#' @name sampleV
#'
#' @importFrom stats rgamma
#' 
#' @param n.sample TBA
#' @param externalvA TBA
#' @param externalvB TBA
#' @param externalXi TBA
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
sampleV <- function(n.sample, externalvA, externalvB, externalXi) {
  vA.post <- externalvA + 0.5 * sum(externalXi[-1] != 0)
  vB.post <- externalvB + 0.5 * sum(externalXi[-1]^2)
  v <- 1 / rgamma(n.sample, shape = vA.post, rate = vB.post)

  return(v)
}
