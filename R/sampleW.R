#' @title Gibbs sampling
#'
#' @description
#' Gibbs sampling for the variance wSq of nonzero zetas
#'
#' @name sampleW
#'
#' @importFrom stats rgamma
#'
#' @param n.sample TBA
#' @param externalwA TBA
#' @param externalwB TBA
#' @param externalZetas TBA
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
sampleW <- function(n.sample, externalwA, externalwB, externalZetas) {
  wA.post <- externalwA + 0.5 * sum(externalZetas != 0)
  wB.post <- externalwB + 0.5 * sum(externalZetas^2)
  wSq <- 1 / rgamma(n.sample, shape = wA.post, rate = wB.post)

  return(wSq)
}
