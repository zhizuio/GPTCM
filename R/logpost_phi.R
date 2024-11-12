#' @title Log posterior
#'
#' @description
#' This is a log posterior of the dispersion parameter phi of the Dirichlet distribution
#'
#' @name logpost_phi
#'
#' @param x
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
logpost_phi <- function(x, ...) {
  concentrations <- # dat$proportion # Here 'dat$proportion' should be replaced by updated alr
    proportion * x
  normalizingConst <- log(gamma(rowSums(concentrations))) -
    rowSums(log(gamma(concentrations)))

  geometricTerm <- rowSums((concentrations - 1) * # log(dat$proportion) )
                             log(proportion))

  Delta <- 20
  logprior <- log(truncnorm::dtruncnorm(x, a = 0, sd = Delta / 3))
  # alternative is half-cauchy 'extraDistr/src/half-cauchy-distribution.cpp'; or 'abs(rcauchy(n, location, scale))', or 'LaplacesDemon::rhalfcauchy'

  logpost <- sum(normalizingConst + geometricTerm) + logprior

  return(logpost)
}
