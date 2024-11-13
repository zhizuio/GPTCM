#' @title Log posterior
#'
#' @description
#' This is a log posterior of each coefficients (including intercept) in the cure fraction of GPTCM. (This is slightly less efficient than 'logpost_xi2')
#'
#' @name logpost_xi
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
logpost_xi <- function(x, ...) {
  cat("Currently function 'logpost_xi()' is not implemented well!")
  xi.tmp <- xi
  xi.tmp[jj] <- x
  weibull.S.tmp <- weibull.S

  if (jj == 1) {
    vSq <- 10
  } else {
    vA.tmp <- vA + 0.5 * sum(xi.tmp[-1] != 0)
    vB.tmp <- vB + 0.5 * sum(xi.tmp[-1]^2)
    vSq <- 1 / rgamma(1, shape = vA.tmp, rate = vB.tmp)
    # vSq <- 10
  }
  if (is.na(vSq)) browser()
  logprior <- -x^2 / vSq / 2
  # logprior <- 0 # noninformative prior

  thetas.tmp <- exp(dat$x0 %*% xi.tmp)
  logpost.first <- sum(log(thetas.tmp[dat$survObj$event == 1]))

  logpost.second <- 0
  for (ll in 1:3) {
    logpost.second <- logpost.second +
      dat$proportion[, ll] * weibull.S.tmp[, ll]
  }
  logpost.second <- -sum(thetas.tmp * (1 - logpost.second))

  logpost <- logpost.first + logpost.second + logprior

  return(logpost)
}
