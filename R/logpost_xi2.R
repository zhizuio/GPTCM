#' @title Log posterior
#'
#' @description
#' Log posterior of all coefficients (including intercept) in the cure fraction of GPTCM
#'
#' @name logpost_xi2
#'
#' @param x proposal value
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
logpost_xi2 <- function(x) {
  weibull.S.tmp <- weibull.S

  # vA.tmp <- vA + 0.5 * sum(x[-1]!=0)
  # vB.tmp <- vB + 0.5 * sum(x[-1]^2)
  # vSq <- c(10, 1/rgamma(length(x)-1, shape = vA.tmp, rate = vB.tmp))

  logprior <- -sum(x^2 / vSq / 2)

  # thetas.tmp <- exp(dat$x0 %*% x)
  # logpost.first <- sum(log(thetas.tmp[dat$survObj$event == 1]))
  thetas.tmp <- exp(datX0 %*% x)
  logpost.first <- sum(log(thetas.tmp[datEvent == 1]))

  logpost.second <- 0
  for (ll in 1:3) {
    logpost.second <- logpost.second + # dat$proportion[, ll] *
      proportion[, ll] * weibull.S.tmp[, ll]
  }
  logpost.second <- -sum(thetas.tmp * (1 - logpost.second))

  logpost <- logpost.first + logpost.second + logprior

  return(logpost)
}
