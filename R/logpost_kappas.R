#' @title Log posterior
#'
#' @description
#' This is a log posterior of Weibull's kappa
#'
#' @name logpost_kappas
#'
#' @importFrom stats dgamma
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
logpost_kappas <- function(x) {
  mu.tmp <- mu.current
  weibull.S.tmp <- weibull.S

  weibull.lambdas <- mu.tmp / gamma(1 + 1 / x)
  weibull.S.tmp[, l] <- exp(-(dat$survObj$time / weibull.lambdas[, l])^x)

  if (kappaSampler == "Gamma") {
    logprior <- dgamma(x, kappaA, kappaB, log = TRUE)
  } else {
    logprior <- log(1 / dgamma(x, kappaA, kappaB))
  }

  logpost.first <- logpost.second <- 0
  for (ll in 1:3) {
    logpost.first <- logpost.first + # dat$proportion[,ll] *
      proportion[, ll] * weibull.lambdas[, ll]^(-x) * weibull.S.tmp[, ll]
    logpost.second <- logpost.second + # dat$proportion[,ll] *
      proportion[, ll] * weibull.S.tmp[, ll]
  }

  logpost.first <- x * dat$survObj$time^(x - 1) * logpost.first
  logpost.first <- sum(log(logpost.first[dat$survObj$event == 1]))

  logpost.second <- sum(thetas * logpost.second)

  logpost <- logpost.first + logpost.second + logprior

  return(logpost)
}
