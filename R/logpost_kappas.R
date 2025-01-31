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
  # mu.tmp <- mu.current
  # weibull.S.tmp <- weibull.S
  #
  # weibull.lambdas <- mu.tmp / gamma(1 + 1 / x)
  # weibull.S.tmp[, l] <- exp(-(dat$survObj$time / weibull.lambdas[, l])^x)
  
  if (is.na(x)) {
    logpost <- -Inf
  } else {
    # print(head(proportion))
    # print(x)
    
    if (kappaPrior == "Gamma") {
      logprior <- dgamma(x, kappaA, kappaB, log = TRUE)
    } else {
      logprior <- log(1 / dgamma(x, kappaA, kappaB))
    }
    
    logpost.first <- logpost.second <- 0
    for (ll in 1:3) {
      weibull.lambdas.tmp <- mu.current[, ll] / gamma(1 + 1 / x)
      weibull.S.tmp <- exp(-(dat$survObj$time / weibull.lambdas.tmp)^x)
      
      logpost.first <- logpost.first + # dat$proportion[,ll] *
        proportion[, ll] * weibull.lambdas.tmp^(-x) * weibull.S.tmp
      logpost.second <- logpost.second + # dat$proportion[,ll] *
        proportion[, ll] * weibull.S.tmp
    }
    
    logpost.first <- x * dat$survObj$time^(x - 1) * logpost.first
    logpost.first <- sum(log(logpost.first[dat$survObj$event == 1]))
    
    logpost.second <- sum(thetas * logpost.second)
    
    logpost <- logpost.first + logpost.second + logprior
  }

  return(logpost)
}
