#' @title Log posterior
#'
#' @description
#' This is a log posterior of all coefficients of one component in the noncure fraction of GPTCM. (This seems worse for estimating cure fraction's effects, maybe due to non-concavity)
#'
#' @name logpost_beta_l
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
logpost_beta_l <- function(x) {
  # browser()
  betas.tmp <- betas.current
  # if(!is.na(x))
  #  betas.tmp[, l] <- x

  mu.tmp <- mu.current
  mu.tmp[, l] <- exp(dat$XX[, , l] %*% betas.tmp[, l])

  weibull.S.tmp <- weibull.S
  weibull.lambdas <- mu.tmp / gamma(1 + 1 / dat$kappas)
  weibull.S.tmp[, l] <- exp(-(dat$survObj$time / weibull.lambdas[, l])^kappas)

  # update beta's variance tauSq
  # This has been moved before updating beta, since ARMS for each tested beta should keep all others fixed
  # tauA <- 2; tauB <- 2
  # tauA.tmp <- tauA + 0.5 * sum(betas.tmp!=0)#as.numeric(x!=0)
  # tauB.tmp <- tauB + 0.5 * sum(betas.tmp^2)#x^2
  # tauSq <- 1/rgamma(length(x), shape = tauA.tmp, rate = tauB.tmp)

  logprior <- -sum(x^2 / tauSq / 2)

  logpost.first <- 0
  for (ll in 1:3) {
    logpost.first <- logpost.first +
      dat$proportion[, ll] *
        # kappas/(weibull.lambdas[,ll]) *
        # (dat$survObj$time/(weibull.lambdas[,ll]))^(kappas-1) *
        # exp(-(dat$survObj$time/(weibull.lambdas[,ll]))^kappas)
        weibull.lambdas[, ll]^(-kappas) * weibull.S.tmp[, ll]
  }

  logpost.first <- sum(log(logpost.first[dat$survObj$event == 1]))

  thetas.tmp <- exp(dat$x0 %*% xi)
  logpost.second <- sum(thetas.tmp * dat$proportion[, l] * weibull.S.tmp[, l])

  logpost <- logpost.first + logpost.second + logprior

  return(logpost)
}
