#' @title Log posterior
#'
#' @description
#' This is a log posterior of each coefficients in the noncure fraction of GPTCM
#'
#' @name logpost_beta_jl
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
logpost_beta_jl <- function(x, ...) {
  # browser()
  betas.tmp <- betas.current
  if (!is.na(x)) {
    betas.tmp[j, l] <- x
  }

  mu.tmp <- mu.current
  mu.tmp[, l] <- exp(dat$XX[, , l] %*% betas.tmp[, l])

  weibull.S.tmp <- weibull.S
  weibull.lambdas <- mu.tmp / gamma(1 + 1 / kappas)
  weibull.S.tmp[, l] <- exp(-(dat$survObj$time / weibull.lambdas[, l])^kappas)

  # update beta's variance tauSq
  # tauSq <- sampleTau(tauA, tauB, betas.tmp) # This has been moved before updating beta.

  logprior <- -x^2 / tauSq / 2

  logpost.first <- 0
  for (ll in 1:3) {
    logpost.first <- logpost.first +
      # dat$proportion[,ll] *
      proportion[, ll] *
        # kappas/(weibull.lambdas[,ll]) *
        # (dat$survObj$time/(weibull.lambdas[,ll]))^(kappas-1) *
        # exp(-(dat$survObj$time/(weibull.lambdas[,ll]))^kappas)
        weibull.lambdas[, ll]^(-kappas) * weibull.S.tmp[, ll]
  }

  logpost.first <- sum(log(logpost.first[dat$survObj$event == 1]))

  thetas.tmp <- exp(dat$x0 %*% xi)
  logpost.second <- sum(thetas.tmp * # dat$proportion[, l]
                          proportion[, l] * weibull.S.tmp[, l])

  logpost <- logpost.first + logpost.second + logprior

  return(logpost)
}
