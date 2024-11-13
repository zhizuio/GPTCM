#' @title Log posterior
#'
#' @description
#' Log posterior of the proportion-related coefficients in the Dirichlet measure error model
#'
#' @name logpost_zeta_jl
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
logpost_zeta_jl <- function(x) {
  zetas.tmp <- zetas.current
  if (!is.na(x)) {
    zetas.tmp[j, l] <- x
  }

  # update proportions with a new proposal zeta
  proportion.tmp <- proportion
  for (ll in 1:(L - 1)) {
    proportion.tmp[, ll] <- exp(cbind(1, dat$XX[, , ll]) %*% zetas.tmp[, ll]) /
      (1 + rowSums(sapply(1:(L - 1), function(xx) {
        exp(cbind(1, dat$XX[, , xx]) %*% zetas.tmp[, xx])
      })))
  }
  proportion.tmp[, L] <- 1 - rowSums(proportion.tmp[, -L])

  # noncure density related censored part
  logpost.first <- 0
  for (ll in 1:L) {
    logpost.first <- logpost.first +
      proportion.tmp[, ll] *
        lambdas[, ll]^(-kappas) * weibull.S[, ll]
  }

  logpost.first <- sum(log(logpost.first[dat$survObj$event == 1]))

  # cure/noncure density combined the censored and noncensored events
  logpost.second <- sum(thetas * proportion.tmp[, l] * weibull.S[, l])

  # Dirichlet density
  concentrations <- proportion.tmp * phi # some issue here, since the sum of each row is phi
  normalizingConst <- log(gamma(rowSums(concentrations))) -
    rowSums(log(gamma(concentrations)))
  # normalizingConst <- -log(gamma(concentrations[,l])) # to be verified

  geometricTerm <- rowSums((concentrations - 1) * log(dat$proportion))
  # geometricTerm <- (concentrations[,l] - 1) * log(dat$proportion[,l]) # to be verified

  # prior
  logprior <- -x^2 / wSq / 2

  # sum all fractions
  logpost <- logpost.first + logpost.second +
    sum(normalizingConst + geometricTerm) + logprior

  return(logpost)
}
