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
  if (dirichlet) { # use logit/alr-link with the last category as reference
    for (ll in 1:(L - 1)) {
      proportion.tmp[, ll] <- exp(cbind(1, dat$XX[, , ll]) %*% zetas.tmp[, ll]) /
        (1 + rowSums(sapply(1:(L - 1), function(xx) {
          exp(cbind(1, dat$XX[, , xx]) %*% zetas.tmp[, xx])
        })))
    }
    proportion.tmp[, L] <- 1 - rowSums(proportion.tmp[, -L])
  } else { # use log-link for concentration parameters without reference category
    alphas <- matrix(NA, nrow = NROW(proportion), ncol = NCOL(proportion))
    for (ll in 1:L) {
      alphas[, ll] <- exp(cbind(1, dat$XX[, , ll]) %*% zetas.tmp[, ll])
    }
    proportion.tmp <- apply(alphas, 2, function(xx) {
      xx / rowSums(alphas)
    })
  }

  # noncure density related censored part
  logpost.first <- logpost.second <- 0
  for (ll in 1:L) {
    tmp <- proportion.tmp[, ll] * weibull.S[, ll]
    logpost.first <- logpost.first + lambdas[, ll]^(-kappas) * tmp
    logpost.second <- logpost.second + tmp
  }

  logpost.first <- sum(log(logpost.first[dat$survObj$event == 1]))

  # cure/noncure density combined the censored and noncensored events
  # logpost.second <- sum(thetas * proportion.tmp[, l] * weibull.S[, l])
  logpost.second <- sum(thetas * logpost.second)

  # Dirichlet density
  concentrations <- proportion.tmp * phi # some issue here, since the sum of each row is phi
  normalizingConst <- # log(gamma(rowSums(concentrations))) -
    log(gamma(phi)) -
    rowSums(log(gamma(concentrations)))
  # normalizingConst <- -log(gamma(concentrations[,l])) # to be verified

  geometricTerm <- rowSums((concentrations - 1) * log(dat$proportion))
  # geometricTerm <- (concentrations[,l] - 1) * log(dat$proportion[,l]) # to be verified

  # prior
  w.tmp <- ifelse(j == 1, w0Sq, wSq)
  logprior <- -x^2 / w.tmp / 2

  # sum all fractions
  logpost <- logpost.first + logpost.second +
    sum(normalizingConst + geometricTerm) + logprior

  return(logpost)
}
