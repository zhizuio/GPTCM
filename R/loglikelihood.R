#' @title Log-likelihood
#'
#' @description
#' Log-likelihood GTPCM
#'
#' @name loglikelihood
#'
#' @param par TBA
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
loglikelihood <- function(par) {
  # Need data variables: survObj, X0, X, proportion
  L <- dim(X)[3]
  p <- dim(X)[2]
  n <- dim(X)[1]
  d <- dim(X0)[2]

  if (any(is.na(par))) {
    loglik <- -Inf
  } else {
    kappas <- par[1]
    phi <- par[2]
    xi <- par[2 + 1:d]
    betas <- matrix(par[2 + d + 1:(p * L)], ncol = L)
    zetas <- matrix(par[-c(1:(2 + d + p * L))], nrow = p + 1)

    thetas <- as.vector(exp(X0 %*% xi))
    # lambdas <- mu <- matrix(0, nrow = n, ncol = L)
    survival.pop <- rep(0, n)
    f <- 0
    proportion.tmp <- proportion
    for (l in 1:L) {
      mu.l <- exp(X[, , l] %*% betas[, l])
      lambdas.l <- mu.l / gamma(1 + 1 / kappas)

      f <- f + proportion[, l] * (lambdas.l)^(-kappas) *
        exp(-(survObj$time / lambdas.l)^kappas)

      if (l < L) {
        proportion.tmp[, l] <- exp(cbind(1, X[, , l]) %*% zetas[, l]) /
          (1 + rowSums(sapply(1:(L - 1), function(xx) {
            exp(cbind(1, X[, , xx]) %*% zetas[, xx])
          })))
      }

      survival.l <- exp(-(survObj$time / lambdas.l)^kappas)
      survival.pop <- survival.pop + proportion.tmp[, l] * survival.l
    }
    if (any(is.na(proportion.tmp))) browser()
    # if (any(proportion.tmp[, -L] == 1)) {
    #   ones.idx <- which(proportion.tmp[, -L] == 1)
    #   proportion.tmp[ones.idx] <- 1 - 1e-10
    # }
    proportion.tmp[, L] <- 1 - rowSums(proportion.tmp[, -L])
    if (any(proportion.tmp[, L] == 0)) {
      zero.idx <- (proportion.tmp[, L] == 0)
      proportion.tmp[zero.idx, L] <- 1e-10
    }

    # Survival part
    log.survival.pop <- (-thetas * (1 - survival.pop))
    log.f <- log(kappas) + (kappas - 1) * log(survObj$time) + log(f)
    log.f.pop <- log(thetas) + log.f + log.survival.pop

    # Dirichlet part
    concentrations <- proportion.tmp * phi
    if (any(is.na(gamma(concentrations)))) browser()
    normalizingConst <- log(gamma(phi)) - rowSums(log(gamma(concentrations)))
    geometricTerm <- rowSums((concentrations - 1) * log(proportion))
    log.dirichlet <- normalizingConst + geometricTerm

    # log-likelihood
    loglik <- sum(log.f.pop[survObj$event == 1]) +
      sum(log.survival.pop[survObj$event == 0]) +
      sum(log.dirichlet)
  }

  return(-loglik) # - log-likelihood
}
