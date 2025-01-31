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
    xi <- par[1 + 1:d]
    betas <- matrix(par[1 + d + 1:(p * L)], ncol = L)
    
    # set bounds to avoid numeric issues
    eps <- 1e-10
    eps1 <- 170
    eps2 <- 700
    
    proportion.tmp <- proportion.obs
    if (proportion.model) {
      phi <- par[1 + d + p * L + 1]
      zetas <- matrix(par[-c(1:(1 + d + p * L + 1))], nrow = p + 1)
      
      # log(concentrations) model for proportions
      logAlphas <- sapply(1:L, function(ll)
        {cbind(1, X[, , ll]) %*% zetas[, ll]}
      )
      # alphas[alphas > eps1] <- eps1
      # proportion.tmp <- alphas / rowSums(alphas)
      #proportion.tmp <- exp(logAlphas - rowSums(logAlphas))
      proportion.tmp <- exp(logAlphas - log(rowSums(exp(logAlphas))))
      
      # if (any(proportion.tmp < eps) || any(proportion.tmp > 1 - eps)) {
      #   proportion.tmp <- (proportion.tmp * (n - 1) + 1 / L) / n
      # }
    } #else {
      #proportion.tmp <- proportion.obs
    #}
    
    thetas <- as.vector(exp(X0 %*% xi))
    # lambdas <- mu <- matrix(0, nrow = n, ncol = L)
    survival.pop <- rep(0, n)
    f <- 0
    
    for (l in 1:L) {
      xbeta <- X[, , l] %*% betas[, l]
      xbeta[xbeta > eps2] <- eps2
      mu.l <- exp(xbeta)
      lambdas.l <- mu.l / gamma(1 + 1 / kappas)
      
      density.weibull <- exp(-kappas * log(lambdas.l) - (survObj$time / lambdas.l)^kappas)
      
      survival.l <- exp(-(survObj$time / lambdas.l)^kappas)
      # survival.l[survival.l < eps]   <- eps
      # survival.l[survival.l > 1-eps] <- 1-eps
      survival.pop <- survival.pop + proportion.tmp[, l] * survival.l
      
      f <- f + kappas * survObj$time^(kappas - 1) * proportion.tmp[, l] * density.weibull
      
      # f <- f + proportion.tmp[, l] * (lambdas.l)^(-kappas) *
      #   exp(-(survObj$time / lambdas.l)^kappas)
      # if (proportion.model) {
      #   if (l < L) {
      #     proportion.tmp[, l] <- exp(cbind(1, X[, , l]) %*% zetas[, l]) /
      #       (1 + rowSums(sapply(1:(L - 1), function(xx) {
      #         exp(cbind(1, X[, , xx]) %*% zetas[, xx])
      #       })))
      #   }
      # } else {
      #   proportion.tmp <- proportion
      # }

      # survival.l <- exp(-(survObj$time / lambdas.l)^kappas)
      # survival.pop <- survival.pop + proportion.tmp[, l] * survival.l
    }
    
    #if (any(is.na(proportion.tmp))) browser()
    # if (any(proportion.tmp[, -L] == 1)) {
    #   ones.idx <- which(proportion.tmp[, -L] == 1)
    #   proportion.tmp[ones.idx] <- 1 - 1e-10
    # }
    # proportion.tmp[, L] <- 1 - rowSums(proportion.tmp[, -L])
    # if (any(proportion.tmp[, L] == 0)) {
    #   zero.idx <- (proportion.tmp[, L] == 0)
    #   proportion.tmp[zero.idx, L] <- 1e-10
    # }

    # Survival part
    log.survival.pop <- (-thetas * (1 - survival.pop))
    #log.f <- log(kappas) + (kappas - 1) * log(survObj$time) + log(f)
    # f[f < eps] <- eps
    log.f <- log(f)
    log.f.pop <- log(thetas) + log.f + log.survival.pop
    #if(any(is.na(log.f))) browser()
    
    # Dirichlet part
    log.dirichlet <- 0
    if (proportion.model) {
      # concentrations <- proportion.tmp * phi
      # if (any(is.na(gamma(concentrations)))) browser()
      # normalizingConst <- log(gamma(phi)) - rowSums(log(gamma(concentrations)))
      
      ## some numeric issue below
      # alphas <- exp(logAlphas)
      # normalizingConst <- log(gamma(rowSums(alphas))) -
      #   rowSums(log(gamma(alphas)))
      # geometricTerm <- rowSums((alphas - 1) * log(proportion.obs))
      # log.dirichlet <- normalizingConst + geometricTerm
      
      #lmvbeta <- rowSums(lgamma(exp(logAlphas))) - lgamma(rowSums(exp(logAlphas)))
      log.dirichlet <- #(-lmvbeta) +
        lgamma(rowSums(exp(logAlphas))) - rowSums(lgamma(exp(logAlphas))) +
        rowSums((exp(logAlphas) -1) * log(proportion.obs))
      log.dirichlet <- sum(log.dirichlet)
    }
    #if(any(is.na(log.dirichlet))) browser()

    # log-likelihood
    loglik <- sum(log.f.pop[survObj$event == 1]) +
      sum(log.survival.pop[survObj$event == 0]) +
      log.dirichlet
  }
  
if(is.na(loglik))browser()
  return(-loglik) # - log-likelihood
}
