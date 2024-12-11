#' @title Fit Bayesian GPTCM Model
#'
#' @description
#' This is the main function to fit the Bayesian GPTCM model with multiscale data
#' for sparse identification of high-dimensional covariates.
#'
#' @name GPTCM
#'
#' @importFrom stats median nlminb
#' @importFrom HI arms
#' @import remotes
#'
#' @param dat a list containing observed data from \code{n} subjects with
#' components \code{t}, \code{di}, \code{X}. For graphical learning of the
#' Markov random field prior, \code{survObj} should be a list of the list with
#' survival and covariates data. For subgroup models with or without graphical
#' learning, \code{survObj} should be a list of multiple lists with each
#' component list representing each subgroup's survival and covariates data
#' @param n TBA
#' @param p TBA
#' @param L TBA
#' @param proportion.model logical value; should the proportions be modeled or
#' not. If (\code{proportion.model = FALSE}), the argument \code{dirichlet} will
#' be invalid
#' @param dirichlet logical value; should the proportions be modeled via the
#' common (\code{dirichlet = FALSE}) or alternative (\code{dirichlet = TRUE})
#' parametrization of the Dirichlet regression model
#' @param method either `MLE` or `Bayes`
#' @param hyperpar TBA
#' @param kappaPrior TBA
#' @param kappaSampler TBA
#' @param w0Sampler TBA
#' @param initial TBA
#' @param nIter TBA
#' @param burnin TBA
#' @param thin TBA
#' @param tick TBA
#'
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
GPTCM <- function(dat, n, p, L,
                  proportion.model = TRUE,
                  dirichlet = TRUE,
                  method = "Bayes",
                  hyperpar = NULL,
                  kappaPrior = "IGamma",
                  kappaSampler = "arms",
                  w0Sampler = "IGamma",
                  initial = NULL,
                  nIter = 1,
                  burnin = 0,
                  thin = 1,
                  tick = 100) {
  # Validation
  stopifnot(burnin < nIter)

  # check the formula
  cl <- match.call()

  # set hyperparamters of all piors
  if (is.null(hyperpar)) {
    hyperpar$tauA <- 20
    hyperpar$tauB <- 50
    # hyperpar$tauA <- 1.5; hyperpar$tauB <- 50; tauSq <- 1
    # hyperpar$wA <- 20; hyperpar$wB <- 50; wSq <- 1
    hyperpar$wA <- 5
    hyperpar$wB <- 20
    hyperpar$vA <- 10
    hyperpar$vB <- 20
    if (kappaPrior == "Gamma") {
      # hyperpar$kappaA <- 3 # 1
      hyperpar$kappaA <- 1
      hyperpar$kappaB <- 1 # This is for Gamma prior
    } else {
      if (kappaPrior == "IGamma") {
        hyperpar$kappaA <- 5
        # hyperpar$kappaB <- 20 # This is for Inverse-Gamma prior
        hyperpar$kappaB <- 5
      } else {
        stop("Argument 'kappaPrior' is either 'Gamma' or 'IGamma'!")
      }
    }
    if (w0Sampler == "IGamma") {
      hyperpar$w0A <- hyperpar$wA
      hyperpar$w0B <- hyperpar$wB
    }
  }

  # initialization of parameters
  if (is.null(initial)) {
    tauSq <- 1
    wSq <- 1
    if (w0Sampler == "IGamma") {
      w0Sq <- 1
    } else {
      if (is.numeric(w0Sampler)) {
        w0Sq <- w0Sampler
      }
    }
    vSq <- c(10, 1, 1)
    kappas <- 0.9

    betas.current <- matrix(0, nrow = dim(dat$XX)[2], ncol = NCOL(dat$proportion))
    mu.current <- matrix(0, nrow = dim(dat$XX)[1], ncol = NCOL(dat$proportion))
    for (l in 1:L) {
      mu.current[, l] <- exp(dat$beta0[l] + dat$XX[, , l] %*% betas.current[, l])
    }
    lambdas <- mu.current / gamma(1 + 1 / kappas) # lambdas is a parameter in WEI3 distr.

    ## proportion Dirichlet part
    phi <- 1
    zetas.current <- matrix(0, nrow = dim(dat$XX)[2] + 1, ncol = NCOL(dat$proportion)) # include intercept

    if (proportion.model) {
      proportion <- matrix(0, nrow = dim(dat$XX)[1], ncol = NCOL(dat$proportion))
      if (dirichlet) { # the 2 cases can be wrapped into a function
        for (l in 1:(L - 1)) { # using the last cell type as reference
          proportion[, l] <- exp(cbind(1, dat$XX[, , l]) %*% zetas.current[, l]) /
            (1 + rowSums(sapply(1:(L - 1), function(xx) {
              exp(cbind(1, dat$XX[, , xx]) %*% zetas.current[, xx])
            })))
        }
        proportion[, L] <- 1 - rowSums(proportion[, -L])
      } else { # use log-link for concentration parameters without reference category
        alphas <- matrix(NA, nrow = NROW(proportion), ncol = NCOL(proportion))
        for (ll in 1:L) {
          alphas[, ll] <- exp(cbind(1, dat$XX[, , ll]) %*% zetas.current[, ll])
        }
        proportion <- apply(alphas, 2, function(xx) {
          xx / rowSums(alphas)
        })
      }
    } else {
      proportion <- dat$proportion
    }

    betas <- mu.current / gamma(1 + 1 / kappas)
    weibull.S <- matrix(nrow = n, ncol = 3)
    for (l in 1:L) {
      weibull.S[, l] <- exp(-(dat$survObj$time / betas[, l])^kappas)
    }

    xi <- c(0, 0, 0)
    thetas <- exp(dat$x0 %*% xi)
  }

  ## define output variables
  xi.mcmc <- matrix(0, nrow = 1 + nIter, ncol = NCOL(dat$x0))
  xi.mcmc[1, ] <- xi
  kappas.mcmc <- c(kappas, rep(0, nIter))
  phi.mcmc <- c(phi, rep(0, nIter))
  tauSq.mcmc <- c(tauSq, rep(0, nIter))
  # tauSq.mcmc <- matrix(0, nrow = 1+nIter, ncol = 3)
  wSq.mcmc <- c(wSq, rep(0, nIter))
  vSq.mcmc <- matrix(0, nrow = 1 + nIter, ncol = 2)
  betas.mcmc <- matrix(0, nrow = 1 + nIter, ncol = NCOL(dat$proportion) * dim(dat$XX)[2])
  betas.mcmc[1, ] <- as.vector(betas.current)
  zetas.mcmc <- matrix(0, nrow = 1 + nIter, ncol = NCOL(dat$proportion) * (dim(dat$XX)[2] + 1))
  zetas.mcmc[1, ] <- as.vector(zetas.current)

  if (method == "MLE") {
    # X0 <- dat$x0
    # X <- dat$XX
    # survObj <- dat$survObj
    # proportion <- dat$proportion

    globalvariable <- list(
      proportion = dat$proportion,
      proportion.model = proportion.model,
      X0 = dat$x0,
      X = dat$XX,
      survObj = dat$survObj
    )
    list2env(globalvariable, .GlobalEnv)
    if (proportion.model) {
      initial_values <- c(kappas, xi, as.vector(betas.current), phi, as.vector(zetas.current[, -L]))
      # initial_values <- c(dat$kappas, 3, dat$xi, as.vector(dat$betas), as.vector(dat$zetas[, -L]))
    } else {
      initial_values <- c(kappas, xi, as.vector(betas.current))
    }
    bound.l <- rep(-5, length(initial_values))
    bound.l[1] <- 1e-2
    bound.u <- rep(5, length(initial_values))
    bound.u[1] <- 10
    if (proportion.model) {
      bound.l[1 + NCOL(X0) + p * L + 1] <- 1e-2
      bound.u[1 + NCOL(X0) + p * L + 1] <- 500
    }

    mle <- nlminb(
      start = initial_values,
      objective = loglikelihood,
      lower = bound.l,
      upper = bound.u
    )
    par <- mle$par
    kappas.mcmc <- par[1]
    xi.mcmc <- par[1 + 1:NCOL(X0)]
    betas.mcmc <- matrix(par[1 + NCOL(X0) + 1:(p * L)], ncol = L)
    if (proportion.model) {
      phi.mcmc <- par[1 + NCOL(X0) + p * L + 1]
      zetas.mcmc <- matrix(par[-c(1:(1 + NCOL(X0) + p * L + 1))], nrow = p + 1)
    }
  } else {
    globalvariable <- list(
      dat = dat,
      proportion.model = proportion.model,
      dirichlet = dirichlet,
      proportion = proportion,
      kappaPrior = kappaPrior,
      kappas = kappas,
      kappaA = hyperpar$kappaA, kappaB = hyperpar$kappaB,
      lambdas = lambdas,
      weibull.S = weibull.S,
      betas.current = betas.current,
      mu.current = mu.current,
      # thetas = thetas,
      # xi = xi,
      # phi = phi,
      zetas.current = zetas.current,
      vSq = vSq,
      wSq = wSq, w0Sq = w0Sq,
      tauSq = tauSq,
      L = L
    )
    list2env(globalvariable, .GlobalEnv)

    ## MCMC iterations
    for (m in 1:nIter) {
      if (m %% tick == 0) {
        cat("MCMC iteration:", m, "\n")
      }

      # new.env()

      ## update \xi's in cure fraction
      xi.mcmc.internal <- arms(
        y.start = xi,
        myldens = logpost_xi2,
        indFunc = convex_support,
        n.sample = 5
      )
      ## n.sample = 20 will result in less variation
      xi <- colMeans(xi.mcmc.internal) # [-c(1:(nrow(xi.mcmc.internal)/2)),])
      xi <- sapply(xi, function(xx) {
        min(abs(xx), 3 - 0.1) * sign(xx)
      })
      # xi <- xi.mcmc.internal
      xi.mcmc[1 + m, ] <- xi
      globalvariable <- list(xi = xi)
      list2env(globalvariable, .GlobalEnv)

      vSq[-1] <- sampleV(2, hyperpar$vA, hyperpar$vB, xi)
      vSq.mcmc[1 + m, ] <- vSq[-1]
      globalvariable <- list(vSq = vSq)
      list2env(globalvariable, .GlobalEnv)

      thetas <- exp(dat$x0 %*% xi)
      globalvariable <- list(thetas = thetas)
      list2env(globalvariable, .GlobalEnv)

      ## update \zeta_l of p_l in non-cure fraction; the last cell type as reference
      if (proportion.model) {
        ## update phi in Dirichlet measurement error model
        phi.mcmc.internal <- arms(
          y.start = phi,
          myldens = logpost_phi,
          indFunc = convex_support_phi,
          n.sample = 10
        )
        ## n.sample = 20 will result in less variation
        # phi <- mean(phi.mcmc.internal)#[-c(1:(length(phi.mcmc.internal)/2))])#median
        phi <- median(phi.mcmc.internal[-c(1:(length(phi.mcmc.internal) / 2))]) # median
        phi <- median(c(phi, 0.1, 200 - 0.1))
        phi.mcmc[1 + m] <- phi
        globalvariable <- list(phi = phi)
        list2env(globalvariable, .GlobalEnv)

        for (l in 1:ifelse(dirichlet, L - 1, L)) {
          for (j in 1:(p + 1)) {
            # if (j == 1) wSq <- 10 # Fixed value 10 seems not good for Dirichlet's common parametrization
            # globalVariables(c("l", "j")) ## This is not safe
            globalvariableJL <- list(j = j, l = l)
            list2env(globalvariableJL, .GlobalEnv)

            zetas.mcmc.internal <- arms(
              y.start = zetas.current[j, l],
              myldens = logpost_zeta_jl,
              indFunc = convex_support,
              n.sample = 20
            )
            ## n.sample = 20 will result in less variation
            # zeta.l.new <- mean(zetas.mcmc.internal)#[-c(1:(length(zetas.mcmc.internal)/2))]) #median
            zeta.l.new <- median(zetas.mcmc.internal[-c(1:(length(zetas.mcmc.internal) / 2))]) # median
            # zeta.l.new <- zetas.mcmc.internal
            zetas.current[j, l] <- min(abs(zeta.l.new), 3 - 0.1) * sign(zeta.l.new)
            # wSq <- sampleW(1, hyperpar$wA, hyperpar$wB, zetas.current[j,l])
            globalvariable <- list(zetas.current = zetas.current)
            list2env(globalvariable, .GlobalEnv)
          }
          if (dirichlet) {
            for (ll in 1:(L - 1)) { # the l-th subtype proportion is updated, and the all composition should be updated
              proportion[, ll] <- exp(cbind(1, dat$XX[, , ll]) %*% zetas.current[, ll]) /
                (1 + rowSums(sapply(1:(L - 1), function(xx) {
                  exp(cbind(1, dat$XX[, , xx]) %*% zetas.current[, xx])
                })))
            }
            proportion[, L] <- 1 - rowSums(proportion[, -L])
          } else { # use log-link for concentration parameters without reference category
            alphas <- matrix(NA, nrow = NROW(proportion), ncol = NCOL(proportion))
            for (ll in 1:L) {
              alphas[, ll] <- exp(cbind(1, dat$XX[, , ll]) %*% zetas.current[, ll])
            }
            proportion <- apply(alphas, 2, function(xx) {
              xx / rowSums(alphas)
            })
          }
          globalvariable <- list(proportion = proportion)
          list2env(globalvariable, .GlobalEnv)
        }
        zetas.mcmc[1 + m, ] <- as.vector(zetas.current)

        ## update wSq, the variance of zetas
        if (w0Sampler == "IGamma") {
          w0Sq <- sampleW(1, hyperpar$w0A, hyperpar$w0B, zetas.current[1, ])
        }
        wSq <- sampleW(1, hyperpar$wA, hyperpar$wB, zetas.current[-1, ])
        wSq.mcmc[1 + m] <- wSq
        globalvariable <- list(wSq = wSq, w0Sq = w0Sq)
        list2env(globalvariable, .GlobalEnv)
      }

      ## update kappa in noncure fraction
      if (kappaSampler == "arms") {
        kappas.mcmc.internal <- arms(
          y.start = kappas,
          myldens = logpost_kappas,
          indFunc = convex_support_kappas,
          n.sample = 20
        )
      } else {
        kappas.mcmc.internal <- slice(
          n = 20, init_theta = kappas, TARGET = logpost_kappas,
          w = 1, max_steps = 10, k_product = 1, type = "log"
        ) # , lower=-5, upper=5)
      }
      ## n.sample = 20 will result in less variation
      # kappas <- mean(kappas.mcmc.internal)#[-c(1:(length(kappas.mcmc.internal)/2))])#median
      kappas <- median(kappas.mcmc.internal[-c(1:(length(kappas.mcmc.internal) / 2))])
      kappas <- median(c(kappas, 0.1 + 0.1, 10 - 0.1))
      kappas.mcmc[1 + m] <- kappas
      # lambdas <- mu.current / gamma(1 + 1 / kappas)
      globalvariable <- list(kappas = kappas) # , lambdas = lambdas)
      list2env(globalvariable, .GlobalEnv)

      ## update \beta_jl of S_l(t) in non-cure fraction
      for (l in 1:L) {
        for (j in 1:p) {
          # globalVariables(c("l", "j")) ## This is not safe
          globalvariableJL <- list(j = j, l = l)
          list2env(globalvariableJL, .GlobalEnv)

          betas.mcmc.internal <- arms(
            y.start = betas.current[j, l],
            myldens = logpost_beta_jl,
            indFunc = convex_support,
            # y.start = betas.current[, l],
            # myldens = logpost_beta_l,
            n.sample = 20
          )
          ## n.sample = 20 will result in less variation
          # beta.l.new <- mean(betas.mcmc.internal)#[-c(1:(length(betas.mcmc.internal)/2))]) #median
          beta.l.new <- median(betas.mcmc.internal[-c(1:(length(betas.mcmc.internal) / 2))]) # median
          # beta.l.new <- betas.mcmc.internal
          betas.current[j, l] <- min(abs(beta.l.new), 3 - 0.1) * sign(beta.l.new)
          ## force irrelevant betas to be zero
          # betas.current[j, l] <- ifelse(dat$betas[j, l] == 0, 0, betas.current[j, l])

          # beta.l.new <- colMeans(betas.mcmc.internal[-c(1:(nrow(betas.mcmc.internal)/2)),])
          # betas.current[,l] <- sapply(beta.l.new, function(xx){min(abs(xx), 5-0.1) * sign(xx)} )
          # tauSq <- sampleTau(1, hyperpar$tauA, hyperpar$tauB, betas.current[j, l])

          # globalvariable <- list(betas.current = betas.current)
          # list2env(globalvariable, .GlobalEnv)

          mu.current[, l] <- exp(dat$XX[, , l] %*% betas.current[, l]) # should this and following be out of this for-loop-j?
          lambdas[, l] <- mu.current[, l] / gamma(1 + 1 / kappas)

          weibull.S[, l] <- exp(-(dat$survObj$time / lambdas[, l])^kappas)
          globalvariable <- list(
            betas.current = betas.current # ,
            # mu.current = mu.current, # no need to enter GlobalEnv, since it's updated in logpost_beta_jl()
            # weibull.S = weibull.S
          )
          list2env(globalvariable, .GlobalEnv)
        }
        # mu.current[, l] <- exp(dat$XX[, , l] %*% betas.current[, l])
        # lambdas[, l] <- mu.current[, l] / gamma(1 + 1 / kappas) # lambdas is a parameter in WEI3 distr.
        # weibull.S[, l] <- exp(-(dat$survObj$time / lambdas[, l])^kappas)
        # # tauSq[l] <- sampleTau(1, hyperpar$tauA, hyperpar$tauB, betas.current[, l])
        globalvariable <- list(
          mu.current = mu.current,
          weibull.S = weibull.S
        )
        list2env(globalvariable, .GlobalEnv)
      }
      betas.mcmc[1 + m, ] <- as.vector(betas.current)
      globalvariable <- list(lambdas = lambdas, mu.current = mu.current, weibull.S = weibull.S) # do we need mu&S her?
      list2env(globalvariable, .GlobalEnv)

      ## update tauSq, the variance of betas
      # tauSq <- sampleTau(3, hyperpar$tauA, hyperpar$tauB, betas.current)
      # tauSq.mcmc[1+m, ] <- tauSq
      tauSq <- sampleTau(1, hyperpar$tauA, hyperpar$tauB, betas.current)
      tauSq.mcmc[1 + m] <- tauSq
      globalvariable <- list(tauSq = tauSq)
      list2env(globalvariable, .GlobalEnv)
    }
  }
  cat("... Done!\n")

  ret <- list(input = list(), output = list(), call = cl)
  class(ret) <- "GPTCM"

  ret$input$dirichlet <- dirichlet
  ret$input$nIter <- nIter
  ret$input$burnin <- burnin
  ret$input$hyperpar <- hyperpar

  # survival predictions based on posterior mean
  if (method == "Bayes") {
    ret$output$posterior <- list(
      xi = colMeans(xi.mcmc[-c(1:(nIter / 2)), ]),
      kappas = mean(kappas.mcmc[-c(1:(nIter / 2))]),
      phi = mean(phi.mcmc[-c(1:(nIter / 2))]),
      betas = matrix(colMeans(betas.mcmc[-c(1:(nrow(betas.mcmc) / 2)), ]), ncol = L),
      thetas = exp(dat$x0 %*% xi)
    )
  }

  ret$output$mcmc <- list(
    xi = xi.mcmc,
    kappas = kappas.mcmc,
    phi = phi.mcmc,
    betas = betas.mcmc,
    zetas = zetas.mcmc,
    tauSq = tauSq.mcmc,
    wSq = wSq.mcmc,
    vSq = vSq.mcmc
  )

  return(ret)
}
