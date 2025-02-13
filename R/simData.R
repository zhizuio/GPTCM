#' @title Simulate data
#'
#' @description
#' Simulate survival data. (This requires external function 'metropolis_sampler')
#'
#' @name simData
#'
#' @importFrom stats rbinom rnorm runif rexp
#'
#' @param n TBA
#' @param p TBA
#' @param L TBA
#' @param Sigma TBA
#' @param kappas value of the Weibull's shape parameter
#' @param proportion.model One of \code{c("alr", "cloglog", "log", "dirichlet")}
#' @param phi TBA
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
simData <- function(n = 200, p = 10, L = 3,
                    Sigma = 0, kappas = 2,
                    proportion.model = "dirichlet",
                    phi = 10) {
  ## predefined functions
  Expo <- function(times, surv) {
    z1 <- -log(surv[1])
    t1 <- times[1]
    lambda <- z1 / (t1)
    list(rate = lambda)
  }

  ## effects; need to change values of zetas
  beta1 <- c(-1, -.5, .5, .2, -.5, 0, 0, 0, 0, 0)
  # beta2 <- c(0, -.3, .8, 0, .8, -.4, 0, 0, 0, 0)
  beta2 <- c(0, -.3, -.8, 0, .8, -.4, 0, 0, 0, 0)
  beta3 <- c(1, 0, -0.4, -1.5, 0, 0, 0.8, 0, 0, 0)

  zeta1 <- c(0.7, -0.7, 0.5, -0.5, 1, 0, 0, 0, 0, 0)
  zeta2 <- c(-0.5, 0.5, 0, 1, 0, -1, 0, 0, 0, 0)
  zeta3 <- c(0, 0, 1, -0.5, -0.7, 0, 1, 0, 0, 0)

  if (p < 10) {
    beta1 <- beta1[1:p]
    beta2 <- beta2[1:p]
    beta3 <- beta3[1:p]
    zeta1 <- zeta1[1:p]
    zeta2 <- zeta2[1:p]
    zeta3 <- zeta3[1:p]
  }

  ## covariates
  # means <- rep(0, p)
  if (Sigma[1] == 0) {
    rho0 <- 0.3 # correlation for each gene between cell types
    rho1 <- 0.1 # correlation for in 1st cell type
    rho2 <- 0.2 # correlation for in 2nd cell type
    rho3 <- 0.08 # correlation for in 3rd cell type
    Sigma1 <- Sigma2 <- Sigma3 <- diag(1, p)
    for (j in 1:p) {
      for (jj in 1:p) {
        Sigma1[j, jj] <- rho1^abs(jj - j)
        Sigma2[j, jj] <- rho2^abs(jj - j)
        Sigma3[j, jj] <- rho3^abs(jj - j)
      }
    }
    # Sigma1[lower.tri(Sigma1)] <- Sigma1[upper.tri(Sigma1)]
    # Sigma2[lower.tri(Sigma2)] <- Sigma2[upper.tri(Sigma2)]
    # Sigma3[lower.tri(Sigma3)] <- Sigma3[upper.tri(Sigma3)]

    # #set.seed(123)
    # x1 <- scale(MASS::mvrnorm(n, means, Sigma1))
    # x2 <- scale(MASS::mvrnorm(n, means, Sigma2))
    # x3 <- scale(MASS::mvrnorm(n, means, Sigma3))
    # simulate (x1, x2, x3) simultaneously with a big Sigma
    Sigma <- Matrix::bdiag(Sigma1, Sigma2, Sigma3)
    off.diag1 <- list(1:(2 * p), p + 1:(2 * p))
    off.diag2 <- list(1:p, 2 * p + 1:p)
    for (j in 1:(2 * p)) {
      Sigma[off.diag1[[1]][j], off.diag1[[2]][j]] <-
        Sigma[off.diag1[[2]][j], off.diag1[[1]][j]] <- rho0
    }
    for (j in 1:p) {
      Sigma[off.diag2[[1]][j], off.diag2[[2]][j]] <-
        Sigma[off.diag2[[2]][j], off.diag2[[1]][j]] <- rho0
    }
  }
  # Matrix::image(Sigma, sub = "", xlab = "X1     X2     X3", ylab = "")
  # pdf("X_cov.pdf", height = 3, width = 3)
  # par(mar = c(2, 2, 2, 2))
  # graphics::image(data.matrix(Sigma), col = gray.colors(33, 0, 1, rev = T), axes = FALSE)
  # axis(1, at = c(0.33/2, 0.33+0.33/2, 1-0.33/2), labels = paste0("X",1:3), tick=F)
  # axis(2, at = c(0.33/2, 0.33+0.33/2, 1-0.33/2), labels = paste0("X",1:3), tick=F)
  # abline(h = c(0.33, 0.67), v = c(0.33, 0.67), col = "gray90", lty = 3, lwd = 1)
  # box(lwd = 0.5)
  # dev.off()


  X <- scale(MASS::mvrnorm(n, rep(0, p * L), Sigma))
  x1 <- X[, 1:p]
  x2 <- X[, p + 1:p]
  x3 <- X[, 2 * p + 1:p]
  XX <- array(X, dim = c(n, p, L))

  x01 <- rbinom(n, 1, 0.5)
  x02 <- rnorm(n)
  x0 <- cbind(1, x01, x02) # mandatory covariates for cured fraction
  # xi <- c(-0.5, 1.5, 0.2) # effects of mandatory covariates for cured fraction
  xi <- c(0.3, 0.5, -0.9)
  beta0 <- rep(0, 3) # c(-1.5, 0.5, -1.) # intercepts in different cell types linked to survival
  zeta0 <- c(-.5, -1, 1.2) # intercepts in different cell types linked to proportions
  # zeta0 <- rep(0, 3)
  zetas <- rbind(zeta0, cbind(zeta1, zeta2, zeta3))

  # censoring function
  # - follow-up time 12 to 60 months
  # - administrative censoring: uniform data entry (cens1)
  # - loss to follow-up: exponential, 20% loss in 60 months (cens2)
  ACT <- 1
  FUT <- 5
  cens.start <- FUT
  cens.end <- ACT + FUT
  cens1 <- runif(n, cens.start, cens.end)
  loss <- Expo(times = 5, surv = 0.8)
  cens2 <- rexp(n, rate = loss$rate)
  cens <- pmin(cens1, cens2)

  ## survival times
  # p0 <- 1 / (1 + exp(-x0 %*% xi))
  thetas <- exp(x0 %*% xi)
  cure <- exp(-thetas) # cure probabilities

  mu1 <- exp(beta0[1] + x1 %*% matrix(beta1, ncol = 1))
  mu2 <- exp(beta0[2] + x2 %*% matrix(beta2, ncol = 1))
  mu3 <- exp(beta0[3] + x3 %*% matrix(beta3, ncol = 1))
  mu0 <- cbind(mu1, mu2, mu3)

  # simulate proportions from cloglog-link function; this is a bit model misspecification
  if (proportion.model == "cloglog") {
    p1 <- 1 - exp(-exp(zeta0[1] + x1 %*% matrix(zeta1, ncol = 1)))
    p2 <- 1 - exp(-exp(zeta0[2] + x2 %*% matrix(zeta2, ncol = 1)))
    p3 <- 1 - exp(-exp(zeta0[3] + x3 %*% matrix(zeta3, ncol = 1)))
    proportion <- cbind(p1, p2, p3)
    proportion <- t(apply(proportion, 1, function(pp) pp / sum(pp)))
  }
  # simulate proportions from log-link function; this is a bit model misspecification
  if (proportion.model == "log") {
    p1 <- exp(zeta0[1] + x1 %*% matrix(zeta1, ncol = 1))
    p2 <- exp(zeta0[2] + x2 %*% matrix(zeta2, ncol = 1))
    p3 <- exp(zeta0[3] + x3 %*% matrix(zeta3, ncol = 1))
    # p3 <- 1 - p1 - p2 # We cannot use this, since p1+p2 can be > 1
    proportion <- cbind(p1, p2, p3)
    proportion <- t(apply(proportion, 1, function(pp) pp / sum(pp)))
  }

  ## simulate proportions from Dirichlet distribution (n, alpha=1:L)
  if (proportion.model == "dirichlet") {
    alpha1 <- exp(zeta0[1] + x1 %*% matrix(zeta1, ncol = 1))
    alpha2 <- exp(zeta0[2] + x2 %*% matrix(zeta2, ncol = 1))
    alpha3 <- exp(zeta0[3] + x3 %*% matrix(zeta3, ncol = 1))
    alphas <- cbind(alpha1, alpha2, alpha3)
    #alpha0 <- alpha1 + alpha2 + alpha3
    #proportion <- cbind(alpha1 / alpha0, alpha2 / alpha0, alpha3 / alpha0)
    
    ## Generate n-individual Dirichlet proportions
    proportion <- matrix(nrow = n, ncol = L)
    for (i in 1:n) {
      proportion[i, ] <- sapply(1:L, function(l) rgamma(1, alphas[i, l]))
    }
    
    # proportion <- matrix(rgamma(L * n, t(1:L)), ncol = L, byrow=TRUE)
    proportion <- proportion / rowSums(proportion)
  }

  ## the following is via logit/alr-link function
  if (proportion.model == "alr") {
    proportion <- tmp <- matrix(nrow = n, ncol = L)
    for (l in 1:(L - 1)) {
      tmp[, l] <- exp(cbind(1, XX[, , l]) %*% zetas[, l])
    }
    proportion[, L] <- 1 / (1 + rowSums(tmp[, -L]))
    for (l in 1:(L - 1)) {
      proportion[, l] <- tmp[, l] / (1 + rowSums(tmp[, -L]))
    }
    
    alphas <- proportion * phi
    for (i in 1:n) {
      proportion[i, ] <- sapply(1:L, function(l) rgamma(1, alphas[i, l]))
    }
    proportion <- proportion / rowSums(proportion)
    
    # the last category is as reference, no need data for zetas[,L]
    zetas <- zetas[, -L]
  }

  ## fixed proportions for debug estimation of other parameters
  # proportion <- matrix(c(0.13, 0.53, 0.34), nrow = n, ncol = 3, byrow = TRUE)
  colnames(proportion) <- paste0("p", 1:3)

  # kappas <- 2 # 0.9
  ## simulate censoring times
  ## to be updated: dim(mu0)=n x 3; how can dim(cens)=n x 3? No, not use rWEI3 but M-H sampler to get times
  # browser()
  # cens <- as.vector( sapply(1:n, function(i) gamlss.dist::rWEI3(1, mu=mu0[i], sigma=kappas)) )

  ## simulate event times by M-H algorithm
  U <- runif(n)
  T.star <- cens
  T.star[U <= cure] <- Inf
  accepted <- numeric(n)

  for (i in which(U > cure)) {
    ## M-H sampler for event time
    out <- metropolis_sampler( # If the target is set as Gompertz distr., it's a bit model misspecification
      initial_value = 10,
      n = 5,
      proposal_shape = 0.9,
      proposal_scale = mean(cens),
      theta = thetas[i],
      proportion = proportion[i, ],
      mu = mu0[i, ],
      kappas = kappas,
      burnin = 100,
      lag = 10
    )
    T.star[i] <- mean(out$value)
    accepted[i] <- mean(out$accepted)
  }

  # survival object
  event <- ifelse(T.star <= cens, 1, 0) # censoring indicator
  times <- pmin(T.star, cens) # observed times
  survObj <- data.frame(event = event, time = times)

  return(list(
    survObj = survObj, accepted = accepted,
    proportion.model = proportion.model,
    proportion = proportion,
    thetas = thetas,
    kappas = kappas, mu = mu0,
    x0 = x0, # x1 = x1, x2 = x2, x3 = x3,
    XX = XX,
    xi = xi, beta0 = beta0, # zeta0=zeta0,
    betas = cbind(beta1, beta2, beta3),
    zetas = zetas
  ))
}
