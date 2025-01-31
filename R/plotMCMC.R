#' @title MCMC trace-plots
#'
#' @description
#' MCMC
#'
#' @name plotBrier
#'
#' @importFrom ggplot2 ggplot aes geom_step theme element_blank
#'
#' @param dat TBA
#' @param datMCMC TBA
#' @param estimator TBA
#' @param ... TBA
#'
#' @return A \code{ggplot2::ggplot} object. See \code{?ggplot2::ggplot} for more
#' details of the object.
#'
#' @examples
#'
#' x <- 1
#'
#' @export
plotMCMC <- function(dat, datMCMC, estimator = "xi") {
  # n <- dim(dat$XX)[1]
  p <- dim(dat$XX)[2]
  L <- dim(dat$XX)[3]
  # nIter <- datMCMC$input$nIter
  # burnin <- datMCMC$input$burnin

  if ("beta" %in% estimator) {
    betas.mcmc <- datMCMC$output$mcmc$betas

    ylabel <- paste0("expression(beta['", rep(1:p, L), ",", rep(1:L, each = p), "'])")
    layout(matrix(1:NCOL(betas.mcmc), ncol = L))
    par(mar = c(2, 4.1, 2, 2))
    for (j in 1:NCOL(betas.mcmc)) {
      plot(betas.mcmc[, j],
        type = "l", lty = 1, ylab = eval(parse(text = ylabel[j])),
        ylim = summary(c(betas.mcmc, dat$betas))[c(1, 6)]
      )
      abline(h = dat$betas[j], col = "red")
    }
  }

  if ("zeta" %in% estimator) {
    dirichlet <- datMCMC$input$dirichlet
    if (dirichlet) {
      zetas.mcmc <- datMCMC$output$mcmc$zetas
      ylabel <- paste0(
        "expression(zeta['", rep(0:p, L), ",",
        rep(1:L, each = p + 1), "'])"
      )
    } else {
      zetas.mcmc <- datMCMC$output$mcmc$zetas[, 1:((p + 1) * (L - 1))]
      ylabel <- paste0(
        "expression(zeta['", rep(0:p, L - 1), ",",
        rep(1:(L - 1), each = p + 1), "'])"
      )
    }

    ylabel <- paste0(
      "expression(zeta['", rep(0:p, ifelse(dirichlet, L, L - 1)), ",",
      rep(1:ifelse(dirichlet, L, L - 1), each = p + 1), "'])"
    )
    layout(matrix(1:NCOL(zetas.mcmc), nrow = p + 1))
    par(mar = c(2, 4.1, 2, 2))
    for (j in 1:NCOL(zetas.mcmc)) {
      plot(zetas.mcmc[, j],
        type = "l", lty = 1, ylab = eval(parse(text = ylabel[j])),
        ylim = summary(c(zetas.mcmc, dat$zetas))[c(1, 6)]
      )
      abline(h = dat$zetas[j], col = "red")
    }
  }

  if ("xi" %in% estimator) {
    p <- dim(datMCMC$output$mcmc$xi)[2]
    xi.mcmc <- datMCMC$output$mcmc$xi

    ylabel <- paste0("expression(xi[", 0:2, "])")
    layout(matrix(1:p, ncol = 1))
    par(mar = c(2, 4.1, 2, 2))
    for (j in 1:p) {
      plot(xi.mcmc[, j],
        type = "l", lty = 1, ylab = eval(parse(text = ylabel[j])),
        ylim = summary(c(xi.mcmc, dat$xi))[c(1, 6)]
      )
      abline(h = dat$xi[j], col = "red")
    }
  }


  if (any(estimator %in% c("kappa", "tau", "w", "v", "phi"))) {
    layout(matrix(1:length(estimator), ncol = 1))
    par(mar = c(2, 4.1, 2, 2))
    if ("kappa" %in% estimator) {
      kappas.mcmc <- datMCMC$output$mcmc$kappas
      plot(kappas.mcmc,
        type = "l", lty = 1,
        ylab = expression(kappa), xlab = "MCMC iteration",
        ylim = summary(c(kappas.mcmc, dat$kappas))[c(1, 6)]
      )
      abline(h = dat$kappas, col = "red")
    }

    if ("tau" %in% estimator) {
      tauSq.mcmc <- datMCMC$output$mcmc$tauSq
      plot(tauSq.mcmc,
        type = "l", lty = 1,
        ylab = expression(tau^2), xlab = "MCMC iteration",
        ylim = summary(as.vector(tauSq.mcmc))[c(1, 6)]
      )
    }

    if ("w" %in% estimator) {
      wSq.mcmc <- datMCMC$output$mcmc$wSq
      plot(wSq.mcmc,
        type = "l", lty = 1,
        ylab = expression(w^2), xlab = "MCMC iteration",
        ylim = summary(as.vector(wSq.mcmc))[c(1, 6)]
      )
    }

    if ("v" %in% estimator) {
      vSq.mcmc <- datMCMC$output$mcmc$vSq
      plot(vSq.mcmc,
        type = "l", lty = 1,
        ylab = expression(v^2), xlab = "MCMC iteration",
        ylim = summary(as.vector(vSq.mcmc))[c(1, 6)]
      )
    }

    if ("phi" %in% estimator) {
      phi.mcmc <- datMCMC$output$mcmc$phi
      plot(phi.mcmc,
        type = "l", lty = 1,
        ylab = expression(phi), xlab = "MCMC iteration",
        ylim = summary(as.vector(phi.mcmc))[c(1, 6)]
      )
    }
  }
}
