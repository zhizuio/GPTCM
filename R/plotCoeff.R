#' @title Plots for estimated coefficients
#'
#' @description
#' create nice plots for estimated coefficients and 95% credible intervals
#'
#' @name plotCoeff
#'
#' @param dat
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
plotCoeff <- function(dat, datMCMC, estimator = "beta", ...) {
  #n <- dim(dat$XX)[1]
  p <- dim(dat$XX)[2]
  L <- dim(dat$XX)[3]
  #nIter <- datMCMC$input$nIter
  burnin <- datMCMC$input$burnin

  if (estimator == "beta") {
    betas.mcmc <- datMCMC$output$mcmc$betas
    # Final estimates
    beta.mcmc0 <- betas.mcmc[-c(1:burnin), ]
    beta.true <- dat$betas[p:1, ]
    beta_min <- min(beta.mcmc0, beta.true)
    beta_max <- max(beta.mcmc0, beta.true)

    # pdf(paste0("gptcm_betaHat_Sigma",Sigma,"_M",M,".pdf"), height = 4, width = 8)
    layout(matrix(1:L, nrow = 1))
    for (l in 0:(L - 1)) {
      plotCoeff(beta.mcmc0[, l * p + 1:p][, p:1], beta.true[, l + 1],
        xlim = c(beta_min, beta_max),
        main = paste("Cell type", l),
        label.y = paste0("x", p:1)
      )
    }
    # dev.off()
  }

  if (estimator == "zeta") {
    # Final estimates
    zetas.mcmc <- datMCMC$output$mcmc$zetas[-c(1:burnin), ]
    zetas.true <- dat$zetas[(p + 1):1, ]
    zetas_min <- min(zetas.mcmc, zetas.true)
    zetas_max <- max(zetas.mcmc, zetas.true)

    # pdf(paste0("gptcm_zetaHat_Sigma",Sigma,"_M",M,".pdf"), height = 5, width = 7)
    layout(matrix(1:(L - 1), nrow = 1))
    for (l in 0:(L - 2)) {
      plotCoeff(zetas.mcmc[, l * (p + 1) + 1:(p + 1)][, (p + 1):1],
        zetas.true[, l + 1],
        xlim = c(zetas_min, zetas_max),
        main = paste("Cell type", l),
        label.y = c(paste0("x", p:1), "intecept")
      )
    }
    # dev.off()
  }
}
