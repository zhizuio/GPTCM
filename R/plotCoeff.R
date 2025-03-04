#' @title Plots for estimated coefficients
#'
#' @description
#' create nice plots for estimated coefficients and 95% credible intervals
#'
#' @name plotCoeff
#'
#' @importFrom stats quantile
#' @importFrom graphics axis points arrows legend
#' @importFrom utils capture.output
#' @importFrom scales alpha
#'
#' @param dat TBA
#' @param datMCMC TBA
#' @param estimator TBA
#' @param ... TBA
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
plotCoeff <- function(dat, datMCMC, 
                      estimator = "beta", 
                      label.y = NULL,
                      legend.labs = NULL, ...) {
  # n <- dim(dat$XX)[1]
  p <- dim(dat$XX)[2]
  L <- dim(dat$XX)[3]
  # nIter <- datMCMC$input$nIter
  burnin <- datMCMC$input$burnin
  
  legend0 <- ifelse(is.null(legend.labs), 0, 1)
  
  if (estimator == "beta") {
    betas.mcmc <- datMCMC$output$mcmc$betas[-c(1:burnin), ]
    # Final estimates
    betas.true <- dat$betas[p:1, ]
    betas_min <- min(betas.mcmc, betas.true)
    betas_max <- max(betas.mcmc, betas.true)
    
    if(is.null(label.y)){
      label.y <- paste0("x", 1:p)
    }

    # pdf(paste0("gptcm_betaHat_Sigma",Sigma,"_M",M,".pdf"), height = 4, width = 8)
    layout(matrix(1:(L + legend0), nrow = 1))
    for (l in 0:(L - 1)) {
      plotCoeff0(betas.mcmc[, l * p + 1:p][, p:1], betas.true[, l + 1],
        xlim = c(betas_min, betas_max),
        legend.labs = legend.labs,
        main = paste("Cell type", l + 1),
        label.y = label.y[p:1]
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
    dirichlet <- datMCMC$input$dirichlet
    
    if(is.null(label.y)){
      label.y <- c("intecept", paste0("x", 1:p))
    }
    
    # pdf(paste0("gptcm_zetaHat_Sigma",Sigma,"_M",M,".pdf"), height = 5, width = 7)
    layout(matrix(1:(ifelse(dirichlet, L, L - 1) + legend0), nrow = 1))
    for (l in 0:ifelse(dirichlet, L - 1, L - 2)) {
      plotCoeff0(zetas.mcmc[, l * (p + 1) + 1:(p + 1)][, (p + 1):1],
        zetas.true[, l + 1],
        xlim = c(zetas_min, zetas_max),
        legend.labs = legend.labs, 
        pars.name = "zeta",
        main = paste("Cell type", l + 1),
        label.y = label.y[(p + 1):1]
      )
    }
    # dev.off()
  }
  
  if (!is.null(legend.labs)) {
    par(mar = c(5.1, 0, 4.1, 0))
    plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
    legend("topleft",
           bty = "n",
           legend = legend.labs,
           pch = c(4, 5, NA), 
           lty = c(NA, NA, 1),
           col = scales::alpha(c("red", "green", "green"), 0.7)
    )
    par(mar = c(5.1, 4.1, 4.1, 2.1))
  }
}

## create nice plots for estimated coefficients and 95% credible intervals

plotCoeff0 <- function(pars.mcmc, true.param,
                       credible.interval = TRUE,
                       legend.labs = legend.labs, 
                       pars.name = "beta",
                       label.y = NULL,
                       main = "Type 1", ...) {
  if (is.null(label.y)) {
    label.y <- paste0("x", 1:NCOL(pars.mcmc))
  }

  pars_mean <- apply(pars.mcmc, 2, mean)
  pars_Q2.5 <- apply(pars.mcmc, 2, function(xx) {
    quantile(xx, 0.025)
  })
  pars_Q97.5 <- apply(pars.mcmc, 2, function(xx) {
    quantile(xx, 0.975)
  })

  # if(is.null(pars_min))
  #  pars_min <- min(pars_Q2.5, true.param, na.rm = TRUE)
  # if(is.null(pars_max))
  #  pars_max <- max(pars_Q97.5, true.param, na.rm = TRUE)

  # plotting true values
  plot(
    x = true.param,
    y = c(seq_along(true.param)),
    pch = 4,
    col = scales::alpha("red", 0.7),
    xlab = "Effect",
    ylab = "",
    # xlim = c(pars_min, pars_max),
    yaxt = "n",
    main = main, ...
  )

  axis(2, at = 1:ncol(pars.mcmc), labels = label.y, las = 1)
  abline(v = 0, col = "gray95", lty = 5)

  # plotting posterior mean
  points(x = pars_mean, y = 1:ncol(pars.mcmc), pch = 5, col = scales::alpha("green", 0.7))
  if (is.null(legend.labs)) {
    if (pars.name == "beta") {
      legend("topleft",
             legend = c(expression(beta), expression(hat(beta))),
             pch = c(4, 5), col = scales::alpha(c("red", "green"), 0.7)
      )
    }
    if (pars.name == "zeta") {
      legend("topleft",
             legend = c(expression(zeta), expression(hat(zeta))),
             pch = c(4, 5), col = scales::alpha(c("red", "green"), 0.7)
      )
    }
  }

  # plotting 95% credible interval
  if (credible.interval) {
    # suppressWarnings(
    arrows(
      y0 = seq_along(pars_mean), x0 = pars_Q2.5,
      y1 = seq_along(pars_mean), x1 = pars_Q97.5,
      length = 0.02, angle = 90, col = scales::alpha("green", 0.7),
      code = 3
    )
    # )
  }
}
