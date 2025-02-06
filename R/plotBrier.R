#' @title Time-dependent Brier scores
#'
#' @description
#' Predict time-dependent Brier scores based on different survival models
#'
#' @name plotBrier
#'
#' @importFrom survival Surv coxph
#' @importFrom riskRegression Score
#' @importFrom stats median as.formula
#' @importFrom ggplot2 ggplot aes geom_step theme element_blank xlab ylab theme_bw
#' @importFrom graphics layout par abline
#' @importFrom utils globalVariables
#'
#' @param dat TBA
#'
#' @return A \code{ggplot2::ggplot} object. See \code{?ggplot2::ggplot} for more
#' details of the object.
#'
#' @examples
#'
#' x <- 1
#'
#' @export
plotBrier <- function(dat, datMCMC, 
                      time.star = NULL, 
                      xlab = "Time", 
                      ylab = "Brier score", ...) {
  n <- dim(dat$XX)[1]
  p <- dim(dat$XX)[2]
  L <- dim(dat$XX)[3]
  # nIter <- datMCMC$input$nIter
  burnin <- datMCMC$input$burnin

  # survival predictions based on posterior mean
  xi.hat <- colMeans(datMCMC$output$mcmc$xi[-c(1:burnin), ])
  betas.hat <- matrix(colMeans(datMCMC$output$mcmc$betas[-c(1:burnin), ]), ncol = L)
  kappas.hat <- mean(datMCMC$output$mcmc$kappas[-c(1:burnin)])
  thetas.hat <- exp(dat$x0 %*% xi.hat)


  # predict survival probabilities based on GPTCM
  time_eval <- sort(dat$survObj$time)
  Surv.prob <- matrix(nrow = n, ncol = length(time_eval))
  for (j in 1:length(time_eval)) {
    tmp <- 0
    for (l in 1:L) {
      mu <- exp(dat$XX[, , l] %*% betas.hat[, l])
      lambdas <- mu / gamma(1 + 1 / kappas.hat)
      weibull.S <- exp(-(time_eval[j] / lambdas)^kappas.hat)
      tmp <- tmp + dat$proportion[, l] * weibull.S
    }
    Surv.prob[, j] <- exp(-thetas.hat * (1 - tmp))
  }

  pred.prob <- 1 - Surv.prob

  # re-organize clinical variables for classical Cox and PTCM models
  x <- apply(dat$XX, c(1, 2), mean)
  colnames(x) <- paste0("x", 1:p)
  x.median <- apply(dat$XX, c(1, 2), median)
  colnames(x.median) <- paste0("x.median", 1:p)
  survObj <- data.frame(dat$survObj, x01 = dat$x0[, 2], x02 = dat$x0[, 3], x, x.median)
  fitCox.clin <- survival::coxph(survival::Surv(time, event) ~ x01 + x02, data = survObj, y = TRUE, x = TRUE)

  # pred.clin <- survival::survfit(fitCox.clin, data = survObj)
  formula.tmp <- as.formula(paste0("Surv(time, event) ~ ", paste0(colnames(x.median), collapse = "+")))
  fitCox.X.median <- survival::coxph(formula.tmp, data = survObj, y = TRUE, x = TRUE)
  # pred.X.median <- survival::survfit(fitCox.X.median, data = survObj)

  formula.tmp <- as.formula(paste0("Surv(time, event) ~ ", paste0(paste0("x", 1:p), collapse = "+")))
  fitCox.X.mean <- survival::coxph(formula.tmp, data = survObj, y = TRUE, x = TRUE)
  # formula.tmp <- as.formula(paste0("Surv(time, event) ~ x01+x02+", paste0(colnames(x.median), collapse = "+")))
  # fitCox.clin.X.median <- survival::coxph(formula.tmp, data = survObj, y=TRUE, x = TRUE)

  formula.tmp <- as.formula(paste0("Surv(time, event) ~ x01+x02+", paste0(paste0("x", 1:p), collapse = "+")))
  fitCox.clin.X.mean <- survival::coxph(formula.tmp, data = survObj, y = TRUE, x = TRUE)

  # library(miCoPTCM) # good estimation for cure fraction; same BS as Cox.clin
  suppressWarnings(
    resMY <- miCoPTCM::PTCMestimBF(Surv(time, event) ~ x01 + x02,
      data = survObj,
      varCov = matrix(0, nrow = 3, ncol = 3),
      init = rep(0, 3)
    )
  )
  Surv.PTCM <- exp(-exp(dat$x0 %*% resMY$coefficients) %*% t(resMY$estimCDF))
  predPTCM.prob <- 1 - Surv.PTCM

  g <- riskRegression::Score(
    list(
      "Cox.clin" = fitCox.clin,
      "Cox.X.mean" = fitCox.X.mean,
      "Cox.X.median" = fitCox.X.median,
      # "Cox.clin.X.median"=fitCox.clin.X.median,
      "Cox.clin.X.mean" = fitCox.clin.X.mean,
      "PTCM.clin" = predPTCM.prob,
      "GPTCM" = pred.prob
    ),
    formula = Surv(time, event) ~ 1,
    metrics = "brier", summary = "ibs",
    data = survObj,
    conf.int = FALSE, times = time_eval
  )
  g1 <- g$Brier$score
  if(!is.null(time.star)){
    g1 <- g1[g1$times <= time.star, ]
  }
  levels(g1$model)[1] <- "Kaplan-Meier"
  # utils::globalVariables(c("times", "Brier", "model"))
  # NOTE: `aes_string()` was deprecated in ggplot2 3.0.0.
  g2 <- ggplot2::ggplot(g1, aes(
    # x = "times", y = "Brier", group = "model", color = "model"
    x = times, y = Brier, group = model, color = model
  )) +
    xlab(xlab) +
    ylab(ylab) +
    geom_step(direction = "vh") + # , alpha=0.4) +
    theme_bw() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.4, 0.25),
      legend.title = element_blank()
    )

  g2
}
