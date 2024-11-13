#' @title Target density
#'
#' @description
#' Predefined target density for S_pop in GTPCM
#'
#' @name target
#'
#' @param x TBA
#' @param theta TBA
#' @param proportion TBA
#' @param mu TBA
#' @param kappas TBA
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
target <- function(x, theta, proportion, mu, kappas) {
  # survival.function <- exp( -theta * (1 - sum(proportion * exp(-lambdas*x^kappa)))  )
  # # un-normalized pdf
  # pdf <- survival.function * theta *
  #   kappa*x^(kappa-1) *
  #   sum(proportion * lambdas * exp(-lambdas*x^kappa))
  lambdas <- mu / gamma(1 + 1 / kappas)
  survival.function <- exp(-(x / lambdas)^kappas)
  # improper pdf
  pdf <- exp(-theta * (1 - sum(proportion * survival.function))) *
    theta *
    sum(proportion * kappas / lambdas * 
          (x / lambdas)^(kappas - 1) * 
          exp(-(x / lambdas)^kappas))

  ## exponetial survival
  # survival.function <- exp( -lambdas * x  )
  # pdf <- exp(-theta * (1-survival.function)) *
  #   theta * lambdas * exp(-lambdas*x)

  ## Weibull 1
  # survival.function <- exp( - (x/lambdas)^kappa  )
  # pdf <- exp(-theta * (1-survival.function)) *
  #   theta * kappa / lambdas * (x/lambdas)^(kappa - 1) * exp(-(x/lambdas)^kappa)

  ## Weibull 2
  # survival.function <- exp( - lambdas*x^kappa  )
  # pdf <- exp(-theta * (1-survival.function)) *
  #   theta * lambdas * x^(kappa - 1) * exp(-lambdas*x^kappa)

  ## Weibull 3
  # betas <- lambdas / gamma(1 + 1/kappa)
  # survival.function <- exp( -(x/betas)^kappa )
  # pdf <- exp( -theta * (1 - survival.function) ) *
  #   theta *
  #   kappa/betas * (x/betas)^(kappa-1) * exp(-(x/betas)^kappa)
  return(pdf)
}
