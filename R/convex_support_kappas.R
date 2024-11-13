#' @title Convex support
#'
#' @description
#' This is an indicator function of the convex support of Weibull's kappa
#'
#' @name convex_support_kappas
#'
#' @param x TBA
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
convex_support_kappas <- function(x) {
  (min(x) > 0.1) * (max(x) < 10)
}
