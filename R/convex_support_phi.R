#' @title Convex support
#'
#' @description
#' This is an indicator function of the convex support of Dirichlet's dispersion phi
#'
#' @name convex_support_phi
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
convex_support_phi <- function(x) {
  (min(x) > 0.1) * (max(x) < 200)
}
