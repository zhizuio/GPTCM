#' @title Convex support
#'
#' @description
#' This is an indicator function of the convex support of the target density. This is required by R-pkg 'HI'
#'
#' @name convex_support
#'
#' @param x
#'
#' @return An object of ...
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
convex_support <- function(x) {
  (min(x) > -3) * (max(x) < 3)
  # (min(x) > -5) * (max(x) < 5)
}
