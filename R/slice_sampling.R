#' @title Slice sampling
#'
#' @description
#' This file will include functions which implement (product) slice sampling on
#' a target distribution adapted from
#' (i) https://github.com/TJHeaton/NonparametricCalibration/blob/main/SliceUpdateFunsFinal.R,
#' (ii) https://www.jarad.me/teaching/2013/10/29/slice-sampling-for-multimodal-TARGETs
#'
#' We just need to call
#' slice(n, init_theta, TARGET, w, max_steps, k_product, lower, upper, type)
#' with TARGET being either:
#' 1) Function f(x) which is proportional to the TARGET density
#' 2) g(x) = log(f(x)) in which case need extra arg type = "log"
#'
#' @name slice
#'
#' @importFrom stats runif rexp
#'
#' @param n specified number of samples to be generated
#' @param init_theta initial theta value
#' @param TARGET function proportional to what we want to sample from
#' @param w TBA
#' @param max_steps TBA
#' @param k_product TBA
#' @param lower TBA
#' @param upper TBA
#' @param type either "raw" or "log" dependet upon if TARGET is log(f(x))
#' @param ... parameters passed to the TARGET function
#'
#' @return A vector of drawn values via slice sampler
#'
#'
#' @examples
#'
#' x <- 1
#'
#' @export
slice <- function(n, init_theta, TARGET,
                  w, max_steps, k_product,
                  lower = -Inf, upper = +Inf,
                  type = "raw", ...) {
  # u = theta = rep(NA, n)
  theta <- rep(NA, n)
  u <- matrix(NA, nrow = n, ncol = k_product) # k number of product slices
  theta[1] <- init_theta
  for (i in 2:n) {
    # u[i,] = runif(k_product, 0, TARGET(theta[i - 1]))
    # Find height of slice (depends on if TARGET is on raw or log scale)
    u[i, ] <- switch(type,
      raw = runif(k_product, min = 0, max = TARGET(theta[i - 1], ...)),
      log = TARGET(theta[i - 1], ...) - rexp(k_product)
    )

    # Sample with shrinkage and return
    theta[i] <- shrink_and_sample(
      theta = theta[i - 1],
      u = u[i, ],
      TARGET = TARGET,
      int = create_interval(
        theta = theta[i - 1],
        u = u[i, ],
        TARGET = TARGET,
        w = w,
        max_steps = max_steps,
        lower = lower,
        upper = upper
      ),
      ...
    )
  }
  # return(data.frame(theta = theta, u = u))
  return(theta)
}


# Find the slice interval [L,R] from which we will sample
# Arguments:
# theta - the current value of x
# u - current vertical level of slice
# TARGET - the function proportinal to the density
# w - estimate of typical slice size
# max_steps - integer limiting slice size to mw
# Returns:
# [L,R] - the slice interval
create_interval <- function(theta, u, TARGET, w, max_steps, lower, upper, ...) {
  L <- theta - runif(1, 0, w)
  R <- L + w
  # Step out
  # If setting J <- infinity, similar to MfUSampler::MfU.Sample, it requires slice.lower and slice.upper
  if (!is.infinite(max_steps)) {
    J <- floor(max_steps * runif(1))
    K <- (max_steps - 1) - J
    while (any(u < TARGET(L, ...)) && J > 0) { # should here be "all()" or "any()"
      L <- L - w # if any height smaller than the subfunction's value, stepping out; then ok for either density<1 or >1
      J <- J - 1
    }
    while (any(u < TARGET(R, ...)) && K > 0) { # should here be "all()" or "any()"
      R <- R + w
      K <- K - 1
    }
  } else {
    repeat{
      if (L <= lower) break
      if (all(u >= TARGET(L, ...))) break # if all endpoints outside the slice, then ok for either density<1 or >1
      L <- L - w
    }

    repeat{
      if (R >= upper) break
      if (all(u >= TARGET(R, ...))) break
      R <- R + w
    }
  }
  return(list(L = L, R = R))
}

# A function which now samples from the slice interval
# Arguments:
# theta - the current value of x
# u - current vertical level of slice
# TARGET - the function proportinal to the density
# int - the slice interval in form (L,R)
shrink_and_sample <- function(theta, u, TARGET, int, ...) {
  L <- int$L
  R <- int$R
  repeat {
    theta_prop <- runif(1, L, R)
    if (all(u < TARGET(theta_prop, ...))) {
      return(theta_prop)
    }
    # shrink
    if (theta_prop > theta) {
      R <- theta_prop
    }
    if (theta_prop < theta) {
      L <- theta_prop
    }
  }
}

# # target function from Example 8.3 page 327 of the book Robert & Casella (2nd ed, 2004)
# TARGET = function(x, ...) {
#   ##dnorm(x,-2)/2+dnorm(x,2)/2 # simple example
#
#   # multimodal density
#   (1 + sin(3*x)^2) * (1 + cos(5*x)^4) * exp(-x^2/2)  # if k_product=1
#
#   #c( (1 + sin(3*x)^2), (1 + cos(5*x)^4), exp(-x^2/2) ) # if k_product>1
# }
