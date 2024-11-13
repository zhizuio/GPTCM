# simulate data
set.seed(123)
n <- 200 # subjects
p <- 10 # variable selection predictors
L <- 3 # cell types
library(GPTCM)
dat <- simData(n, p, L)

## run Bayesian GPTCM
set.seed(123)
fit <- GPTCM(dat, n, p, L, nIter = 50, burnin = 10)

test_that("fit has properly class and length", {
  expect_s3_class(fit, "GPTCM")
  expect_length(fit, 3L)
  expect_length(fit$input, 3L)
  expect_length(fit$input$hyperpar, 8L)
  expect_length(fit$output, 2L)
  expect_length(fit$output$posterior, 5L)
  expect_length(fit$output$mcmc, 8L)
})
