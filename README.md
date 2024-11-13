# GPTCM
<!-- 
[![CRAN
status](https://www.r-pkg.org/badges/version/GPTCM)](https://cran.r-project.org/package=GPTCM)
[![r-universe](https://ocbe-uio.r-universe.dev/badges/GPTCM)](https://ocbe-uio.r-universe.dev/GPTCM)
[![License](https://img.shields.io/badge/License-GPLv3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/doi-10.32614%2FCRAN.package.GPTCM-brightgreen)](https://doi.org/10.32614/CRAN.package.GPTCM)
-->
[![R-CMD-check](https://github.com/zhizuio/GPTCM/workflows/R-CMD-check/badge.svg)](https://github.com/zhizuio/GPTCM/actions)


This is an R package **GPTCM** implementing Bayesian hierarchical modeling for a generalized promotion time cure model (GPTCM) for the identification of cell-type-specific tumor driver genes and survival prognosis. 

## Installation

Install the latest development version from [GitHub](https://github.com/zhizuio/GPTCM)

```r
#install.packages("remotes")
remotes::install_github("zhizuio/GPTCM")
```

## Examples

```{r}
rm(list=ls())

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

# draw time-dependent Brier scores
plotBrier(dat, datMCMC = fit)

# MCMC diagnosis
plotMCMC(dat, datMCMC = fit, estimator = "xi")
plotMCMC(dat, datMCMC = fit, estimator = c("kappa", "tau", "w", "phi"))

# Final estimates
plotCoeff(dat, datMCMC = fit, estimator = "beta")
plotCoeff(dat, datMCMC = fit, estimator = "zeta")

# MCMC traceplots of betas
plotMCMC(dat, datMCMC = fit, estimator = "beta")

# MCMC traceplots of zetas
plotMCMC(dat, datMCMC = fit, estimator = "zeta")

```


## References

Chen, M.-H., Ibrahim, J. G., and Sinha, D. (1999). A new Bayesian model for survival data with
a surviving fraction. Journal of the American Statistical Association, 94(447):909--919.

Yakovlev, A. (1996). Threshold models of tumor recurrence. Mathematical and Computer Modelling,
23(6):153--164.
