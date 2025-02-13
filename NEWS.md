# GPTCM 0.0.8

* Test C/C++ (derivative-free) ARMS from CRAN pkg 'genscore' that is close to Gilks' original C code

# GPTCM 0.0.7

* Add C++ ARS (not derivative-free) adapted from GitHub 'kuperov/adsample' for generating $\xi$, but it is strange that it fails when the lower bound is negative. 

# GPTCM 0.0.6

* Add Bayesian variable selection (BVS), not yet finished
* When using `HI::arms()`, setting `n.sample=2`, but somehow not work, not know why

# GPTCM 0.0.5

* Add Bayesian variable selection (BVS), not yet finished
* Using log(concentrations) for modeling proportions by default
* In function `GPTCM()`, transforming proportions data with values (very close to) 0 or 1 in the same way as DirichletReg::DR_data()
* Using log-scale for every fraction of the Dirichlet density to avoid numeric issues

# GPTCM 0.0.4

* Bring back R-package `HI` 
* Add function of slice sampler

# GPTCM 0.0.3

* Remove R-package `HI` and integrate the R and C functions from `HI`

# GPTCM 0.0.2

* Add modeling of the proportions through the common parametrization of Dirichlet regression, i.e. without reference category
* Add options of IGamma and Gamma prior for Weibull's shape parameter kappa
* Change the prior variance of zeta's intercept from fixed value 10 to an IGamma prior

# GPTCM 0.0.1

* First beta version
