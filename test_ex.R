# === === === === === === === === === === === === === === === === 
# Test ARMS sampler for our new promotion time cure model ####
# Use 'miCoPTCM' R-pkg for comparison
# Use 'HI' R-pkg to update each beta separately
# Simulation code was adapted from the example in 'miCoPTCM::PTCMestimBF()'
# === === === === === === === === === === === === === === === === 

rm(list=ls())

set.seed(123)
n <- 200 # subjects
p <- 10 # variable selection predictors
L <- 3 # cell types
#dat <- simSurvObj(n, p, L)
Sigma <- 0 # "Sigma=0" means with dependent X; "Sigma=1" means with independent X
library(GPTCM)
dat <- simData(n, p, L, Sigma = diag(p*L) * Sigma)

set.seed(123)
fit <- GPTCM(dat, n, p, L,
             nIter = 50, burnin = 10)

# survival predictions based on posterior mean
library(survival)
b <- plotBrier(dat, datMCMC = fit)
#pdf(paste0("sim_brierSigma",Sigma,".pdf"), height = 4, width = 5)
b
#dev.off()

# MCMC diagnosis
#pdf(paste0("gptcm_xiHat_diagnosis_Sigma",Sigma,"_M",M,".pdf"), height = 5, width = 5)
plotMCMC(dat, datMCMC = fit, estimator = "xi")
#dev.off()

#pdf(paste0("gptcm_tauKappaPhi_diagnosis_Sigma",Sigma,"_M",M,".pdf"), height = 5, width = 5)
plotMCMC(dat, datMCMC = fit, estimator = c("kappa", "tau", "w", "phi"))
#dev.off()

# Final estimates

#pdf(paste0("gptcm_betaHat_Sigma",Sigma,"_M",M,".pdf"), height = 4, width = 8)
##pdf("gptcm_betaHat_noGamma.pdf", height = 4, width = 8)
plotCoeff(dat, datMCMC = fit, estimator = "beta")
#dev.off()

#pdf(paste0("gptcm_zetaHat_Sigma",Sigma,"_M",M,".pdf"), height = 5, width = 7)
plotCoeff(dat, datMCMC = fit, estimator = "zeta")
#dev.off()

# MCMC traceplots of betas
#pdf(paste0("gptcm_betaHat_diagnosis_Sigma",Sigma,"_M",M,".pdf"), height = 18, width = 8)
plotMCMC(dat, datMCMC = fit, estimator = "beta")
#dev.off()

# MCMC traceplots of zetas
#pdf(paste0("gptcm_zetaHat_diagnosis_Sigma",Sigma,"_M",M,".pdf"), height = 18, width = 12)
plotMCMC(dat, datMCMC = fit, estimator = "zeta")
#dev.off()
