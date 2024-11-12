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
             nIter = 50, burnin = 10, tick = 100)

#library(HI) 
#library(truncnorm)

set.seed(123)

# hyperparameters
tauA <- 20; tauB <- 50; tauSq <- 1 
#tauA <- 1.5; tauB <- 50; tauSq <- 1 
#wA <- 20; wB <- 50; wSq <- 1 
wA <- 5; wB <- 20; wSq <- 1 
vA <- 10; vB <- 20; v <- c(10, 1, 1)
#kappaA <- 3; kappaB <- 1 # This is for Gamma prior
#kappaA <- 3; kappaB <- 10#5#10 # This is for Inverse-Gamma prior
kappaA <- 5; kappaB <- 20

# initialization of parameters
kappas <- 0.9
## Survival mean part
betas.current <- matrix(0, nrow = NCOL(dat$x1), ncol = NCOL(dat$proportion))
mu.current <- matrix(0, nrow = NROW(dat$x1), ncol = NCOL(dat$proportion))
for(l in 1:L)
  mu.current[, l] <- exp(dat$beta0[l] + dat$XX[,,l] %*% betas.current[,l])
lambdas <- mu.current / gamma(1+1/kappas) # lambdas is a parameter in WEI3 distr.

## proportion Dirichlet part
phi <- 1
zetas.current <- matrix(0, nrow = NCOL(dat$x1)+1, ncol = NCOL(dat$proportion)-1) # include intercept
proportion <- matrix(0, nrow = NROW(dat$x1), ncol = NCOL(dat$proportion))
for(l in 1:(L-1)) # using the last cell type as reference
  proportion[, l] <- exp(cbind(1,dat$XX[, ,l]) %*% zetas.current[, l]) / 
    ( 1 + rowSums(sapply(1:(L-1), function(xx){exp(cbind(1,dat$XX[, ,xx]) %*% zetas.current[, xx])} )))
proportion[, L] <- 1 - rowSums(proportion[, -L])

betas <- mu.current / gamma(1+1/kappas)
weibull.S <- matrix(nrow = n, ncol = 3)
for(l in 1:L)
  weibull.S[, l] <- exp(-(dat$survObj$time / betas[, l] )^kappas)

xi <- c(0, 0, 0)
#dat$x0 <- cbind(1, dat$x01, dat$x02)
thetas <- exp(dat$x0 %*% xi)
# mode <- function(x) {
#   event <- density(x)
#   event$x[which.max(event$y)]
# }

M <- 50 # outside MCMC iterations
xi.mcmc <- matrix(0, nrow = 1+M, ncol = NCOL(dat$x0))
xi.mcmc[1, ] <- xi
kappas.mcmc <- c(kappas, rep(0, M))
phi.mcmc <- c(phi, rep(0, M))
tauSq.mcmc <- c(tauSq, rep(0, M))
#tauSq.mcmc <- matrix(0, nrow = 1+M, ncol = 3)
wSq.mcmc <- c(wSq, rep(0, M))
v.mcmc <- matrix(0, nrow = 1+M, ncol = 2)
betas.mcmc <- matrix(0, nrow = 1+M, ncol = NCOL(dat$proportion) * NCOL(dat$x1))
betas.mcmc[1, ] <- as.vector(betas.current)
zetas.mcmc <- matrix(0, nrow = 1+M, ncol = (NCOL(dat$proportion)-1) * (NCOL(dat$x1)+1))
zetas.mcmc[1, ] <- as.vector(zetas.current)
for(m in 1:M){
  if(m %% 100 == 0)
    cat("iteration:", m, "\n")
  
  ## update \xi's in cure fraction
  xi.mcmc.internal <- HI::arms(
    y.start = xi,
     myldens = logpost_xi2, 
     indFunc = convex_support, 
     n.sample = 5)
  ## n.sample = 20 will result in less variation
  xi <- colMeans(xi.mcmc.internal)#[-c(1:(nrow(xi.mcmc.internal)/2)),])
  xi <- sapply(xi, function(xx){min(abs(xx), 3-0.1) * sign(xx)} )
  #xi <- xi.mcmc.internal
  
  v <- sampleV(2, vA, vB, xi)
  v.mcmc[1+m, ] <- v[-1]
  
  thetas <- exp(dat$x0 %*% xi)
  xi.mcmc[1+m, ] <- xi
  
  ## update phi in Dirichlet measurement error model
  phi.mcmc.internal <- HI::arms(
    y.start = phi, 
    myldens = logpost_phi, 
    indFunc = convex_support_phi, 
    n.sample = 10)
  ## n.sample = 20 will result in less variation
  #phi <- mean(phi.mcmc.internal)#[-c(1:(length(phi.mcmc.internal)/2))])#median
  phi <- median(phi.mcmc.internal[-c(1:(length(phi.mcmc.internal)/2))])#median
  #phi <- phi.mcmc.internal
  phi.mcmc[1+m] <- median(c(phi, 0.1, 200-0.1))
  
  ## update \zeta_l of p_l in non-cure fraction; the last cell type as reference
  for(l in 1:(L-1)){
    for(j in 1:(p+1)) {
      if(j == 1) wSq <- 10
      zetas.mcmc.internal <- HI::arms(
        y.start = zetas.current[j, l],
        myldens = logpost_zeta_jl,  
        indFunc = convex_support, 
        n.sample = 20)
      ## n.sample = 20 will result in less variation
      #zeta.l.new <- mean(zetas.mcmc.internal)#[-c(1:(length(zetas.mcmc.internal)/2))]) #median
      zeta.l.new <- median(zetas.mcmc.internal[-c(1:(length(zetas.mcmc.internal)/2))]) #median
      #zeta.l.new <- zetas.mcmc.internal
      zetas.current[j, l] <- min(abs(zeta.l.new), 3-0.1) * sign(zeta.l.new)
      #wSq <- sampleW(1, wA, wB, zetas.current[j,l])
    } 
    for(ll in 1:(L-1)) # the l-th subtype proportion is updated, and the all composition should be updated
      proportion[, ll] <- exp(cbind(1,dat$XX[, ,ll]) %*% zetas.current[, ll]) / 
      ( 1 + rowSums(sapply(1:(L-1), function(xx){exp(cbind(1,dat$XX[, ,xx]) %*% zetas.current[, xx])} )))
    proportion[, L] <- 1 - rowSums(proportion[, -L])
  }
  zetas.mcmc[1+m, ] <- as.vector(zetas.current)
  
  ## update wSq, the variance of zetas
  wSq <- sampleW(1, wA, wB, zetas.current[-1,])
  wSq.mcmc[1+m] <- wSq
  
  ## update kappa in noncure fraction
  kappas.mcmc.internal <- HI::arms(
    y.start = kappas, 
    myldens = logpost_kappas, 
    indFunc = convex_support_kappas, 
    n.sample = 20)
  ## n.sample = 20 will result in less variation
  #kappas <- mean(kappas.mcmc.internal)#[-c(1:(length(kappas.mcmc.internal)/2))])#median
  kappas <- median(kappas.mcmc.internal[-c(1:(length(kappas.mcmc.internal)/2))])
  #kappas <- kappas.mcmc.internal
  kappas.mcmc[1+m] <- median(c(kappas, 0.1+0.1, 10-0.1))
  
  ## update \beta_l of S_l(t) in non-cure fraction
  for(l in 1:3){
    for(j in 1:10) {
      betas.mcmc.internal <- HI::arms(
        y.start = betas.current[j, l],
        myldens = logpost_beta_jl,  
        indFunc = convex_support, 
        #y.start = betas.current[, l],
        #myldens = logpost_beta_l,  
        n.sample = 20)
      ## n.sample = 20 will result in less variation
      #beta.l.new <- mean(betas.mcmc.internal)#[-c(1:(length(betas.mcmc.internal)/2))]) #median
      beta.l.new <- median(betas.mcmc.internal[-c(1:(length(betas.mcmc.internal)/2))]) #median
      #beta.l.new <- betas.mcmc.internal
      betas.current[j, l] <- min(abs(beta.l.new), 3-0.1) * sign(beta.l.new)
      #beta.l.new <- colMeans(betas.mcmc.internal[-c(1:(nrow(betas.mcmc.internal)/2)),])
      #betas.current[,l] <- sapply(beta.l.new, function(xx){min(abs(xx), 5-0.1) * sign(xx)} )
      #tauSq <- sampleTau(1, tauA, tauB, betas.current[j, l])
    } 
    mu.current[, l] <- exp(dat$XX[,,l] %*% betas.current[,l])
    lambdas <- mu.current / gamma(1+1/kappas) # lambdas is a parameter in WEI3 distr.
    weibull.S[, l] <- exp(-(dat$survObj$time / lambdas[, l] )^kappas)
    #tauSq[l] <- sampleTau(1, tauA, tauB, betas.current[, l])
  }
  betas.mcmc[1+m, ] <- as.vector(betas.current)
  
  ## update tauSq, the variance of betas
  #tauSq <- sampleTau(3, tauA, tauB, betas.current)
  #tauSq.mcmc[1+m, ] <- tauSq
  tauSq <- sampleTau(1, tauA, tauB, betas.current)
  tauSq.mcmc[1+m] <- tauSq
}
cat("... Done!\n")

max(xi.mcmc[-c(1:ceiling(nrow(xi.mcmc)/2)), ])
mean(xi.mcmc[, 1])
colMeans(xi.mcmc[-c(1:(M/2)),])
dat$xi
mean(kappas.mcmc[-c(1:(M/2))])
mean(phi.mcmc[-c(1:(M/2))])

# survival predictions based on posterior mean
xi.hat <- colMeans(xi.mcmc[-c(1:(M/2)), ])
beta.hat <- matrix(colMeans(betas.mcmc[-c(1:(nrow(betas.mcmc)/2)), ]), ncol=3)
kappa.hat <- mean(kappas.mcmc[-c(1:(M/2))])
thetas.hat <- exp(dat$x0 %*% xi.hat)

time_eval <- sort(dat$survObj$time)
Surv.prob <- matrix(nrow=n, ncol=length(time_eval))
for(j in 1:length(time_eval)) {
    tmp <- 0
    for(l in 1:3){
      mu <- exp(dat$XX[,,l] %*% beta.hat[,l])
      lambdas <- mu / gamma(1+1/kappa.hat) 
      weibull.S <- exp(-(time_eval[j] / lambdas )^kappa.hat)
      tmp <- tmp + dat$proportion[,l] * weibull.S
    }
    Surv.prob[,j] <- exp(-thetas.hat * (1-tmp))
}

pred.prob <- 1 - Surv.prob

library(survival)
library(riskRegression)
library(ggplot2)
# Cox model with clinical variables
x <- apply(dat$XX, c(1,2), mean); colnames(x) <- paste0("x", 1:10)
x.median <- apply(dat$XX, c(1,2), median); colnames(x.median) <- paste0("x.median", 1:10)
survObj <- data.frame(dat$survObj, x01=dat$x0[,2], x02=dat$x0[,3], x, x.median)
fitCox.clin <- survival::coxph(Surv(time, event) ~ x01 + x02, data = survObj, y=TRUE, x = TRUE)
#pred.clin <- survival::survfit(fitCox.clin, data = survObj)
formula.tmp <- as.formula(paste0("Surv(time, event) ~ ", paste0(colnames(x.median), collapse = "+")))
fitCox.X.median <- survival::coxph(formula.tmp, data = survObj, y=TRUE, x = TRUE)
#pred.X.median <- survival::survfit(fitCox.X.median, data = survObj)
formula.tmp2 <- as.formula(paste0("Surv(time, event) ~ ", paste0(paste0("x",1:10), collapse = "+")))
fitCox.X.mean <- survival::coxph(formula.tmp2, data = survObj, y=TRUE, x = TRUE)
#formula.tmp <- as.formula(paste0("Surv(time, event) ~ x01+x02+", paste0(colnames(x.median), collapse = "+")))
#fitCox.clin.X.median <- survival::coxph(formula.tmp, data = survObj, y=TRUE, x = TRUE)
formula.tmp <- as.formula(paste0("Surv(time, event) ~ x01+x02+", paste0(paste0("x",1:10), collapse = "+")))
fitCox.clin.X.mean <- survival::coxph(formula.tmp, data = survObj, y=TRUE, x = TRUE)

library(miCoPTCM) # good estimation for cure fraction; same BS as Cox.clin
resMY <- PTCMestimBF(Surv(time, event) ~ x01 + x02, data=survObj, varCov=matrix(0, nrow=3, ncol=3), init=rep(0, 3))
Surv.PTCM <- exp( - exp(dat$x0 %*% resMY$coefficients) %*% t(resMY$estimCDF) )
predPTCM.prob <- 1 - Surv.PTCM

a <- Score(list("Cox.clin"=fitCox.clin,
                "Cox.X.mean"=fitCox.X.mean,
                "Cox.X.median"=fitCox.X.median,
                #"Cox.clin.X.median"=fitCox.clin.X.median,
                "Cox.clin.X.mean"=fitCox.clin.X.mean,
                #"PTCM.clin"=predPTCM.prob,
                "GPTCM"=pred.prob),
           formula=Surv(time,event)~1,
           metrics="brier", summary = "ibs", 
           data=survObj,
           conf.int=F,times=time_eval)
a0 <- a$Brier$score
#a0 <- a0[a0$times <= 71, ]
levels(a0$model)[1] <- "Kaplan-Meier"
b <- ggplot(a0, aes(times, Brier, group = model, color = model)) +
  xlab("Time") +
  ylab("Brier score") +
  geom_step(direction = "vh") + #, alpha=0.4) +
  theme_bw() + 
  theme(legend.position = "inside", 
        legend.position.inside = c(0.4, 0.25), 
        legend.title = element_blank()) 
pdf(paste0("sim_brierSigma",Sigma,".pdf"), height = 4, width = 5)
b
dev.off()


# MCMC diagnosis
pdf(paste0("gptcm_xiHat_diagnosis_Sigma",Sigma,"_M",M,".pdf"), height = 5, width = 5)
ylabel <- paste0("expression(xi[", 0:2, "])")
layout(matrix(1:ncol(xi.mcmc), ncol = 1))
par(mar = c(2, 4.1, 2, 2))
for(j in 1:ncol(xi.mcmc)) {
  plot(xi.mcmc[, j], type = "l", lty = 1, ylab = eval(parse(text=ylabel[j])),
       ylim = summary(c(xi.mcmc, dat$xi))[c(1,6)])
  abline(h = dat$xi[j], col = "red")
}
dev.off()


pdf(paste0("gptcm_tauKappaPhi_diagnosis_Sigma",Sigma,"_M",M,".pdf"), height = 5, width = 5)
layout(matrix(1:4, ncol = 1))
par(mar = c(2, 4.1, 2, 2))
plot(kappas.mcmc, type = "l", lty = 1, 
     ylab = expression(kappa), xlab = "MCMC iteration",
     ylim = summary(c(kappas.mcmc, dat$kappas))[c(1,6)])
abline(h = dat$kappas, col = "red")
#for(l in 1:3)
plot(tauSq.mcmc, type = "l", lty = 1, 
     ylab = expression(tau^2), xlab = "MCMC iteration",
     ylim = summary(as.vector(tauSq.mcmc))[c(1,6)])
plot(wSq.mcmc, type = "l", lty = 1, 
     ylab = expression(w^2), xlab = "MCMC iteration",
     ylim = summary(as.vector(wSq.mcmc))[c(1,6)])
plot(phi.mcmc, type = "l", lty = 1, 
     ylab = expression(phi), xlab = "MCMC iteration",
     ylim = summary(as.vector(phi.mcmc))[c(1,6)])
dev.off()

# Final estimates
beta.mcmc0 <- betas.mcmc[-c(1:(nrow(betas.mcmc)/2)), ]
beta.true <- dat$betas[p:1, ]
beta_min <- min(beta.mcmc0, beta.true)
beta_max <- max(beta.mcmc0, beta.true)

pdf(paste0("gptcm_betaHat_Sigma",Sigma,"_M",M,".pdf"), height = 4, width = 8)
#pdf("gptcm_betaHat_noGamma.pdf", height = 4, width = 8)
layout(matrix(1:L, nrow=1))
plotCoeff(beta.mcmc0[, 1:p][,p:1], beta.true[, 1], xlim=c(beta_min, beta_max), main = "Cell type 1", label.y=paste0("x",p:1))
plotCoeff(beta.mcmc0[, p+1:p][,p:1], beta.true[, 2], xlim=c(beta_min, beta_max), main = "Cell type 2", label.y=paste0("x",p:1))
plotCoeff(beta.mcmc0[, 2*p+1:p][,p:1], beta.true[, 3], xlim=c(beta_min, beta_max), main = "Cell type 3", label.y=paste0("x",p:1))
dev.off()

zeta.mcmc0 <- zetas.mcmc[-c(1:(nrow(zetas.mcmc)/2)), ]
zeta.true <- dat$zetas[(p+1):1, ]
zeta_min <- min(zeta.mcmc0, zeta.true)
zeta_max <- max(zeta.mcmc0, zeta.true)
pdf(paste0("gptcm_zetaHat_Sigma",Sigma,"_M",M,".pdf"), height = 5, width = 7)
layout(matrix(1:(L-1), nrow=1))
plotCoeff(zeta.mcmc0[, 1:(p+1)][,(p+1):1], zeta.true[, 1], xlim=c(zeta_min, zeta_max), main = "Cell type 1", label.y=c(paste0("x",p:1), "intecept"))
plotCoeff(zeta.mcmc0[, 1+p+1:(p+1)][,(p+1):1], zeta.true[, 2], xlim=c(zeta_min, zeta_max), main = "Cell type 2", label.y=c(paste0("x",p:1), "intecept"))
dev.off()

# MCMC traceplots of betas
pdf(paste0("gptcm_betaHat_diagnosis_Sigma",Sigma,"_M",M,".pdf"), height = 18, width = 8)
ylabel <- paste0("expression(beta['", rep(1:p,3),",", rep(1:3,each=p), "'])")
layout(matrix(1:ncol(betas.mcmc), ncol = 3))
par(mar = c(2, 4.1, 2, 2))
for(j in 1:ncol(betas.mcmc)) {
  plot(betas.mcmc[, j], type = "l", lty = 1, ylab = eval(parse(text=ylabel[j])),
       ylim = summary(c(betas.mcmc, dat$betas))[c(1,6)])
  abline(h = dat$betas[j], col = "red")
}
dev.off()

# MCMC traceplots of zetas
pdf(paste0("gptcm_zetaHat_diagnosis_Sigma",Sigma,"_M",M,".pdf"), height = 18, width = 12)
ylabel <- paste0("expression(zeta['", rep(0:p,L-1),",", rep(1:(L-1),each=p+1), "'])")
layout(matrix(1:ncol(zetas.mcmc), nrow = p+1))
par(mar = c(2, 4.1, 2, 2))
for(j in 1:ncol(zetas.mcmc)) {
  plot(zetas.mcmc[, j], type = "l", lty = 1, ylab = eval(parse(text=ylabel[j])),
       ylim = summary(c(zetas.mcmc, dat$zetas))[c(1,6)])
  abline(h = dat$zetas[j], col = "red")
}
dev.off()
