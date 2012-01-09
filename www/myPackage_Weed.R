
###############################
##### Example for Weed Data
###############################

rm(list = ls(all = TRUE))

#################
# Load myPackage
require(myPackage)
myPackage_hello_world()
data(package = "myPackage")

# load data and plot
data(datWeed)
dat <- Weed
n <- nrow(dat)
YRaw <- dat[,3]
locRaw <- dat[,1:2]
locRaw.u <- unifLoc(locRaw)
plotData(YRaw, locRaw, xlab="", ylab="")

# randomly choose 60 observations as training data
set.seed(111)
ind <- sample(1:n, 60, replace=FALSE)
L <- rep(1, length(ind))
Y <- YRaw[ind]
loc <- locRaw.u[ind,]
plotData(Y/L, loc, size=c(0.3, 3.7), xlab="", ylab="")

# tuning and run MCMC algorithm
input0 <- MCMCinput( run=2000, run.S=1, rho.family="rhoPowerExp", 
          Y.family = "Poisson", ifkappa=0,
          scales=c(0.5, 3.5, 0.9, 0.6, 0.5), 
          phi.bound=c(0.005, 1), 
          initials=list(c(-1), 1, 0.1, 1) )

# Single chain MCMC
res <- runMCMC(Y, L = L, loc=loc, X=NULL, MCMCinput = input0 )
# Cut chains
res.m <- cutChain(res, chain.ind=1:4, burnin=500, thining=1)

##################################### NOT RUN
# Mltiple chain MCMC 
# !!!CAUTION!!!
# runMCMC.multiChain() won't works on R with optimized BLAS
# runMCMC.sf works but won't' gain significant speed

# parallel with {multicore}
require(multicore)
res.prl <- runMCMC.multiChain(Y, L = L, loc=loc, X=NULL, 
            MCMCinput = input0, n.chn = 2, n.cores = 2)
# parallel with {snowfall}
require(snowfall)
res.prl <- runMCMC.sf(Y, L = L, loc=loc, X=NULL, 
            MCMCinput = input0, n.chn = 8, n.cores = 8)
# Cut chains                    
res.m.prl <- lapply(res.prl, cutChain, chain.ind=1:4, burnin=200, thining=10)
res.mix <- mixChain(res.m.prl)
res.m <- res.mix                    
##################################### NOT RUN

# Examine chains
require(coda)
## for covariance matrix parameters
chn1 <- cbind(sigma=res.m$s, phi=res.m$a)
chn1.mcmc <- mcmc(chn1); dim(chn1.mcmc)
summary(chn1.mcmc)
plot(chn1.mcmc, auto.layout = TRUE)
crosscorr.plot(chn1.mcmc)
autocorr.plot(chn1.mcmc)
effectiveSize(chn1.mcmc)
geweke.diag(chn1.mcmc, frac1=0.1, frac2=0.5)
heidel.diag(chn1.mcmc, eps=0.1, pvalue=0.05)
## for coefficients 
if(is.matrix(res.m$m)){
  beta.mcmc <- mcmc(t(res.m$m))
  } else beta.mcmc <- mcmc(res.m$m)
plot(beta.mcmc, auto.layout = TRUE)
crosscorr.plot(beta.mcmc)
autocorr.plot(beta.mcmc)
## for latent variables 
S.mcmc <- mcmc(t(res.m$S))
dim(S.mcmc)
plotACF(S.mcmc)
crosscorr.plot(S.mcmc)

# select the remaining locations to be predicted
ind.left <- setdiff(1:n, ind)
ind.p <- ind.left
locp <- locRaw.u[ind.p,] 
Lp <- rep(1, nrow(locp)) 
# prediction
Ypred <- predY(res.m, loc, locp, X=NULL, Xp=NULL, Lp=Lp, k=1, 
               rho.family="rhoPowerExp", Y.family="Poisson"
               #, mc = TRUE, n.cores = 4
               )
Ypred.avg <- rowMeans(Ypred$Y); 
Ypred.mode <- apply(Ypred$Y, 1, findMode)
# plot observed, predicted and actual data
#op <- par(mfrow=c(2,1))
plotData(Y/L, locRaw[ind,], Ypred.avg/Lp, locRaw[ind.p,], 
          YRaw[ind.p]/Lp, locRaw[ind.p,], xlab="", ylab="")
plotData(Y/L, locRaw[ind,], Ypred.mode/Lp, locRaw[ind.p,], 
          YRaw[ind.p]/Lp, locRaw[ind.p,], xlab="", ylab="")
plot(Ypred.avg/Lp, YRaw[ind.p]/Lp, xlab="Predicted values", ylab="Actual values")
points(Ypred.mode/Lp, YRaw[ind.p]/Lp, pch=2, col=3)
abline(c(0,1), col=2)
#par(op)
