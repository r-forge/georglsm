##########################
### Case 1 Comparison
##########################

rm(list=ls(all=T))
#########################

########### Data Initializing
#############################
n2 <- read.table("no2.txt")
no2.dat <- data.frame(no2=n2$V1,cars=n2$V2,temp=n2$V3,wind=n2$V4,
tempdiff=n2$V5,wdirection=n2$V6,hour=n2$V7,day=n2$V8)
#pairs(no2.dat)
no2.dat <- scale(no2.dat,center=T,scale=T)

y <- no2.dat[,c(1)]
X <- cbind(rep(1,nrow(no2.dat)),no2.dat[,c(2,3,4)]) #n*(p+1) matrix
Z <- no2.dat[,c(5,6,7)] #n*K matrix
X <- as.matrix(X)
Z <- as.matrix(Z)
C <- as.matrix(cbind(X,Z))
p <- ncol(X)-1
K <- ncol(Z)

########### MCMC iteration
##########################
library(mvtnorm)

n <- nrow(X)
Au <- 0.001
Bu <- 0.001
Ae <- 0.001
Be <- 0.001

sigma.e <- 0.5	
sigma.u <- 0.5
iter <- 2000
D <- diag(c(rep(0,p+1),rep(1,K)))

### "sample"
############
res <- NULL
for(i in 1:iter){
mu <- solve(t(C)%*%C+sigma.e/sigma.u*D)%*%t(C)%*%y
Va <- sigma.e*solve(t(C)%*%C+sigma.e/sigma.u*D)
b_u <- rmvnorm(1,mu,Va)
b <- b_u[1:(p+1)]
u <- b_u[-(1:(p+1))]
sigma.u <- 1/rgamma(1,Au+0.5*K,Bu+0.5*sum(u^2))
sigma.e <- 1/rgamma(1,Ae+0.5*n,Be+0.5*sum((y-X%*%b-Z%*%u)^2))
res <- rbind(res,c(b,u,sigma.u,sigma.e))
}
dim(res)
colnames(res) <- c("beta[1]",  "beta[2]",  "beta[3]", 
	 "beta[4]",  "b[1]",  "b[2]", "b[3]", "sigma.b",  "sigma.e")

### "maximize"
##############
res1 <- NULL
library(MCMCpack)
for(i in 1:iter){
mu <- solve(t(C)%*%C+sigma.e/sigma.u*D)%*%t(C)%*%y
Va <- sigma.e*solve(t(C)%*%C+sigma.e/sigma.u*D)
b_u <- rmvnorm(1,mu,Va)
b <- b_u[1:(p+1)]
u <- b_u[-(1:(p+1))]
f <- function(x=0.1) -sum(dinvgamma(x,Au+0.5*K,Bu+0.5*sum(u^2)))
#f <- function(x=1) -sum(dnorm(x,1.5,0.05))
sigma.u <- coef(mle(f))
sigma.e <- 1/rgamma(1,Ae+0.5*n,Be+0.5*sum((y-X%*%b-Z%*%u)^2))
res1 <- rbind(res1,c(b,u,sigma.u,sigma.e))
}
dim(res1)
colnames(res1) <- c("beta[1]",  "beta[2]",  "beta[3]", 
	 "beta[4]",  "b[1]",  "b[2]", "b[3]", "sigma.b",  "sigma.e")

############# Compare
#####################
apply(res,2,summary)
apply(res1,2,summary)
matplot(res,type='l')
matplot(res1,type='l')

op <- par(mfrow=c(3,3))
for(i in 1:ncol(res)) matplot(cbind(res[,i],res1[,i]),type='l')
par(op)

op <- par(mfrow=c(3,3))
for(i in 1:ncol(res)){
plot(density(res[,i]),col=i,xlim=range(c(res[,i],res1[,i]))	
	,ylim=range(c(density(res[,i])$y,density(res1[,i])$y)),ann=F)
par(new=T)
plot(density(res1[,i]),col=i+3,xlim=range(c(res[,i],res1[,i]))	
	,ylim=range(c(density(res[,i])$y,density(res1[,i])$y)),ann=F)
title(xlab=paste("",colnames(res)[i]),ylab="density")
}
par(op)

i <- 8
plot(density(res[,i]),col=i,xlim=c(0,.1)	
	,ylim=range(c(density(res[,i])$y,density(res1[,i])$y)),ann=F)
par(new=T)
plot(density(res1[,i]),col=i+3,xlim=c(0,.1)	
	,ylim=range(c(density(res[,i])$y,density(res1[,i])$y)),ann=F)
title(xlab=paste("",colnames(res)[i]),ylab="density")


