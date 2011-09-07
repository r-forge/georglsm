##########################
### Porject 1 for STA 7113
### Case 1
##########################

rm(list=ls(all=T))
#########################

########### Data Initializing
#############################
#n2 <- read.table("D:/Study/STA 7113/project/no2.txt")
n2 <- read.table("X:/liang.jing/STA 7113/GLMM/no2.txt")
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
res <- NULL
iter <- 2000
D <- diag(c(rep(0,p+1),rep(1,K)))

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

########### Result Analysis
###########################
apply(res,2,median)
#matplot(res,type='l')
#matplot(res[,7],type='l')
matplot(res[,-8],type='l',ann=F)
title(xlab="iteration",ylab="results",main="Case 1")
matplot(res[,8],type='l',col=2,ann=F)
title(xlab="iteration",ylab="results",main="Case 1")

for(i in 1:9){
plot(density(res[,i]),col=i,xlim=c(-0.5,1),ylim=c(0,16),ann=F)
par(new=T)
}
par(new=F)
title(xlab="x",ylab="Density",main="Posterior Density")

s.b <- res[,(p+2+K)]
s.e <- res[,(p+3+K)]
windows()
plot(s.b,type="l",ann=F)
points(s.e,type="l",col=3)
sb.c <- s.b[s.b<quantile(s.b,0.95)]
points(sb.c,type="l",col=2)

########### using WINBUGS
#########################
library(R2WinBUGS)

data <- list(n=n,p=p,K=K,y=y,X=X,Z=Z)
para <- c("beta","b","taue","taub")
inits <- function(){list(beta=rep(0.1,p+1),b=rep(0.1,K),taub=0.001,taue=0.001)}
#no2.sim <- bugs(data,inits,para,"X:/liang.jing/STA 7113/GLMM/no2.bug",n.iter=2000,n.chains=1,
#	bugs.directory="c:/Program Files (x86)/WinBugs/")
no2.sim <- bugs(data,inits,para,"X:/liang.jing/STA 7113/GLMM/no2.bug",n.iter=2000,n.chains=1,
	bugs.directory="C:/Program Files/WinBUGS14")
#print(no2.sim)
res.sim <- no2.sim$sims.array
dim(res.sim)
res.sim[,,(p+2+K):(p+3+K)] <- 1/res.sim[,,(p+2+K):(p+3+K)]
dimnames(res.sim)[[3]][(p+2+K):(p+3+K)] <- c("sigma.e","sigma.b")
apply(res.sim,c(2,3),median)

b.sim <- res.sim[,,1:(p+1)]
u.sim <- res.sim[,,(p+2):(p+1+K)]
sigma.e.sim <- res.sim[,,(p+2+K)]
sigma.b.sim <- res.sim[,,(p+3+K)]
windows()
matplot(cbind(b.sim,u.sim),type="l")
windows()
plot(sigma.b.sim,type="l",ann=F)
points(sigma.e.sim,type="l",col=3)
sb.sim <- sigma.b.sim[sigma.b.sim<quantile(sigma.b.sim,0.95)]
points(sb.sim,type="l",col=2)

########### Result Comparison
#############################
windows()
plot(1:2000,rep(0,2000),type="n",ylim=c(0,3),ann=F)
points(1:1000,s.b[1:1000],type="l",col=1)
points(1001:2000,sigma.b.sim[1:1000],type="l",col=1)
points(1:1000,s.e[1:1000],type="l",col=3)
points(1001:2000,sigma.e.sim[1:1000],type="l",col=2)
points(1:1000,sb.c[1:1000],type="l",col=2)
points(1001:2000,sb.sim[1:1000],type="l",col=3)
abline(v=1000,ylim=c(0,3),col=4,lwd=3,lty=4)
title(xlab="",ylab="",main="Result Comprison")
legend(300,3,c("R result:"),text.col=1,bty="n")
legend(300,2.85,c(expression(sigma[b]),expression(sigma[b]^{"w/o outliers"}),expression(sigma[e])),bty="n",
	lty=rep(1,3),col=1:3,text.col=1)
legend(1100,3,c("WinBUGS result:"),text.col=1,bty="n")
legend(1100,2.85,c(expression(sigma[b]),expression(sigma[b]^{"w/o outliers"}),expression(sigma[e])),bty="n",
	lty=rep(1,3),col=c(1,3,2),text.col=1)


res2 <- res.sim[,1,]
res2 <- cbind(res2[,1:(p+1+K)], sigma.b=res2[,(p+3+K)], sigma.e=res2[,(p+2+K)])

op <- par(mfrow=c(3,3))
for(i in 1:ncol(res)){
plot(density(res[,i]),col=i,xlim=range(c(res[,i],res2[,i]))	
	,ylim=range(c(density(res[,i])$y,density(res2[,i])$y)),ann=F)
par(new=T)
plot(density(res2[,i]),col=i+3,xlim=range(c(res[,i],res2[,i]))	
	,ylim=range(c(density(res[,i])$y,density(res2[,i])$y)),ann=F)
title(xlab=paste("",colnames(res)[i]),ylab="density")
}
par(op)




