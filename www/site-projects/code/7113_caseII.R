##########################
### Porject 1 for STA 7113
### Case 2
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


sigma.u <- NULL
Au <- NULL
Bu <- NULL
l <- 3
q <- rep(1,l)
for(i in 1:l){
	sigma.u[i] <- round(runif(1,0.1,1),1)
	Au[i] <- 0.1
	Bu[i] <- 0.1
}
Ae <- 0.01
Be <- 0.01
sigma.e <- 0.5
n <- nrow(X)
res <- NULL
iter <- 2000

for(i in 1:iter){

G <- diag(sigma.u)
B <- diag(c(rep(0,p+1),diag(solve(G))))
mu <- solve(t(C)%*%C+sigma.e*B)%*%t(C)%*%y
Va <- sigma.e*solve(t(C)%*%C+sigma.e*B)
b_u <- rmvnorm(1,mu,Va)
b <- b_u[1:(p+1)]
u <- b_u[-(1:(p+1))]

end <- 0
for(i in 1:l){
sigma.u[i] <- 1/rgamma(1,Au[i]+0.5*q[i],Bu[i]+0.5*sum(u[(end+1):(end+q[i])]^2))
end <- end+q[i]
}

sigma.e <- 1/rgamma(1,Ae+0.5*n,Be+0.5*sum((y-X%*%b-Z%*%u)^2))
res <- rbind(res,c(b,u,sigma.u,sigma.e))
}

########### Result Analysis
###########################
dim(res)
colnames(res) <- c("beta[1]",  "beta[2]",  "beta[3]", 
	 "beta[4]",  "b[1]",  "b[2]", "b[3]", 
	"sigma.b1",  "sigma.b2", "sigma.b3", "sigma.e")
apply(res,2,median)
matplot(res[,-(8:10)],type='l',ann=F)
title(xlab="iteration",ylab="results",main="Case 2")
matplot(res[,c(8:10)],type='l',ylim=c(0,500),ann=F)
title(xlab="iteration",ylab="results",main="Case 2")

#res2 <- res
#rm(list=ls(all=T)[ls(all=T)!="res2"])



