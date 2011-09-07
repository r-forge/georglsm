##########################
### Case 3
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

y <- sample(c(0,1),500,replace=T)


########### MCMC iteration
##########################
library(mvtnorm)
library(msm)

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
n <- nrow(X)
res <- NULL
iter <- 2000
b <- runif(p+1,0.1,0.5)
u <- runif(K,0.1,0.5)
a <- t(rmvnorm(1,X%*%b+Z%*%u,diag(1,length(y))))

for(i in 1:iter){

G <- diag(sigma.u)
B <- diag(c(rep(0,p+1),diag(solve(G))))
mu <- solve(t(C)%*%C+B)%*%t(C)%*%a
Va <- solve(t(C)%*%C+B)
b_u <- rmvnorm(1,mu,Va)
b <- b_u[1:(p+1)]
u <- b_u[-(1:(p+1))]

end <- 0
for(i in 1:l){
sigma.u[i] <- 1/rgamma(1,Au[i]+0.5*q[i],Bu[i]+0.5*sum(u[(end+1):(end+q[i])]^2))
end <- end+q[i]
}

index.0 <- c(1:length(y))[y==0]
index.1 <- c(1:length(y))[y==1]
#sort(unique(c(index.1,index.0)))
a[index.0] <- rtnorm(length(index.0),c(X%*%b+Z%*%u)[index.0],1,lower=-Inf,upper=0)
a[index.1] <- rtnorm(length(index.1),c(X%*%b+Z%*%u)[index.1],1,lower=0,upper=Inf)

res <- rbind(res,c(b,u,sigma.u))
}

########### Result Analysis
###########################
dim(res)
colnames(res) <- c("beta[1]",  "beta[2]",  "beta[3]", 
	 "beta[4]",  "b[1]",  "b[2]", "b[3]", 
	"sigma.b1",  "sigma.b2", "sigma.b3")
apply(res,2,median)
matplot(res,type='l', ylim=c(-0.5,100))
matplot(res[,-(8:10)],type='l',ann=F)
title(xlab="iteration",ylab="results",main="Case 3")
matplot(res[,c(8:10)],type='l',ylim=c(0,500),ann=F)
title(xlab="iteration",ylab="results",main="Case 3")




