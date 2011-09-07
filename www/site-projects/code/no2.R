
#########################################
rm(list=ls(all=T))
#########################################


#########################################
#### Import and initialize data set
#########################################
# setwd("...")
n2 <- read.table("no2.txt")
no2.dat <- data.frame(no2=n2$V1,cars=n2$V2,temp=n2$V3,wind=n2$V4,
tempdiff=n2$V5,wdirection=n2$V6,hour=n2$V7,day=n2$V8)

no2.cor <- cor(no2.dat)    #examine the correlations
diag(no2.cor) <- 0
range(no2.cor)

#tt <- no2.dat[,-6]
#write.table(tt, file = "tt.txt", sep = "   ",row.names = F, col.names = F) 

Y <- no2.dat[,1]
hist(Y,br=33,prob=T)
curve(dnorm(x,mean(Y),sd(Y)),add=T,col=2,lwd=2)
curve(dgamma(x,shape=mean(Y)^2/var(Y),scale=var(Y)/mean(Y)),add=T,lty=2,col=3,lwd=2)
legend(x="topleft", c("normal","gamma"), col=2:3, lty=1:2, ncol = 2, cex = 0.8)

pairs(no2.dat,cex=0.6)

cor(no2.dat[,-1])

#########################################
##### Bayesian Model Averaging Functions
#########################################
library(rrcov)
library(BMA)

#Function to identify potential outliers
outliers <- out.ltsreg(no2.dat[,-1],no2.dat[,1], delta=2)
length(outliers)

#BMA for linear regression models
no2.reg <- bicreg(no2.dat[,-1],no2.dat[,1],strict = FALSE, OR = 20)
su.reg <- summary(no2.reg)
plot(no2.reg,mfrow=c(4,2),include.intercept=T)
imageplot.bma(no2.reg)

#BMA for generalized linear models
no2.glm <- bic.glm(no2.dat[,-1],no2.dat[,1],strict = FALSE, OR = 20,glm.family="Gamma")
summary(no2.glm)
plot(no2.glm,mfrow=c(4,2),include.intercept=T)
imageplot.bma(no2.glm)

#Bayesian simultaneous variable selection and outlier identification
no2.mc3 <-MC3.REG(no2.dat[,1],no2.dat[,-1],num.its=1000,outliers = FALSE)
summary(no2.mc3)
no2.mc3[1:24,]

#Iterated Bayesian Model Averaging variable selection for GLM, linear models
no2.iBMA.bicreg <- iBMA.bicreg(no2.dat[,-1],no2.dat[,1], thresProbne0 = 5,
verbose = TRUE, maxNvar = 30, nIter =2000)
summary(no2.iBMA.bicreg)
orderplot(no2.iBMA.bicreg)

no2.iBMA.glm <- iBMA.glm(no2.dat[,-1],no2.dat[,1], glm.family="gaussian",nIter = 1000)
summary(no2.iBMA.glm)
orderplot(no2.iBMA.glm)

#Model uncertainty in generalized linear models using Bayes factors 
############ (glib function NOT WORKING YET!!! it is not necessary to use though
no2.glib <- glib(no2.dat[,-1],no2.dat[,1],error="gaussian", link="identity",
glimvar=TRUE, output.priorvar=TRUE, output.postvar=TRUE)
summary(no2.glib)



#########################################
#### Frequentist Linear model
#########################################
no2.1 <- lm(no2~.,no2.dat)
Y.flm <- predict(no2.1)
Y <- no2.dat[,1]
e.flm <- mean(Y-Y.flm)
e1.flm <- mean(abs(Y-Y.flm))
s.flm <- sum((Y-Y.flm)^2)/(length(Y)-1)

#res <- residuals(no2.1)
#plot(Y.flm,res)
#su.FLM <- summary(no2.1)
#qqnorm(res)
#qqline(res,col=2)

#########################################
#### Frequentist Generalized Linear model
#########################################
no2.fglm.n <- glm(no2~.,no2.dat,family=gaussian)
Y.fglm.n <- predict(no2.fglm.n,no2.dat)
e.fglm.n <- mean(Y-Y.fglm.n)
e1.fglm.n <- mean(abs(Y-Y.fglm.n))
s.fglm.n <- sum((Y-Y.fglm.n)^2)/(length(Y)-1)

no2.fglm.ga <- glm(no2~.,no2.dat,family=Gamma)
Y.fglm.ga <- 1/predict(no2.fglm.ga,no2.dat)
e.fglm.ga <- mean(Y-Y.fglm.ga)
e1.fglm.ga <- mean(abs(Y-Y.fglm.ga))
s.fglm.ga <- sum((Y-Y.fglm.ga)^2)/(length(Y)-1)

no2.fglm.ig <- glm(no2~.,no2.dat,family=inverse.gaussian)
Y.fglm.ig <- (predict(no2.fglm.ig,no2.dat))^(-1/2)
e.fglm.ig <- mean(Y-Y.fglm.ig)
e1.fglm.ig <- mean(abs(Y-Y.fglm.ig))
s.fglm.ig <- sum((Y-Y.fglm.ig)^2)/(length(Y)-1)

#######################################################
#### Prediction Power Comparison (for full data set)(bicreg)
#######################################################
no2.reg <- bicreg(no2.dat[,-1],no2.dat[,1],strict = FALSE, OR = 20)
model.reg <- no2.reg[[7]] #top models
P <- no2.reg[[1]]         #model posterio probability
Beta <- no2.reg[[14]]     #model beta_hat
SE <- no2.reg[[15]]       #model sd
summary(no2.reg)

X <- as.matrix(no2.dat[,-1])
X <- cbind(rep(1,nrow(X)),X)
Y <- no2.dat[,1]

Y.reg <- X%*%t(Beta)      #predict value for all models
P.reg <- matrix(rep(P,nrow(Y.reg)),nrow(Y.reg),,byrow=T)
Y.p <- Y.reg*P.reg
Y.a <- apply(Y.p,1,sum)   #average predict for all models
Y.reg.a <- Y.a

if(F){
###trash
library(lattice)
n <- length(Y.reg.a)
temp <- data.frame(Y=c(Y.reg.a,Y.reg[,1],Y.o),id=c(rep(1,n),rep(2,n),rep(3,n)))
bwplot(Y|id,group=id,temp)
###trash end
}

y.lim <- c(min(Y.a,Y.reg[,1]),max(Y.a,Y.reg[,1]))
plot(Y,Y.a,col=3,pch=2,ylab="Y_hat",ylim=y.lim,main="BMA LM Averaged vs Best")
par(new=T)
plot(Y,Y.reg[,1],col=4,pch=3,ylim=y.lim,ann=F)
abline(c(0,1),col=2,lwd=2)
legend(x="topleft", c("Averaged","Best"), col=c(3,4), pch=c(2,3), ncol = 2, cex = 0.8)

resid.a <- Y-Y.a
resid.1 <- Y-Y.reg[,1]
plot(resid.a,col=3,ylim=c(-2.5,2),pch=2,ylab="Residual",main="BMA LM Averaged vs Best")
par(new=T)
plot(resid.1,col=4,ylim=c(-2.5,2),pch=3,ann=F)
abline(h=0,col=2,lwd=2)
abline(h=3*sqrt(var(resid.a)),col=3,lwd=2)
abline(h=-3*sqrt(var(resid.a)),col=3,lwd=2)
abline(h=3*sqrt(var(resid.1)),col=4,lwd=2,lty=2)
abline(h=-3*sqrt(var(resid.1)),col=4,lwd=2,lty=2)
legend(x="topleft", c("Averaged","Best"), col=c(3,4), pch=c(2,3), ncol = 2, cex = 0.8)

e.a <- mean(Y-Y.a)
e.1 <- mean(Y-Y.reg[,1])
e.2 <- mean(Y-Y.reg[,2])
e.3 <- mean(Y-Y.reg[,3])
e.4 <- mean(Y-Y.reg[,4])
e.a1 <- mean(abs(Y-Y.a))
e.11 <- mean(abs(Y-Y.reg[,1]))
e.21 <- mean(abs(Y-Y.reg[,2]))
e.31 <- mean(abs(Y-Y.reg[,3]))
e.41 <- mean(abs(Y-Y.reg[,4]))
s.a <- sum((Y-Y.a)^2)/(length(Y)-1)
s.1 <- sum((Y-Y.reg[,1])^2)/(length(Y)-1)
s.2 <- sum((Y-Y.reg[,2])^2)/(length(Y)-1)
s.3 <- sum((Y-Y.reg[,3])^2)/(length(Y)-1)
s.4 <- sum((Y-Y.reg[,4])^2)/(length(Y)-1)
ss <- c(s.a,s.1,s.2,s.3,s.4)
ee <- c(e.a,e.1,e.2,e.3,e.4)
ee1 <- c(e.a1,e.11,e.21,e.31,e.41)

ss.r <- ss
ee.r <- ee
ee1.r <- ee1

#plot(ss,ee1,col=c(2,3,4,5,6),xlim=c(min(ss),max(ss)))
#windows()
#plot(ss,ee,col=c(2,3,4,5,6),xlim=c(min(ss),max(ss)))
#abline(h=0)

#cor(Y.a,Y)^2
#cor(Y.reg[,1],Y)^2
#cor(Y.reg[,2],Y)^2
#cor(Y.reg[,3],Y)^2
#cor(Y.reg[,4],Y)^2
#cor(no2.pred,Y)^2

#######################################################
### Comparison between p-value and posterior probability
###      for x-s (shortage of p-value)
#######################################################
no2.1 <- lm(no2~.,no2.dat)
su.FLM <- summary(no2.1)
p.value <- round(su.FLM[[4]][,4],5)
no2.reg <- bicreg(no2.dat[,-1],no2.dat[,1],strict = FALSE, OR = 20)
su.reg <- summary(no2.reg)
PEP.reg <- as.numeric(su.reg[1:8,1])
cpp <- cbind(p.value,PEP.reg)
cpp
plot(p.value[order(p.value)],PEP.reg[order(p.value)],"b",las=1,col=1:4,pch=1:8,ann=F)
title(main="P-value v.s. Post Prob",xlab="P-value",ylab="Post Prob")
legend(x="topright", c("Intercept","day","wdirection"), col=c(2,3,4), pch=c(6,7,8), ncol = 2, cex = 0.8)
grid()


#######################################################
#### Prediction Power Comparison (K folders case)(bicreg)
#######################################################

K.sample <- function(dat,K){    
  #funtion for randomly dividing Data Set into K folders, return an index of observations
    nobs <- nrow(dat)
    set.seed(123)
    K1=floor(nobs/K)
    rand <- c(rep(1:K,K1),sample(1:K,nobs-K*K1))
    rand <- sample(rand)
    rand
}

k <- 50
ind <- K.sample(no2.dat,k)

Y.o <- NULL
cor.r <- NULL
for(i in 1:k){
###iter starts
no2.train <- no2.dat[ind!=i,]
no2.test <- no2.dat[ind==i,]

no2.reg <- bicreg(no2.train[,-1],no2.train[,1],strict = FALSE, OR = 20)

P <- no2.reg[[1]]         #model posterio probability
Beta <- no2.reg[[14]]     #model beta_hat
X <- as.matrix(no2.test[,-1])
X <- cbind(rep(1,nrow(X)),X)
Y <- no2.test[,1]

Y.reg <- X%*%t(Beta)      #predict values for all models
P.reg <- matrix(rep(P,nrow(Y.reg)),nrow(Y.reg),,byrow=T)
Y.p <- Y.reg*P.reg
Y.a <- apply(Y.p,1,sum)   #average predict for all models

Y.o[ind==i] <- Y.a
cor.r <- rbind(cor.r, 
			c(i, cor(Y.a,Y)^2, cor(Y.reg[,1],Y)^2,
			cor(Y.reg[,2],Y)^2, cor(Y.reg[,3],Y)^2 )
		)
###iter ends
}
Y <- no2.dat[,1]

y.lim <- c(min(Y.reg.a,Y.o),max(Y.reg.a,Y.o))
plot(Y,Y.reg.a,col=3,pch=2,ylab="Y_hat",ylim=y.lim,main="BMA LM Full Data vs 50-Fold CV")
par(new=T)
plot(Y,Y.o,col=4,pch=3,ylim=y.lim,ann=F)
abline(c(0,1),col=2,lwd=2)
legend(x="topleft", c("Averaged","50-Fold CV"), col=c(3,4), pch=c(2,3), ncol = 2, cex = 0.8)

resid.a <- Y-Y.reg.a
resid.1 <- Y-Y.o
plot(resid.a,col=3,ylim=c(-2.5,2),pch=2,ylab="Residual",main="BMA LM Full Data vs 50-Fold CV")
par(new=T)
plot(resid.1,col=4,ylim=c(-2.5,2),pch=3,ann=F)
abline(h=0,col=2,lwd=2)
abline(h=3*sqrt(var(resid.a)),col=3,lwd=2)
abline(h=-3*sqrt(var(resid.a)),col=3,lwd=2)
abline(h=3*sqrt(var(resid.1)),col=4,lwd=2,lty=2)
abline(h=-3*sqrt(var(resid.1)),col=4,lwd=2,lty=2)
legend(x="topleft", c("Averaged","50-Fold CV"), col=c(3,4), pch=c(2,3), ncol = 2, cex = 0.8)

colnames(cor.r) <- c("","r2.average","r2.1st","r2.2nd","r2.3rd")
#cor.r
#apply(cor.r,2,mean)

Y <- no2.dat[,1]
e.o <- mean(Y-Y.o)
e1.o <- mean(abs(Y-Y.o))
s.o <- sum((Y-Y.o)^2)/(length(Y)-1)
#cor(Y,Y.o)^2


if(F){
### BACKUP
ss <- c(s.a,s.1,s.2,s.3,s.4,s.o)
ee <- c(e.a,e.1,e.2,e.3,e.4,e.o)
ee1 <- c(e.a1,e.11,e.21,e.31,e.41,e.o1)
plot(ss,ee1,col=c(2,3,4,5,6,7),xlim=c(min(ss),max(ss)))
windows()
plot(ss,ee,col=c(2,3,4,5,6,7),xlim=c(min(ss),max(ss)))
abline(h=0)



plot(cor.r[,1],cor.r[,2],col=2,"b",lwd=2,
	xlab="folder",ylab="r2",ylim=c(min(cor.r[,-1]),max(cor.r[,-1]))
	,main="R2 for Averaging Model and Individual Models")
par(new=T)
plot(cor.r[,1],cor.r[,3],col=3,"b",lty=2,
	xlab="",ylab="",ylim=c(min(cor.r[,-1]),max(cor.r[,-1])))
par(new=T)
plot(cor.r[,1],cor.r[,4],col=4,"b",lty=2,
	xlab="",ylab="",ylim=c(min(cor.r[,-1]),max(cor.r[,-1])))
par(new=T)
plot(cor.r[,1],cor.r[,5],col=5,"b",lty=2,
	xlab="",ylab="",ylim=c(min(cor.r[,-1]),max(cor.r[,-1])))
### BACKUP END
}


#######################################################
#### Prediction Power Comparison (for full data set)(Gamma GLM)
#######################################################
no2.glm <- bic.glm(no2.dat[,-1],no2.dat[,1],strict = FALSE, OR = 20,glm.family="Gamma")
summary(no2.glm)
model.glm <- no2.glm[[7]] #top models
P <- no2.glm[[1]]         #model posterio probability
Beta <- no2.glm[[16]]     #model beta_hat
SE <- no2.glm[[17]]       #model sd

X <- as.matrix(no2.dat[,-1])
X <- cbind(rep(1,nrow(X)),X)
Y <- no2.dat[,1]

Y.glm <- 1/(X%*%t(Beta))      #predict value for all models
P.glm <- matrix(rep(P,nrow(Y.glm)),nrow(Y.glm),,byrow=T)
Y.p <- Y.glm*P.glm
Y.a <- apply(Y.p,1,sum)   #average predict for all models

y.lim <- c(min(Y.a,Y.glm[,1]),max(Y.a,Y.glm[,1]))
plot(Y,Y.a,col=3,pch=2,ylab="Y_hat",ylim=y.lim,main="BMA LM Averaged vs Best")
par(new=T)
plot(Y,Y.glm[,1],col=4,pch=3,ylim=y.lim,ann=F)
abline(c(0,1),col=2,lwd=2)
legend(x="topleft", c("Averaged","Best"), col=c(3,4), pch=c(2,3), ncol = 2, cex = 0.8)

resid.a <- Y-Y.a
resid.1 <- Y-Y.glm[,1]
plot(resid.a,col=3,ylim=c(-2.5,2),pch=2,ylab="Residual",main="BMA LM Averaged vs Best")
par(new=T)
plot(resid.1,col=4,ylim=c(-2.5,2),pch=3,ann=F)
abline(h=0,col=2,lwd=2)
abline(h=3*sqrt(var(resid.a)),col=3,lwd=2)
abline(h=-3*sqrt(var(resid.a)),col=3,lwd=2)
abline(h=3*sqrt(var(resid.1)),col=4,lwd=2,lty=2)
abline(h=-3*sqrt(var(resid.1)),col=4,lwd=2,lty=2)
legend(x="topleft", c("Averaged","Best"), col=c(3,4), pch=c(2,3), ncol = 2, cex = 0.8)

e.a <- mean(Y-Y.a)
e.1 <- mean(Y-Y.glm[,1])
e.2 <- mean(Y-Y.glm[,2])
e.3 <- mean(Y-Y.glm[,3])
e.4 <- mean(Y-Y.glm[,4])
e.a1 <- mean(abs(Y-Y.a))
e.11 <- mean(abs(Y-Y.glm[,1]))
e.21 <- mean(abs(Y-Y.glm[,2]))
e.31 <- mean(abs(Y-Y.glm[,3]))
e.41 <- mean(abs(Y-Y.glm[,4]))
s.a <- sum((Y-Y.a)^2)/(length(Y)-1)
s.1 <- sum((Y-Y.glm[,1])^2)/(length(Y)-1)
s.2 <- sum((Y-Y.glm[,2])^2)/(length(Y)-1)
s.3 <- sum((Y-Y.glm[,3])^2)/(length(Y)-1)
s.4 <- sum((Y-Y.glm[,4])^2)/(length(Y)-1)
ss <- c(s.a,s.1,s.2,s.3,s.4)
ee <- c(e.a,e.1,e.2,e.3,e.4)
ee1 <- c(e.a1,e.11,e.21,e.31,e.41)

ss.g <- ss
ee.g <- ee
ee1.g <- ee1

#windows()
#plot(ss,ee1,col=c(2,3,4,5,6),xlim=c(min(ss),max(ss)))
#windows()
#plot(ss,ee,col=c(2,3,4,5,6),xlim=c(min(ss),max(ss)))
#abline(h=0)


#######################################################
#### Prediction Power Comparison (for full data set)(Normal GLM)
#######################################################
no2.glm <- bic.glm(no2.dat[,-1],no2.dat[,1],strict = FALSE, OR = 20,glm.family="gaussian")
summary(no2.glm)
model.glm <- no2.glm[[7]] #top models
P <- no2.glm[[1]]         #model posterio probability
Beta <- no2.glm[[16]]     #model beta_hat
SE <- no2.glm[[17]]       #model sd

X <- as.matrix(no2.dat[,-1])
X <- cbind(rep(1,nrow(X)),X)
Y <- no2.dat[,1]

Y.glm <- (X%*%t(Beta))      #predict value for all models
P.glm <- matrix(rep(P,nrow(Y.glm)),nrow(Y.glm),,byrow=T)

Y.p <- Y.glm*P.glm
Y.a <- apply(Y.p,1,sum)   #average predict for all models

e.a <- mean(Y-Y.a)
e.1 <- mean(Y-Y.glm[,1])
e.2 <- mean(Y-Y.glm[,2])
e.3 <- mean(Y-Y.glm[,3])
e.4 <- mean(Y-Y.glm[,4])
e.a1 <- mean(abs(Y-Y.a))
e.11 <- mean(abs(Y-Y.glm[,1]))
e.21 <- mean(abs(Y-Y.glm[,2]))
e.31 <- mean(abs(Y-Y.glm[,3]))
e.41 <- mean(abs(Y-Y.glm[,4]))
s.a <- sum((Y-Y.a)^2)/(length(Y)-1)
s.1 <- sum((Y-Y.glm[,1])^2)/(length(Y)-1)
s.2 <- sum((Y-Y.glm[,2])^2)/(length(Y)-1)
s.3 <- sum((Y-Y.glm[,3])^2)/(length(Y)-1)
s.4 <- sum((Y-Y.glm[,4])^2)/(length(Y)-1)
ss <- c(s.a,s.1,s.2,s.3,s.4)
ee <- c(e.a,e.1,e.2,e.3,e.4)
ee1 <- c(e.a1,e.11,e.21,e.31,e.41)

ss.n <- ss
ee.n <- ee
ee1.n <- ee1

#windows()
#plot(ss,ee1,col=c(2,3,4,5,6),xlim=c(min(ss),max(ss)))
#windows()
#plot(ss,ee,col=c(2,3,4,5,6),xlim=c(min(ss),max(ss)))
#abline(h=0)


#######################################################
### Plot of Error Comparison for All Models
#######################################################
x.lim <- c(min(c(ss.n,ss.g,ss.r,s.flm,s.fglm.ga,s.fglm.ig)),max(c(ss.n,ss.g,ss.r,s.flm,s.fglm.ga,s.fglm.ig)))
y.lim <- c(min(c(ee1.n,ee1.g,ee1.r,e1.flm,e1.fglm.ga,e1.fglm.ig)),max(c(ee1.n,ee1.g,ee1.r,e1.flm,e1.fglm.ga,e1.fglm.ig)))
windows()
plot(ss.r,ee1.r,col=c(2,3,4,5,6),xlim=x.lim,ylim=y.lim,las=1,
	xlab="square error mean",ylab="absolute error mean",main="Error Comparison for All Models")
points(ss.g,ee1.g,pch=2,col=c(2,3,4,5,6),xlim=x.lim,ylim=y.lim)
#points(ss.n,ee1.n,pch=3,col=c(2,3,4,5,6),xlim=x.lim,ylim=y.lim)
points(s.o,e1.o,pch=5)
points(s.flm,e1.flm,pch=4)
points(s.fglm.ga,e1.fglm.ga,pch=6)
points(s.fglm.ig,e1.fglm.ig,pch=7)
legend(x="topleft", bty="n", pch=c(1,2,5,4,6,7),
       legend=c("BMA LM (5)", "BMA GLM(gamma) (5)","BMA LM 50-Fold CV","LM","GLM(gamma)","GLM(inverse gamma)"))




