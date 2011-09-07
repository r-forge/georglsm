
rm(list=ls(all=T))
##################

# setwd("/...")

#### Data Initializing
######################
gr.dat <- read.table("eye tracking.txt")
id.f <- scan("eye tracking_idfemale.txt")

colnames(gr.dat) <- c("id","type",paste("cycle", 1:11, sep = ""))
id.blink <- which(apply(gr.dat==0,1,sum)!=0)
gr.dat <- gr.dat[-id.blink,]
#gr.dat$id
sex <- rep(0,length(gr.dat$id))
for(i in 1:length(gr.dat$id)){
	sex[i] <- sum(gr.dat$id[i]==id.f)
}
gr.dat <- cbind(gr.dat,sex)

nrow(gr.dat)             #60
length(unique(gr.dat$id))  #[1] 31
length(id.f)             #22
sum(gr.dat$type=="TR")   #[1] 9
sum(gr.dat$type=="CS")   #[1] 17
sum(gr.dat$type=="PS")   #[1] 34

### Three moving types
######################
temp <- seq(0,3*pi,length=300)
plot(temp,sin(temp),"l",axes=F,xlab="time",ylab="velocity of target",main="PS/CS")
abline(h=0)
plot(c(1,2,2,3,3,4),c(1,1,-1,-1,1,1),axes=F,"l",xlab="time",ylab="velocity of target",main="TR")
abline(h=0)

### First Subject
#################
k <- 3
y.lim <- c(min(gr.dat[1:k,3:13]),max(gr.dat[1:k,3:13]))
for(i in 1:k){
plot(1:11,gr.dat[i,3:13],"b",col=i,ylim=y.lim,ann=F)
par(new=T)
}
title(main="First Patient",xlab="Cycle",ylab="Gain Ratio")


### Data Transformation
#######################
temp <- gr.dat$id
temp <- t(matrix(rep(temp,11),length(gr.dat$id),))
library(gdata)
id <- unmatrix(temp)
gr <- unmatrix(t(as.matrix(gr.dat[,3:13])))
cycle <- rep(1:11,60)
sex <-  unmatrix(t(matrix(rep(gr.dat$sex,11),60,)))
type <-  unmatrix(t(matrix(rep(gr.dat$type,11),60,)))
eye.dat <- data.frame(id=id,sex=sex,type=type,cycle=cycle,gr=gr)

library(lattice)
xyplot(gr~type|sex,eye.dat,group=id,type="p")
xyplot(gr~cycle|sex,eye.dat[1:121,],group=id,type="l")
dotplot(gr~type|sex,eye.dat,group=id)
stripplot(gr~type,eye.dat,group=id)
barchart(gr~type,eye.dat,group=id)
xyplot(gr~cycle|type*sex,eye.dat,group=id,type="l")
bwplot(gr~type|sex,eye.dat,group=id)


### LM
######
gr.lm1 <- lm(gr~sex+type,eye.dat)
#summary(gr.lm1)
gr.lm2 <- lm(gr~sex*type,eye.dat)
gr.lm3 <- lm(gr~sex*type+cycle,eye.dat)
gr.lm4 <- lm(gr~sex+type*cycle,eye.dat)
gr.lm5 <- lm(gr~type+sex*cycle,eye.dat)
gr.lm6 <- lm(gr~type*sex*cycle,eye.dat)
aic.lm <- AIC(gr.lm1,gr.lm2,gr.lm3,gr.lm4,gr.lm5,gr.lm6)
llk.lm1 <- logLik(gr.lm1)
llk.lm2 <- logLik(gr.lm2)
llk.lm3 <- logLik(gr.lm3)
llk.lm4 <- logLik(gr.lm4)
llk.lm5 <- logLik(gr.lm5)
llk.lm6 <- logLik(gr.lm6)
lik.lm <- c(llk.lm1,llk.lm2,llk.lm3,llk.lm4,llk.lm5,llk.lm6)
#plot(1:6,aic.lm[,2],"b",col=3)


yhat.lm3 <- predict(gr.lm3)
yhat.lm6 <- predict(gr.lm6)
cor(eye.dat$gr,yhat.lm3)
cor(eye.dat$gr,yhat.lm6)

xy.lim <- c(min(eye.dat$gr,yhat.lm6),max(eye.dat$gr,yhat.lm6))
plot(eye.dat$gr,yhat.lm6,xlim=xy.lim,ylim=xy.lim,ann=F)
title(xlab="True Y",ylab="Predicted Y")
abline(coef=c(0,1),col=2,lwd=2)

plot(resid(gr.lm6),ylim=c(-0.4,0.5))
abline(h=0,col=2,lwd=2)
abline(h=3*sqrt(var(resid(gr.lm6))),col=3,lwd=2)
abline(h=-3*sqrt(var(resid(gr.lm6))),col=3,lwd=2)


### Mixed Models
################

### nlme
########
library(nlme)
gr.nm1 <- lme(gr~sex+type,data=eye.dat,random=~(1)|id,method="ML")
summary(gr.nm1)
gr.nm2 <- lme(gr~sex+type,data=eye.dat,random=~(cycle-1)|id,method="ML")
summary(gr.nm2)
gr.nm3 <- lme(gr~sex+type,data=eye.dat,random=~(cycle)|id,method="ML")
summary(gr.nm3)
gr.nm4 <- lme(gr~sex*type,data=eye.dat,random=~(1)|id,method="ML")
summary(gr.nm4)
gr.nm5 <- lme(gr~sex*type,data=eye.dat,random=~(cycle-1)|id,method="ML")
summary(gr.nm5)
gr.nm6 <- lme(gr~sex*type,data=eye.dat,random=~(cycle)|id,method="ML")
summary(gr.nm6)

gr.nm1 <- lme(gr~type,data=eye.dat,random=~(1)|id,method="ML")
summary(gr.nm1)
gr.nm2 <- lme(gr~type,data=eye.dat,random=~(cycle-1)|id,method="ML")
summary(gr.nm2)
gr.nm3 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML")
summary(gr.nm3)


gr.nm1 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corExp())
summary(gr.nm1)
gr.nm2 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corCompSymm())
summary(gr.nm2)
gr.nm3 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corAR1())
summary(gr.nm3)
gr.nm4 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corGaus())
summary(gr.nm4)
gr.nm5 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corLin())
summary(gr.nm5)
gr.nm6 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corSpher())
summary(gr.nm6)
#gr.nm6 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corSymm())
#summary(gr.nm6)


gr.nm1 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML")
summary(gr.nm1)
gr.nm2 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corAR1())
summary(gr.nm2)


#mcom <- read.table("D:/Study/STA 7723/Final/AIC.txt")
mcom <- read.table("X:/liang.jing/STA 7723/Final/AIC.txt")

y.lim <- c(min(aic.lm[,2],mcom[,2]),max(aic.lm[,2],mcom[,2]))
plot(1:7,mcom[,2],"b",col=2,pch=c(rep(1,3),rep(2,4)),ann=F,ylim=y.lim)
title(main="",xlab="model",ylab="AIC")
points(1:6,aic.lm[,2],"b",col=3,pch=4)

y.lim <- c(min(lik.lm,mcom[,1]),max(lik.lm,mcom[,1]))
plot(1:7,mcom[,1],"b",col=2,pch=c(rep(1,3),rep(2,4)),ann=F,ylim=y.lim)
title(main="",xlab="model",ylab="Loglik")
points(1:6,lik.lm,"b",col=3,pch=4)


gr.ar1 <- lme(gr~type,data=eye.dat,random=~(cycle)|id,method="ML",corr=corAR1())
VarCorr(gr.ar1)
yhat.ar1 <- fitted(gr.ar1)
res.ar1 <- resid(gr.ar1)
cor(eye.dat$gr,yhat.ar1)
e <- c(mean(abs(res.ar1)),mean(abs(resid(gr.lm3))),mean(abs(resid(gr.lm6))))
s <- c(mean(res.ar1^2),mean(resid(gr.lm3)^2),mean(resid(gr.lm6)^2))
plot(s,e,col=1:3,pch=1:3,cex=2,type="b")

xy.lim <- c(min(eye.dat$gr,yhat.ar1),max(eye.dat$gr,yhat.ar1))
plot(eye.dat$gr,yhat.ar1,xlim=xy.lim,ylim=xy.lim,ann=F)
title(xlab="True Y",ylab="Predicted Y")
abline(coef=c(0,1),col=2,lwd=2)

plot(resid(gr.ar1),ylim=c(-0.4,0.5))
abline(h=0,col=2,lwd=2)
abline(h=3*sqrt(var(resid(gr.ar1))),col=3,lwd=2)
abline(h=-3*sqrt(var(resid(gr.ar1))),col=3,lwd=2)

bwplot(yhat.ar1~type|sex,eye.pre,group=id)
bwplot(res.ar1~type|sex,eye.pre,group=id)
xyplot(res.ar1~cycle|type,eye.pre,group=id,type="l")
xyplot(res.ar1~cycle|sex,eye.pre,group=id,type="l")




### lme4
########
library(lme4)
gr.mm0 <- lmer(gr~sex+type+(1|id),eye.dat)
#summary(gr.mm0)
gr.mm1 <- lmer(gr~sex+type+(cycle-1|id),eye.dat)
gr.mm2 <- lmer(gr~sex+type+(cycle|id),eye.dat)
gr.mm3 <- lmer(gr~sex+type+(cycle|id)+(sex|id),eye.dat)
gr.mm4 <- lmer(gr~sex+type+(sex-1|id),eye.dat)
gr.mm5 <- lmer(gr~sex+type+(type-1|id),eye.dat)
gr.mm6 <- lmer(gr~type+(sex-1|id),eye.dat)
anova(gr.mm0,gr.mm1,gr.mm2,gr.mm3,gr.mm4,gr.mm5,gr.mm6)
yhat.mm6 <- fitted(gr.mm6)
plot(resid(gr.mm6))
abline(h=0,col=2)

e <- NULL
s <- NULL

temp <- gr.mm6
e <- c(e,mean(abs(resid(temp))))
s <- c(s,sum(resid(temp)^2)/(nrow(eye.dat)-1))

e.abs <- e
s.2 <- s

plot(e.abs,s.2,col=1:2,pch=1:2)
text(e.abs,s.2,labels=as.character(0:6),adj=c(0,0.002))
#library(calibrate)
#textxy(e.abs,s.2,labs=1:7,cx=1)

lmer(gr~sex+type+(cycle|id),eye.dat)
library(RLRsim)
extract.lmerDesign(gr.mm3)
















