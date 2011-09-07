
########################################################
### 1.satisfying rating ~ time + expense + gender 
###   + age + education + family size + income 
########################################################
type <- factor(trans[,2])
gender <- factor(trans[,ind.col(17)])
inc <- trans[,ind.col(22)]
inc[!is.na(inc) & inc==12] <- 0
inc <- inc*500
dat1 <- data.frame(rate=rate.t,trans[,c(3,4,sapply(c(19,20),ind.col))],
gender=gender,age=age,income=inc)
isna <- apply(is.na(dat1),1,sum)
tt <- dat1[isna==0,]
dim(tt)
temp <- div(tt,0.7)
d1.train <- temp[[1]]
d1.test <- temp[[2]]

#################
### Training data
#################
#LM Modeling
library(MASS)
d1.lm <-lm(rate~.,data=d1.train)
pred <- predict(d1.lm,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.1882451
mean((d1.train$rate-pred)^2) #[1] 173.2422

d1.lms <-stepAIC(d1.lm,trace=F)
#summary(d1.lms)
pred <- predict(d1.lms,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.1878193
mean((d1.train$rate-pred)^2) #[1] 173.3331

d1.lm2 <-lm(rate~.^2,data=d1.train)
pred <- predict(d1.lm2,d1.train)
cor(d1.train$rate,pred)^2 #0.2615466
mean((d1.train$rate-pred)^2) #[1] 157.5984

d1.lm2s <-stepAIC(d1.lm2,trace=F)
pred <- predict(d1.lm2s,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.2461347
mean((d1.train$rate-pred)^2) #[1] 160.8876

#GLM Modeling
d1.glm <-glm(rate~.,family=poisson,data=d1.train)
pred <- predict(d1.glm,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.1874426
mean((d1.train$rate-pred)^2) #[1] 1395.557

d1.glms <-stepAIC(d1.glm,trace=F)
summary(d1.glms)
pred <- predict(d1.glms,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.1874426
mean((d1.train$rate-pred)^2) #[1] 1395.557

d1.glm2 <-glm(rate~.^2,family=poisson,data=d1.train)
pred <- predict(d1.glm2,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.2586693
mean((d1.train$rate-pred)^2) #[1] 1395.178

d1.glm2s <-stepAIC(d1.glm2,trace=F)
pred <- predict(d1.glm2s,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.2586693
mean((d1.train$rate-pred)^2) #[1] 1395.178

#GAM Modeling
library(splines)
library(gam)
d1.gam <- gam(rate~lo(Q2)+lo(Q3)+lo(Q19)+lo(Q20)+gender+lo(age)
              +lo(income),data=d1.train)
pred <- predict(d1.gam,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.2827699
mean((d1.train$rate-pred)^2) #[1] 153.1906

#TREE Modeling
library(rpart)
d1.tree <- rpart(rate~.,d1.train,method="anova",cp=0.001)
printcp(d1.tree)
pred <- predict(d1.tree,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.3919258
mean((d1.train$rate-pred)^2) #[1] 129.7733              
              
d1.treep <- prune(d1.tree,cp=0.02)
pred <- predict(d1.treep,d1.train)
cor(d1.train$rate,pred)^2 #[1] 0.2423955
mean((d1.train$rate-pred)^2) #[1] 161.6856
              
#NNET Modeling
library(nnet)
pred <- 0
iter <- 20
for(i in 1:iter){
d1.nnet <- nnet(rate~.,d1.train,size=10,decay=0.0001,linout=T,
                skip=T,maxit=1000,Hess=T)
pred <- pred+predict(d1.nnet,d1.train)
}
pred <- pred/iter
cor(d1.train$rate,pred)^2 #[1,] 0.266148
mean((d1.train$rate-pred)^2) #[1] 157.7281
                
#################
### Testing data
#################
### LM
pred <- predict(d1.lm,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

pred <- predict(d1.lms,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

pred <- predict(d1.lm2,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

pred <- predict(d1.lm2s,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

### GLM
pred <- predict(d1.glm,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

pred <- predict(d1.glms,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

pred <- predict(d1.glm2,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

pred <- predict(d1.glm2s,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

### GAM
pred <- predict(d1.gam,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

### Tree
pred <- predict(d1.tree,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

pred <- predict(d1.treep,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

### NNET
pred <- predict(d1.nnet,d1.test)
cor(d1.test$rate,pred)^2
mean((d1.test$rate-pred)^2)

#####################
### Cross-validation
#####################
CV.d1 <- function(dat,K)
{
    res <- dat$rate
    set.seed(111)
    K1=floor(length(res)/K)
    rand=c(rep(1:K,K1),sample(1:K,length(res)-K*K1))
    rand=sample(rand)
  
  fitfn<-function(x) lm(rate~., data=dat[x,])    
#  fitfn<-function(x) lm(rate~.^2, data=dat[x,])
#  fitfn<-function(x) gam(rate ~ lo(Q2) + lo(Q3) + lo(Q19) + lo(Q20) + gender + lo(age) , data=dat[x,])    
#  fitfn<-function(x) rpart(rate~.,data=dat[x,],method="anova",cp=0.001)
    
  predfn<-function(obj, x) predict(obj, dat[x, ])

    for (i in sort(unique(rand))) {
        cat("fold ", i, "\n", sep = "")
        learn <- fitfn(rand != i)
        res[rand == i] <- predfn(learn, rand == i)
    }
	c(cor(dat$rate,res)^2, mean((dat$rate-res)^2))
}
CV.d1(tt,10)
#lm(~.):   [1]   0.1235876 196.0953042
#lm(~.^2): [1]   0.096604 213.077082    
#gam:      [1]   0.07793132 215.92547026               
#TREE:     [1]   0.1064169 225.1646964
