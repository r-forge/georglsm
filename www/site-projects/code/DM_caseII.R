

########################################################
### 2.driving habits ~ gender + age + education + income
########################################################
dat2 <- data.frame(dhabit=score.dh,gender=gender,age=age,
        edu=trans[,ind.col(19)],income=inc)
isna <- apply(is.na(dat2),1,sum)
tt <- dat2[isna==0,]
dim(tt)
temp <- div(tt,0.7)
d2.train <- temp[[1]]
d2.test <- temp[[2]]

#################
### Training data
#################
#LM Modeling
d2.lm <- lm(dhabit~.,d2.train)
pred <- predict(d2.lm,d2.train)
cor(d2.train$dhabit,pred)^2 #[1] 0.02133738
mean((d2.train$dhabit - pred)^2) #[1] 436.95
                                                     
d2.lms <- stepAIC(d2.lm,trace=F)
pred <- predict(d2.lms,d2.train)
cor(d2.train$dhabit,pred)^2 #0.01204088
mean((d2.train$dhabit - pred)^2) #[1] 441.1007
                   
d2.lm2 <- lm(dhabit~.^2,d2.train)
pred <- predict(d2.lm2,d2.train)
cor(d2.train$dhabit,pred)^2 #[1] 0.05917191
mean((d2.train$dhabit - pred)^2) #[1] 420.0578
                   
#GLM Modeling
d2.glm1 <-glm(dhabit~.,family=Gamma,data=d2.train)
pred <- predict(d2.glm1,d2.train)
cor(d2.train$dhabit,pred)^2 #[1] 0.02133034
mean((d2.train$dhabit - pred)^2) #[1] 4092.436

d2.glm2 <-glm(dhabit~.^2,family=poisson,data=d2.train)
pred <- predict(d2.glm2,d2.train)
cor(d2.train$dhabit,pred)^2 #[1] 0.05887669
mean((d2.train$dhabit - pred)^2) #[1] 3615.434
                   
d2.glm2s <-stepAIC(d2.glm2,trace=F)
pred <- predict(d2.glm2s,d2.train)
cor(d2.train$dhabit,pred)^2 #[1] 0.05887669
mean((d2.train$dhabit - pred)^2) #[1] 3615.434
                   
#GAM Modeling
d2.gam <- gam(dhabit~lo(age)+lo(edu)+lo(income)+gender,data=d2.train)
pred <- predict(d2.gam,d2.train)
cor(d2.train$dhabit,pred)^2 #[1] 0.06865203
mean((d2.train$dhabit - pred)^2) #[1] 416.9481
                   
#TREE Modeling
d2.tree <- rpart(dhabit~.,d2.train,method="anova",cp=0.0001)
printcp(d2.tree)
pred <- predict(d2.tree,d2.train)
cor(d2.train$dhabit,pred)^2 #[1] 0.1865511
mean((d2.train$dhabit - pred)^2) #[1] 363.186
                   
d2.treep <- prune(d2.tree,cp=0.03)
pred <- predict(d2.treep,d2.train)
cor(d2.train$dhabit,pred)^2 #[1] 0.031824
mean((d2.train$dhabit - pred)^2) #[1] 432.268
                   
#NNET Modeling
library(nnet)
pred <- 0
iter <- 20
for(i in 1:iter){
d2.nnet <- nnet(dhabit~.,d2.train,size=9,decay=0.1,linout=T,
                skip=T,maxit=1000,Hess=T)
pred <- pred+predict(d2.nnet,d2.train)
}
pred <- pred/iter
cor(d2.train$dhabit,pred)^2 #0.2990196
mean((d2.train$dhabit - pred)^2) #[1] 343.3025


#################
### Testing data
#################
### LM
pred <- predict(d2.lm,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)
pred <- predict(d2.lms,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)
pred <- predict(d2.lm2,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)

### GLM
pred <- predict(d2.glm2,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)
pred <- predict(d2.glm2s,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)

### GAM
pred <- predict(d2.gam,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)

### Tree
pred <- predict(d2.tree,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)
pred <- predict(d2.treep,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)

### NNET
pred <- predict(d2.nnet,d2.test)
cor(d2.test$dhabit,pred)^2; mean((d2.test$dhabit - pred)^2)



