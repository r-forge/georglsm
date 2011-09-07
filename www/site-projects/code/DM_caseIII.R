
###############################################################
### 3.accident loss ~ driving habits + driving time + frequency
###############################################################
dat3 <- data.frame(loss=as.factor(acc),habit=score.dh,time=trans[,ind.col(11)],
freq=trans[,ind.col(12)])
isna <- apply(is.na(dat3),1,sum)
tt <- dat3[isna==0,]
dim(tt)
temp <- div(tt,0.7)
d3.train <- temp[[1]]
d3.test <- temp[[2]]

#################
### Training data
#################
### TREE 
d3.tree <- rpart(loss~.,d3.train,method="class",cp=0.001)
pred <- predict(d3.tree,d3.train,type = "class")
con(true=d3.train$loss,predicted=pred) # error rate =  40.19 %
printcp(d3.tree)
d3.treep <- prune(d3.tree,cp=0.01)
pred <- predict(d3.treep,d3.train,type = "class")
con(true=d3.train$loss,predicted=pred) # error rate =  41.12 %

pred <- CV.d3(d3.train,10)
con(true=d3.train$loss,predicted=pred) # error rate =  55.14 %

### NNET 
pred <- 0
iter <- 20
for(i in 1:iter){
d3.nnet <- nnet(loss~.,d3.train,size=9,decay=0.001,entropy=T,maxit=1000)
pred <- pred+predict(d3.nnet,d3.train,type="raw")
}
pred <- pred/iter
dim.p <- dim(pred)
pm <- apply(pred,1,max)
rr <- rep(pm,dim.p[2])
res <- matrix(rr,dim.p[1],dim.p[2])
res <- pred-res
pred.n <- apply(res,1,function(x){ names(x[x==0]) })
tru <- as.character(d3.train$loss)
( error <- 1-sum(tru==pred.n)/length(tru) ) # error rate = [1] 0.2336449

#################
### Testing data
#################
### TREE
pred <- predict(d3.tree,d3.test,type = "class")
con(true=d3.test$loss,predicted=pred) # error rate =  68.09 %
pred <- predict(d3.treep,d3.test,type = "class")
con(true=d3.test$loss,predicted=pred) # error rate =  57.45 %

pred <- CV.d3(d3.test,10)
con(true=d3.test$loss,predicted=pred) # error rate =  53.19 %

### NNET
pred <- 0
iter <- 20
for(i in 1:iter){
d3.nnet <- nnet(loss~.,d3.test,size=9,decay=0.001,entropy=T,maxit=1000)
pred <- pred+predict(d3.nnet,d3.test,type="raw")
}
pred <- pred/iter
dim.p <- dim(pred)
pm <- apply(pred,1,max)
rr <- rep(pm,dim.p[2])
res <- matrix(rr,dim.p[1],dim.p[2])
res <- pred-res
pred.n <- apply(res,1,function(x){
names(x[x==0])
})
tru <- as.character(d3.test$loss)
error <- 1-sum(tru==pred.n)/length(tru) # error rate = [1] 0.06382979

#################
### Bootstrapping
#################
### TREE
for(nBS in c(500, 1000, 1500, 2000, 2500, 3000)){
  d3.trbs <- BS(tt, nBS)
  pred <- CV.d3(d3.trbs,10)
  con(true=d3.trbs$loss,predicted=pred)
} # error rate =  35 %
