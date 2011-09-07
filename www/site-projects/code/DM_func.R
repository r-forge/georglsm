

#########################
### Customized functions 
#########################

ind.col <- function(colname){
#find the column number for column with name 'colname'
temp <- names(trans)
ind <- 1:length(temp)
colname <- paste('Q',colname,sep='')
ind[temp==colname]
}

div <- function(dat,p){
#divide 'dat' into training set as p and test set as (1-p)
nr <- nrow(dat)
set.seed(123)
dat.sample <- c(rep(0,floor(p*nr)),rep(1,nr-floor(p*nr)))
dat.sample <- sample(dat.sample,nr)
dat.train <- dat[dat.sample==0,]
dat.test <- dat[dat.sample==1,]
list(dat.train,dat.test)
}

con <- function(...)   {
#confusion matrix
    print(tab <- table(...))
    diag(tab) <- 0
    cat("error rate = ",
        round(100*sum(tab)/length(list(...)[[1]]), 2), "%\n")
    invisible()
}

CV.d3 <- function(dat,K){
#K cross-validation in 'dat' for data set 3
    res <- dat$loss
    set.seed(111)
    K1=floor(length(res)/K)
    rand=c(rep(1:K,K1),sample(1:K,length(res)-K*K1))
    rand=sample(rand)
  
	fitfn<-function(x) rpart(loss~.,data=dat[x,],method="class",cp=0.001)
	predfn<-function(obj, x) predict(obj, dat[x, ],type="class")

    for (i in sort(unique(rand))) {
        cat("fold ", i, "\n", sep = "")
        learn <- fitfn(rand != i)
        res[rand == i] <- predfn(learn, rand == i)
    }
	res
}

BS <- function(dat,n){
#Bootstrap 'dat' n times
	set.seed(101)
	ind <- sample(1:nrow(dat),n,replace=T)
	dat[ind,]
}
#### END OF FUNCTIONS ########