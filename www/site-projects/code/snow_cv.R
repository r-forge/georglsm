
#############################
### snow for cross-validation
#############################

require(snow)

data(cars)
str(cars)
fitlm <- lm(speed ~ ., cars)
summary(fitlm)
plot(cars$dist, cars$speed)
abline(fitlm$coef, col=2)

n <- nrow(cars)
k <- 8
d <- ceiling(n/k)
tmp <- sample(1:n)
ind <- lapply(0:(d-1), function(x) tmp[(x*k+1):min(x*k+k,n)] )
# Step 1: write code draft by using for loop.
res <- numeric(d)
for(i in 1:d){
fitlm <- lm(speed~., cars[-ind[[i]],])
res[i] <- sum(fitlm$resid^2)
}
# Step 2: rewrite code in terms of a sequential mapping function such as lapply().
res2 <- unlist(lapply(ind, function(x){
  fitlm <- lm(speed~.,cars[-x,])
  sum(fitlm$resid^2)
  }
))
res2 <- unlist(lapply(1:d, function(i){ 
  fitlm <- lm(speed~., cars[-ind[[i]],])
  sum(fitlm$resid^2)
  }
))
# Step 3: turn the code into a parallel function by adding a cluster argument and replacing the call to lapply by a call to parLapply. 
cl <- makeCluster(4, "MPI")
res3 <- unlist( parLapply(cl, ind, function(x){
  fitlm <- lm(speed~.,cars[-x,])
  sum(fitlm$resid^2)
  }
))
stopCluster(cl)
# object 'ind' should be included in function.
cl <- makeCluster(4, "MPI")
res4 <- unlist( parLapply(cl, 1:d, function(i, ind){ 
  fitlm <- lm(speed~., cars[-ind[[i]],])
  sum(fitlm$resid^2)
  }, ind
))
stopCluster(cl)

gc(verbose = getOption("verbose"), reset=TRUE)
gc()


