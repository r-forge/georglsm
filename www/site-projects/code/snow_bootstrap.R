

######################
### snow for Bootstrap
######################

require(snow)
library(boot)
#  In this example we show the use of boot in a prediction from 
#  regression based on the nuclear data.  This example is taken 
#  from Example 6.8 of Davison and Hinkley (1997).  Notice also 
#  that two extra arguments to statistic are passed through boot.
data(nuclear)
nuke <- nuclear[,c(1,2,5,7,8,10,11)]
nuke.lm <- glm(log(cost)~date+log(cap)+ne+ ct+log(cum.n)+pt, data=nuke)
nuke.diag <- glm.diag(nuke.lm)
nuke.res <- nuke.diag$res*nuke.diag$sd
nuke.res <- nuke.res-mean(nuke.res)

#  We set up a new dataframe with the data, the standardized 
#  residuals and the fitted values for use in the bootstrap.
nuke.data <- data.frame(nuke,resid=nuke.res,fit=fitted(nuke.lm))

#  Now we want a prediction of plant number 32 but at date 73.00
new.data <- data.frame(cost=1, date=73.00, cap=886, ne=0,
                       ct=0, cum.n=11, pt=1)
new.fit <- predict(nuke.lm, new.data)

nuke.fun <- function(dat, inds, i.pred, fit.pred, x.pred) {
     assign(".inds", inds, envir=.GlobalEnv)
     lm.b <- glm(fit+resid[.inds] ~date+log(cap)+ne+ct+
                 log(cum.n)+pt, data=dat)
     pred.b <- predict(lm.b,x.pred)
     remove(".inds", envir=.GlobalEnv)
     c(coef(lm.b), pred.b-(fit.pred+dat$resid[i.pred]))
}


## sequential version
R <- 1000
system.time(nuke.boot <-
            boot(nuke.data, nuke.fun, R=R, m=1,
                 fit.pred=new.fit, x.pred=new.data))
# elapsed: 19.392

## parallel version
cl <- makeCluster(4, "MPI")
clusterEvalQ(cl,library(boot))
clusterSetupRNG(cl)
system.time(cl.nuke.boot <-
            clusterCall(cl,boot,nuke.data, nuke.fun, R=R/length(cl), m=1,
                        fit.pred=new.fit, x.pred=new.data))
stopCluster(cl)
# elapsed: 7.123
class(cl.nuke.boot); length(cl.nuke.boot)
class(cl.nuke.boot[[1]]); names(cl.nuke.boot[[1]])
dim(cl.nuke.boot[[1]]$t)
cl.nuke.boot[[2]]$pred.i
## Merging the results:
# The need to write a merging function is a common pattern 
# in parallelizing computations using the intermediate-level 
# functions clusterCall() and clusterApply(). 
fixboot <- function(bootlist) {
    # let returning list have the same structure as lists from workers
    boot <- bootlist[[1]] 
    # merge attribute "t" among the results of all workers
    boot$t <- do.call(rbind,lapply(bootlist, function(x) x$t))
    # merge attribute "R"
    boot$R <- sum(sapply(bootlist, function(x) x$R))
    # merge attribute "pred.i" if exist
    if (! is.null(boot$pred.i))
        boot$pred.i <- do.call(rbind,lapply(bootlist,
                                            function(x) x$pred.i))
    # other attributes have same values as 1st worker's result
    boot
}
cl.nuke.boot.fixed <- fixboot(cl.nuke.boot)
class(cl.nuke.boot.fixed)
dim(cl.nuke.boot.fixed$t)
cl.nuke.boot.fixed$pred.i

##  The bootstrap prediction error would then be found by
mean(nuke.boot$t[,8]^2)
mean(cl.nuke.boot.fixed$t[,8]^2)

##  Basic bootstrap prediction limits would be
new.fit-sort(nuke.boot$t[,8])[c(975,25)]
new.fit-sort(cl.nuke.boot.fixed$t[,8])[c(975,25)]