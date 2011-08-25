
####################################
#### Calculate transformed residuals
####################################
cdfU <- function(Y.obs,Y.rep, discrete=TRUE){
  n <- length(Y.obs)
  nc <- ncol(Y.rep)
  Y.tmp <- matrix(rep(Y.obs,times=nc),,nc)
  if(discrete){
    res <- rowMeans(Y.rep <= (Y.tmp-1)) + runif(n)*rowMeans(Y.rep == Y.tmp)
  } else res <- rowMeans(Y.rep <= Y.tmp)
  res
}
tranR <- function(Y.obs,Y.rep, discrete=TRUE){
  u <- cdfU(Y.obs,Y.rep, discrete)
  u[u<1e-04] <- 1e-04; u[u>(1-1e-04)] <- 1-1e-04
  qnorm(u)
}
####################################
#### Plot transformed residuals
####################################
plot_etran <- function(e.tran, fig=1:4)
  {
e.tran <- e.tran[e.tran>-1e5 & e.tran<1e5]
if(length(fig)==4){
  op <- par(mfrow=c(2,2))
} else if(length(fig)==3){
  op <- par(mfrow=c(3,1))
} else op <- par(mfrow=c(1,length(fig)))

if(any(fig==1)) plot(e.tran)
if(any(fig==2)){
  qqnorm(e.tran); qqline(e.tran,col=2)
}
if(any(fig==3)){
  plot(density(e.tran))
  lines(x <- seq(-5,5,by=0.1),dnorm(x),col=2)
}
if(any(fig==4)){
  hist(e.tran,freq=FALSE,ylim=c(0,dnorm(0)+0.07))
  lines(x <- seq(-5,5,by=0.1),dnorm(x),col=2)
}
par(op)
}
####################################
#### Calculate three distances
####################################
e2dist <- function(e.tran)
{
    x <- e.tran
    x <- x[x>-1e5 & x<1e5]
    res <- c(HellingerDist(Norm(), x),HellingerDist(x, Norm(), asis.smooth.discretize = "smooth"),KolmogorovDist(Norm(),x) )
    names(res) <- c("Hellinger.discre","Hellinger.smooth","Kolmogorov")
    res
  }
####################################
#### Generate baseline of distances
####################################
baseline.dist <- function(n, iter){
res <- matrix(0, iter, 3)
for(i in 1:iter){
x <- rnorm(n)
res[i,] <- c(HellingerDist(x, Norm()),HellingerDist(x, Norm(), asis.smooth.discretize = "smooth"),KolmogorovDist(x, Norm()) )
}
colnames(res) <- c("Hellinger.discre","Hellinger.smooth","Kolmogorov")
res
}
baseline.parallel <- function(n, iter){
res0 <- mclapply(1:iter, function(t){
        x <- rnorm(n)
        c(HellingerDist(x, Norm()),HellingerDist(x, Norm(), asis.smooth.discretize = "smooth"),KolmogorovDist(x, Norm()) )
      }
    )
res <- matrix(unlist(res0),,3,byrow=TRUE)
colnames(res) <- c("Hellinger.discre","Hellinger.smooth","Kolmogorov")
res
}
####################################
#### Plot baseline of distances
####################################
plot_baseline <- function(res.in, dist.name){
val.crit <- quantile(res.in,prob=c(0.025,0.5,0.975))
plot(density(res.in),main=paste("Baseline of ",dist.name," Distance",sep=""))
abline(v=val.crit,lty=2,col=c(3,2,3))
#text(val.crit,0,paste(c("2.5%","50%","97.5%"),round(val.crit,4),sep="="),col=c(3,2,3),cex=0.8)
legend("topright",paste(c("2.5%","50%","97.5%"),round(val.crit,4),sep="="),lty=2,col=c(3,2,3))
}
####################################
#### Calculate one-side p-value
####################################
pOne <- function(d.obs, d.base){
  if(length(d.obs)==1){
  p <- mean(d.base<=d.obs)
  p <- ifelse(p<=0.5, p, 1-p)
  } else{
      tmp <- rbind(d.obs, d.base)
      p <- apply(tmp, 2, function(t){
              pp <- mean(t[-1]<=t[1]) 
              ifelse(pp<=0.5, pp, 1-pp)
              } 
            )
    }
  p
  }
####################################
#### END
####################################