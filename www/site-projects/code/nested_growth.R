
##### Cell growth program for Nested Sampling ####
#In 2-D parameter space, expand from a starting point 
#until all points on boundary reach threshold likelihood.


###################
rm(list=ls(all=T))
###################

###### Functions (run once)#######
LL <- function(mu,sigma){sum(log(dnorm(dat,mu,sigma)))}
#log-likelihood function

walk <- function(p){
p.x <- p.nx <- p.y <- p.ny <- NULL
if(p$x==0) p.x <- c(x=0,y=0,nx=1,ny=0,a=p$a+da,b=p$b)
if(p$nx==0) p.nx <- c(x=1,y=0,nx=0,ny=0,a=p$a-da,b=p$b)
if(p$y==0) p.y <- c(x=0,y=0,nx=0,ny=1,a=p$a,b=p$b+db)
if(p$ny==0) p.ny <- c(x=0,y=1,nx=0,ny=0,a=p$a,b=p$b-db)
as.data.frame(rbind(p.x,p.nx,p.y,p.ny))
}

expand <- function(dat){
res <- NULL
for(i in 1:nrow(dat)){
	res <- rbind(res,walk(dat[i,]))
	}
res <- aggregate(res[,c('x','y','nx','ny')],res[,c('a','b')],sum)
res <- data.frame(res,step=dat$step[1]+1,open=0)
res
}

updt.IB <- function(I,B,L.m){
nB <- expand(B[B$open==0,])
L.res <- NULL
for(i in 1:nrow(nB)){
	L.res <- c(L.res,LL(nB$a[i],nB$b[i]))
	}
if(sum(L.res<=L.m)!=0){
	ind <- which(L.res<=L.m)
	nB[ind,'open'] <- 1
	}
I <- rbind(I,B[B$open==0,])
B <- rbind(nB,B[B$open==1,])
list(I=I,B=B)
}

rplot <- function(R){
#output plotting 
r <- 0.5
plot(R$I[,'a'],R$I[,'b'],col=1,pch='o',
	xlim=range(R$I[,'a'])+c(-r,r),ylim=range(R$I[,'b'])+c(-r,r),
	xlab="parameter a",ylab="parameter b")
points(B0$a,B0$b,pch="+",col=2,cex=2)
points(R$B[R$B$open==1,'a'],R$B[R$B$open==1,'b'],col=4,pch='o')
points(R$B[R$B$open==0,'a'],R$B[R$B$open==0,'b'],col=3,pch='o')
}
####### End Function #######


### data set (run once)###
set.seed(122)
dat <- rnorm(3,3,1)

### initializing (run once)###
#adjust starting point and threshold likelihood
da <- 0.1 #how far parameter "a " will expand in each iteration
db <- 0.1 #how far parameter "b " will expand in each iteration
B0 <- data.frame(a=3.2,b=0.9,step=0,x=0,y=0,nx=0,ny=0,open=0)
	#starting point
I0 <- NULL #internal point at the beginning
L.m <- -7 #threshold likelihood (boundary condition)
R <- updt.IB(I0,B0,L.m)


### manually expanding and check output ###
#repeat running this line and see how it expand
#"red cross"=starting point; "black circle"=internal point
#"green circle"=point on current boundary but still expandable(likelihood 
#	doesn't reach threshold yet)
#"blue circle"=point on current boundary but not expandable any more
R <- updt.IB(R$I,R$B,L.m) ; rplot(R)


### auto-expanding until boundary ###
t0 <- proc.time()
while( sum(R$B$open)<nrow(R$B) ){ #keep expanding until all points hit threshold
	R <- updt.IB(R$I,R$B,L.m) #expanding
	}
cat("Done!\nDone!\nDone!\n")
(time <- (proc.time()-t0)[1]) #running time
rplot(R)




