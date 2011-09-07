
###################
rm(list=ls(all=T))
###################


t0 <- proc.time()

set.seed(122)
dat <- rnorm(10,3,1)
range <- c(-10,10)+mean(dat)

LH <- function(mu){sum(log(dnorm(dat,mu,1)))}

Rprior <- function(n,range){runif(n,range[1],range[2])}


### Initializing
################

X <- NULL
w <- NULL
L.m <- NULL
L.res <- NULL
l.range <- NULL

N <- 50
j <- 300

lamda.s <- matrix(Rprior(N,range),,1)
L.res <- apply(lamda.s,1,LH)
L.m[1] <- min(L.res)
ind <- which(L.res!=L.m[1])

X[1] <- exp(-1/N)
w[1] <- 1-X[1]

### Main Loop
#############
for(i in 2:j){

l.cur <- lamda.s[sample(ind,1)]

d <- 0.005
direction <- 1
while(direction!= 0){
l.next <- l.cur+d
if(direction==1&&LH(l.next)<=L.m[i-1]){
	l.max <- l.next
	d <- -d
	direction <- -1
	l.next <- l.cur+d
	}
if(direction==-1&&LH(l.next)<=L.m[i-1]){
	l.min <- l.next
	direction <- 0
	}
l.cur <- l.next
}
l.range <- rbind(l.range,c(l.min,l.max))
l.cur <- runif(1,l.min,l.max)

ind.l <- which(L.res==L.m[i-1])
L.res[ind.l] <- LH(l.cur)
lamda.s[ind.l] <- l.cur

L.m[i] <- min(L.res)
ind <- which(L.res!=L.m[i])
X[i] <- exp(-i/N)
w[i] <- X[i-1]-X[i]
}
time <- (proc.time()-t0)[1]

### Result
##########
Z <- sum(exp(L.m)*w)+sum(exp(L.res))*X[j]/N
cat("  CPU time=",time,"; Estimate of Z=",Z,"\n")
plot(c(1,X), exp(c(L.m,mean(L.res))), xlab="X.i", ylab="L.i", 'l',col=3)
matplot(l.range,type="l", xlab="Step", ylab="L.range")
#plot(l.range[,2]-l.range[,1],type='l')








