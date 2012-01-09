
################ Function ##############
plot_dart <- function(x0, divide=FALSE)
  { ### Input.1: x0 = matrix, data.frame, vector;
    ### with nrow(x0) observations and ncol(x0) attributes.
    ### Input.2: divide = FALSE/TRUE; determine if observations
    ### should be plotted in seperate plots.
if(is.data.frame(x0)){
  x <- as.matrix(x0)
} else if(is.vector(x0)){
  x <- matrix(x0,1,)
  dimnames(x)[[2]] <- names(x0)
  if(divide) stop("Can't divide for one object!")
} else x <- x0
n <- nrow(x)
N <- ncol(x)
x.name <- dimnames(x)[[2]]
x <- apply(x, 2, function(t){
  tmax <- max(t); tmin <- min(t)
  0.9*(t-tmin)/(tmax-tmin)+0.1
})
agl <- (0:(N-1))*(2*pi/N)
agl.base <- seq(0,2*pi,len=500)
loc <- loc.n <- vector(mode="list", length=n)
for(i in 1:n)
  {
    tmp <- cbind( x[i,]*cos(agl), x[i,]*sin(agl) )
    loc[[i]] <- tmp
    loc.n[[i]] <- rbind(tmp[-1,], tmp[1,])
 }
if(!divide)
  {
plot(cos(agl.base), sin(agl.base), asp=1,type="l",bty="n", axes=FALSE, ann=FALSE, lty=3)
if(n>1){
  title(main=paste("Comparison of ",n," Objects",sep=""))
if(n<11)
  legend("bottomright",paste("Object.",1:n,sep=""),col=1:n,lty=1, bty="n", cex=ifelse(n<6, 1, 0.8) )
}
text(1.05*cos(agl),1.05*sin(agl),x.name)
lines(0.5*cos(agl.base), 0.5*sin(agl.base), lty=3)
segments(0,0,cos(agl),sin(agl),lty=3 )
for(i in 1:n) segments(loc[[i]][,1],loc[[i]][,2], loc.n[[i]][,1],loc.n[[i]][,2], col=i, lwd=5)
} else{
op <- par(mfrow=c(ceiling(n/2), 2))
for(i in 1:n){
plot(cos(agl.base), sin(agl.base), asp=1,type="l",bty="n", axes=FALSE, ann=FALSE, lty=3)
title(main=paste("Object.",i,sep=""))
text(1.05*cos(agl),1.05*sin(agl),x.name)
lines(0.5*cos(agl.base), 0.5*sin(agl.base), lty=3)
segments(0,0,cos(agl),sin(agl),lty=3 )
segments(loc[[i]][,1],loc[[i]][,2], loc.n[[i]][,1],loc.n[[i]][,2], col=i, lwd=5)
}
par(op)
}
}
############# End of Function ##############


