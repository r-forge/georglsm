
####################################
#### Plot data
####################################
plotData <- function(Y, loc, Yp=NULL, locp=NULL, Yt=NULL, loct=NULL, col=1:2, colt=3, pch=1, size=c(0.3, 2.7)){
  
plot( c(loc[,1],locp[,1]), c(loc[,2],locp[,2]), type="n", xlab="", ylab="")
points(loc[,1], loc[,2], col=col[1], pch=pch,
    cex = size[1]+size[2]*(Y-min(Y))/(max(Y)-min(Y)) )
points(locp[,1], locp[,2], col=col[2], pch=pch,
    cex = size[1]+size[2]*(Yp-min(Y))/(max(Y)-min(Y)) )
points(loct[,1], loct[,2], col=colt, pch=pch,
    cex = size[1]+size[2]*(Yt-min(Y))/(max(Y)-min(Y)) )

}
####################################
#### Plot boundary
####################################
plotDataBD <- function(bdry, Y=NULL, loc=NULL, Yp=NULL, locp=NULL, Yt=NULL, loct=NULL, col=1:2, colt=3, pch=1, size=c(0.3, 2.7)){  

n <- length(bdry)
tt <- sapply(bdry, function(t) apply(t, 2, range))
tt1 <- apply(tt, 1, range)
xrange <- range(tt1[1:4])
yrange <- range(tt1[5:8])

plot(xrange, yrange, type="n", xlab = "Longitude", ylab = "Latitude")
for(i in 1:n) lines(bdry[[i]])
if(is.null(Y)) Y <- 0
points(loc[,1], loc[,2], col=col[1], pch=pch,
    cex = size[1]+size[2]*(Y-min(Y))/(max(Y)-min(Y)) )
points(locp[,1], locp[,2], col=col[2], pch=pch,
    cex = size[1]+size[2]*(Yp-min(Y))/(max(Y)-min(Y)) )
points(loct[,1], loct[,2], col=colt, pch=pch,
    cex = size[1]+size[2]*(Yt-min(Y))/(max(Y)-min(Y)) )
}
####################################
#### Plot Texas
####################################
plotTexas <- function(bdry=TexasCounty.boundary, ind.col=NULL){  
n <- length(bdry)
tt <- sapply(bdry, function(t) apply(t, 2, range))
tt1 <- apply(tt, 1, range)
xrange <- range(tt1[1:4])
yrange <- range(tt1[5:8])
plot(xrange, yrange, type="n", xlab = "Longitude", ylab = "Latitude")
for(i in 1:n) polygon(bdry[[i]], col=ind.col[i])
}
####################################
#### END
####################################
