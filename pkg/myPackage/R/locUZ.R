
####################################
#### Uniformize locations
####################################
unifLoc <- function(loc, scale=1){
    res <- apply( loc, 2, function(x) (x - min(x))/(max(x) - min(x)) )
    res*scale
  }
####################################
#### locations to Distance matrix
####################################
loc2U <- function(loc_){
  .Call( "loc2Ucpp", loc_, PACKAGE = "myPackage" )
}
# R version
loc2U_R <- function(loc){ 
n <- nrow(loc)
dst <- function(t){sqrt((t[1]-t[3])^2+(t[2]-t[4])^2)}
U <- matrix(0, n, n)
for(i in 1:n){
  U[,i] <- apply(cbind(matrix(rep(loc[i,],each=n),n,),loc),1,dst)
}
U
}
####################################
#### locations to Distance matrix
####################################
locUloc <- function(loc_, locp_){
  .Call( "locUloccpp", loc_, locp_, PACKAGE = "myPackage" )
}
# R version
locUloc_R <- function(loc, locp){
  n <- nrow(loc); np <- nrow(locp)
  dst <- function(t){sqrt((t[1]-t[3])^2+(t[2]-t[4])^2)}
  res <- matrix(0, np,n )
  for(i in 1:n){    
    res[,i] <- apply(cbind(matrix(rep(loc[i,], each=np),np,), locp),1,dst)
  }
  res
}
####################################
#### Correlation functions
####################################
rhoPowerExp <- function(u, a, k) { exp(-(u/a)^k ) }
rhoMatern <- function(u, a, k){
  phi <- a; kappa <- k
  res <- ifelse(u > 0,   ( 2^(1-kappa) / gamma(kappa) ) * (u/phi)^kappa * 
besselK(x = u/phi, nu=kappa, expon.scaled = FALSE) , 1)
  res
}
#########################################
#### Distance matrix to Covariance matrix
#########################################
U2Z <- function(U, cov.par, rho.family = "rhoPowerExp"){
  s <- cov.par[1]; a <- cov.par[2]; k <- cov.par[3]
  if(rho.family=="rhoPowerExp"){
      Z <- s^2* rhoPowerExp(U, a, k)
    } else if(rho.family=="rhoMatern"){
        Z <- s^2* rhoMatern(U, a, k)
      } else {
          cat("Notice: rho.family=", rho.family, " doesn't exist! rho.family=rhoPowerExp will be used.\n", sep="")
          Z <- s^2* rhoPowerExp(U, a, k)
        }
  Z
}
####################################
#### END
####################################