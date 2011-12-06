

####################################
#### runMCMC_: internal
#### Call runMCMCBPcpp() in C++
####################################
runMCMC_ <- function(Y_, L_, T_, D_, run_, nmLan_, fam_, famY_, ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, alow_, aup_, mini_, sini_, aini_, kini_){
  .Call( "runMCMCBPcpp", Y_, L_, T_, D_, run_, nmLan_, fam_, famY_, ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, alow_, aup_, mini_, sini_, aini_, kini_, PACKAGE = "myPackage" )
}
runMCMCpartialPois_ <- function(Y_, L_, T_, D_, run_, nmLan_, fam_, famY_, famT_, ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, alow_, aup_, mini_, sini_, aini_, kini_){
  .Call( "runMCMCpartialPoiscpp", Y_, L_, T_, D_, run_, nmLan_, fam_, famY_, famT_, ifkappa_, scale_, mscale_, sscale_, ascale_, kscale_, alow_, aup_, mini_, sini_, aini_, kini_, PACKAGE = "myPackage" )
}
####################################
#### set up MCMCinput
####################################
MCMCinput <- function( run = 200, run.S = 1,
    rho.family = "rhoPowerExp", Y.family = "Poisson", ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1) ){
    list( run=run, run.S=run.S, rho.family=rho.family, 
          Y.family=Y.family, ifkappa=ifkappa,
          scales=scales, phi.bound=phi.bound, initials=initials )  
    }
####################################
#### runMCMC(): base 
####################################
runMCMC <- function( Y, L=0, loc, X=NULL, 
    run = 200, run.S = 1,
    rho.family = "rhoPowerExp", Y.family = "Poisson", ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1), 
    MCMCinput=NULL, partial = FALSE, famT=1 ){
      
if(!is.null(MCMCinput)){
  run <- MCMCinput$run; run.S <- MCMCinput$run.S
  rho.family <- MCMCinput$rho.family; Y.family <- MCMCinput$Y.family;
  ifkappa <- MCMCinput$ifkappa
  scales <- MCMCinput$scales; 
  phi.bound <- MCMCinput$phi.bound
  initials <- MCMCinput$initials
  }
  
  if(any(Y<=0)){
      warning("\nY contains non-positive element! \n0.1 is added to all elements!")
      Y <- Y+0.1
    } 
  Y <- matrix(Y,,1)
  if(any(L==0)){
    L <- matrix(rep(1,nrow(Y)),,1)
    warning("\nL contains zero!\nL is set to 1 for all locations!")
    } else { L <- matrix(L,,1)}
  U <- loc2U(loc)
  D <- cbind( rep(1, nrow(Y)), X )
  
  if(rho.family=="rhoPowerExp"){
      fam = 1
    } else if(rho.family=="rhoMatern"){
        fam = 2
      } else {
          warning( paste("\nrho.family=", rho.family, " doesn't exist! \nrho.family=rhoPowerExp is used!", sep="") )
          fam = 1
        }
  if(Y.family=="Poisson"){
      famY = 1
    } else if(Y.family=="Binomial"){
        famY = 2
      } else {
          warning( paste("\nY.family=", Y.family, " doesn't exist! \nY.family=Poisson is used!", sep="") )
          famY = 1
        }
        
  scale <- scales[1]; mscale <- scales[2]; sscale <- scales[3]; 
  ascale <- scales[4]; kscale <- scales[5]
  alow <- phi.bound[1]; aup <- phi.bound[2]
  mini <- matrix(initials[[1]],,1); sini <- initials[[2]]; 
  aini <- initials[[3]]; kini <- initials[[4]];
  if(ncol(D) != nrow(mini)) stop("The number of covariates is not equal to the number of coefficients!")
  
  message("### MCMC Starts!\n")
  t0 <- proc.time()
  if(partial){
    res <- runMCMCpartialPois_(Y, L, U, D, run, run.S, fam, famY, famT, ifkappa, scale, mscale, sscale, ascale, kscale, alow, aup, mini, sini, aini, kini)
  } else {
    res <- runMCMC_(Y, L, U, D, run, run.S, fam, famY, ifkappa, scale, mscale, sscale, ascale, kscale, alow, aup, mini, sini, aini, kini)
  }
  run.time <- proc.time() - t0
  message("### MCMC Done!\n")
  
  message("### MCMC Running Time: ")
  print(run.time)
  message("### MCMC Acceptance Rate: ")
  print(res$Acc)

  res[[1]] <- matrix(res[[1]], nrow(Y), , )
  if(nrow(mini)!=1) res[[2]] <- matrix(res[[2]], nrow(mini), , )
  return(res)
  }
####################################
#### runMCMC.multiChain(): multiple
# require pre-load library{multicore}
####################################
runMCMC.multiChain <- function(Y, L=0, loc=loc, X=NULL, 
    run = 200, run.S = 1,
    rho.family = "rhoPowerExp", Y.family = "Poisson", ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1), 
    MCMCinput=NULL, partial = FALSE, famT=1,
    n.chn = 2, n.cores = getOption("cores")) {
      
if(!is.null(MCMCinput)){
  run <- MCMCinput$run; run.S <- MCMCinput$run.S
  rho.family <- MCMCinput$rho.family; Y.family <- MCMCinput$Y.family;
  ifkappa <- MCMCinput$ifkappa
  scales <- MCMCinput$scales; 
  phi.bound <- MCMCinput$phi.bound
  initials <- MCMCinput$initials
  }
  
## set different starting points    
s.ini <- 0.1*c( mean(initials[[1]]), initials[[2]], initials[[3]], initials[[4]])
ini <- lapply( 1:n.chn, function(tt)
          list( initials[[1]]+s.ini[1]*rnorm(length(initials[[1]])), 
            abs(initials[[2]]+s.ini[2]*rnorm(1)), 
            abs(initials[[3]]+s.ini[3]*rnorm(1)),
            ifelse(ifkappa==0, initials[[4]],
                  abs(initials[[4]]+s.ini[4]*rnorm(1)) )
            )
          )
          
  message("### MCMC Starts!\n")
  t0 <- proc.time()
  res.prl <- mclapply(1:n.chn, function(t) 
              runMCMC(Y, L, loc, X, run, run.S, rho.family, Y.family, 
                      ifkappa, scales, phi.bound, ini[[t]], MCMCinput=NULL,
                      partial, famT), 
                      mc.cores = n.cores
          )
  run.time <- proc.time() - t0
  message("### MCMC Done!\n")
  
  message("### MCMC Running Time: ")
  print(run.time)
  message("### MCMC Acceptance Rate: ")
  print(sapply(res.prl, function(tt) tt$Acc))
  
  res.prl
}

####################################
#### runMCMC.sf(): multiple-chain
# require pre-load library{snowfall}
####################################
runMCMC.sf <- function(Y, L=0, loc=loc, X=NULL, 
    run = 200, run.S = 1,
    rho.family = "rhoPowerExp", Y.family = "Poisson", ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1), 
    MCMCinput=NULL, partial = FALSE, famT=1,
    n.chn = 2, cluster.type="SOCK", n.cores = getOption("cores")) {
      
if(!is.null(MCMCinput)){
  run <- MCMCinput$run; run.S <- MCMCinput$run.S
  rho.family <- MCMCinput$rho.family; Y.family <- MCMCinput$Y.family;
  ifkappa <- MCMCinput$ifkappa
  scales <- MCMCinput$scales; 
  phi.bound <- MCMCinput$phi.bound
  initials <- MCMCinput$initials
  }
  
## set different starting points    
s.ini <- 0.1*c( mean(initials[[1]]), initials[[2]], initials[[3]], initials[[4]])
ini <- lapply( 1:n.chn, function(tt)
          list( initials[[1]]+s.ini[1]*rnorm(length(initials[[1]])), 
            abs(initials[[2]]+s.ini[2]*rnorm(1)), 
            abs(initials[[3]]+s.ini[3]*rnorm(1)),
            ifelse(ifkappa==0, initials[[4]],
                  abs(initials[[4]]+s.ini[4]*rnorm(1)) )
            )
          )
  
sfInit(parallel=TRUE, cpus= n.cores, type=cluster.type)
sfExportAll( except=NULL, debug=FALSE )
sfLibrary(myPackage)
sfClusterSetupRNG()
  message("### MCMC Starts!\n")
  t0 <- proc.time()
  res.prl <- sfLapply(1:n.chn, function(t) 
              runMCMC(Y, L, loc, X, run, run.S, rho.family, Y.family, 
                      ifkappa, scales, phi.bound, ini[[t]], MCMCinput=NULL,
                      partial, famT)
              )
  run.time <- proc.time() - t0
  message("### MCMC Done!\n")
sfStop()

  message("### MCMC Running Time: ")
  print(run.time)
  message("### MCMC Acceptance Rate: ")
  print(sapply(res.prl, function(tt) tt$Acc))

  res.prl
}

####################################
#### Prediction at new locations
####################################
predY <- function(res.m, loc, locp, X=NULL, Xp=NULL, Lp=0, k=1, rho.family="rhoPowerExp", Y.family="Poisson", mc=FALSE){

t0 <- proc.time()

n <- nrow(loc); np <- nrow(locp); ns <- ncol(res.m$S)
Sp.post <- Yp <- matrix(0,np,ns)

  if(any(Lp==0)){
    Lp <- matrix(rep(1,np),,1)
    warning("\nLp contains zero!\nLp is set to 1 for all locations!")
    } else { Lp <- matrix(Lp,,1)}
    
Uxx <- loc2U(loc)
Uyy <- loc2U(locp)
Uyx <- locUloc(loc, locp)
Dx <- cbind( rep(1, n), X )
Dy <- cbind( rep(1, np), Xp )

if( !is.null(res.m$k) ){
  k.post <- res.m$k
  } else  k.post <- rep(k, ns) 


  message("### Prediction Starts!\n")

if(!mc){
  res.prl <- lapply(1:ns, function(i){
  
x <-res.m$S[,i]; 
s <- res.m$s[i]; a <- res.m$a[i]; k <- k.post[i]
if(is.matrix(res.m$m)){
  m <- res.m$m[,i]
  } else m <- res.m$m[i]

mu.x <- Dx%*%m; mu.y <- Dy%*%m

  if(rho.family=="rhoPowerExp"){
      Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
      Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
      Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
    } else if(rho.family=="rhoMatern"){
        Z.xx <- s^2* rhoMatern(Uxx, a, k)
        Z.yy <- s^2* rhoMatern(Uyy, a, k)
        Z.yx <- s^2* rhoMatern(Uyx, a, k)
      } else {
          cat("Notice: rho.family=", rho.family, " doesn't exist! rho.family=rhoPowerExp will be used.\n", sep="")
                Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
                Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
                Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
        }
E <- mu.y + Z.yx%*%solve(Z.xx)%*%(x-mu.x)
V <- Z.yy - Z.yx%*%solve(Z.xx)%*%t(Z.yx)
z <- rnorm(np)
y <- E + chol(V)%*%z
Sp.post[,i] <- y; 
if(Y.family=="Poisson"){
  Yp[,i] <- rpois(np, Lp*exp(y))
  } else if(Y.family=="Binomial"){
          Yp[,i] <- rbinom(np, Lp, exp(y)/(1+exp(y)))
          }
c(Sp.post[,i], Yp[,i])
} 
  )
} else {
  res.prl <- mclapply(1:ns, function(i){
  
x <-res.m$S[,i]; 
s <- res.m$s[i]; a <- res.m$a[i]; k <- k.post[i]
if(is.matrix(res.m$m)){
  m <- res.m$m[,i]
  } else m <- res.m$m[i]

mu.x <- Dx%*%m; mu.y <- Dy%*%m

  if(rho.family=="rhoPowerExp"){
      Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
      Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
      Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
    } else if(rho.family=="rhoMatern"){
        Z.xx <- s^2* rhoMatern(Uxx, a, k)
        Z.yy <- s^2* rhoMatern(Uyy, a, k)
        Z.yx <- s^2* rhoMatern(Uyx, a, k)
      } else {
          cat("Notice: rho.family=", rho.family, " doesn't exist! rho.family=rhoPowerExp will be used.\n", sep="")
                Z.xx <- s^2* rhoPowerExp(Uxx, a, k)
                Z.yy <- s^2* rhoPowerExp(Uyy, a, k)
                Z.yx <- s^2* rhoPowerExp(Uyx, a, k)
        }
E <- mu.y + Z.yx%*%solve(Z.xx)%*%(x-mu.x)
V <- Z.yy - Z.yx%*%solve(Z.xx)%*%t(Z.yx)
z <- rnorm(np)
y <- E + chol(V)%*%z
Sp.post[,i] <- y; 
if(Y.family=="Poisson"){
  Yp[,i] <- rpois(np, Lp*exp(y))
  } else if(Y.family=="Binomial"){
          Yp[,i] <- rbinom(np, Lp, exp(y)/(1+exp(y)))
          }
c(Sp.post[,i], Yp[,i])
} 
  )
}

run.time <- proc.time() - t0
  message("### Prediction Done!\n")
  message("### Prediction Running Time: ")
  print(run.time)
res <- matrix(unlist(res.prl),,ns)
list(Sp.posterior=res[1:np,], Y.predict=res[(np+1):(2*np),])

}







          