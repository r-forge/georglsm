\name{MCMCinput}
\alias{MCMCinput}

\title{Settings for GLSM MCMC Algorithm}

\usage{
MCMCinput(run = 200, run.S = 1, rho.family = "rhoPowerExp", 
        Y.family = "Poisson", ifkappa = 0, 
        scales = c(0.5, 1.65^2 + 0.8, 0.8, 0.7, 0.15), 
        phi.bound = c(0.005, 1), 
        initials = list(c(1), 1.5, 0.2, 1))
}

\arguments{
  \item{run}{the number of iterations.}
  \item{run.S}{the number of internal iterations for latern variables.}
  \item{rho.family}{\code{"rhoPowerExp"} indicates powered exponential correlation function.}
  \item{Y.family}{\code{"Poisson"} indicates Poisson distribution for response variables}
  \item{ifkappa}{indicates whether \eqn{\kappa} should be sampled.}
  \item{scales}{the tuning parameters for \eqn{(S, \beta, \sigma,\phi,\kappa)}. }
  \item{phi.bound}{the upper and lower bound for \eqn{\phi}. }
  \item{initials}{the initials values.}
}

\description{
  This function sets up the parameters used for GLSM MCMC algorithm.
}

\details{
  ...
}

\value{
 A list of setting parameters.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  input <- MCMCinput( run = 200, run.S = 1,
    rho.family = "rhoPowerExp", Y.family = "Poisson", ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1) )
  res <- runMCMC(Y, L=0, loc=loc, X=loc, MCMCinput = input )
}
}

\seealso{
   \code{\link{runMCMC}}, \code{\link{runMCMC.multiChain}}.
}

\keyword{MCMC}
