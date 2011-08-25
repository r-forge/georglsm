\name{runMCMC.multiChain}
\alias{runMCMC.multiChain}

\title{Generate Posterior Samples for GLSM with MCMC Algorithm}

\usage{
runMCMC.multiChain(Y, L = 0, loc, X = NULL, run = 200, run.S = 1, 
        rho.family = "rhoPowerExp", Y.family = "Poisson", ifkappa = 0, 
        scales = c(0.5, 1.65^2 + 0.8, 0.8, 0.7, 0.15), 
        phi.bound = c(0.005, 1), initials = list(c(1), 1.5, 0.2, 1), 
        MCMCinput = NULL, partial = FALSE, famT = 1, n.chn = 2)
}

\arguments{
  \item{Y}{the response variables at given locations.}
  \item{L}{the time duration for each location.}
  \item{loc}{\eqn{n \times 2} matrix that indicates the coordinates of locations.}
  \item{X}{a \eqn{n \times p} covariate matrix.}
  \item{run}{the number of iterations.}
  \item{run.S}{the number of internal iterations for latern variables.}
  \item{rho.family}{\code{"rhoPowerExp"} indicates powered exponential correlation function.}
  \item{Y.family}{\code{"Poisson"} indicates Poisson distribution for response variables}
  \item{ifkappa}{indicates whether \eqn{\kappa} should be sampled.}
  \item{scales}{the tuning parameters for \eqn{(S, \beta, \sigma,\phi,\kappa)}. }
  \item{phi.bound}{the upper and lower bound for \eqn{\phi}. }
  \item{initials}{the initials values.}
  \item{MCMCinput}{the list of alternative settings.}
  \item{partial}{indicats whether partial posterior sampling should be used.}
  \item{famT}{indicate the type of posterior sampleing.}
  \item{n.chn}{the number of running chains.}
}

\description{
  This function generates multiple Markov chains of posterior samples for GLSM.
}

\details{
  ...
}

\value{
 A list containing the posterior samples for \eqn{(S, \beta, \sigma,\phi,\kappa)} and the acceptance rates.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  require(multicore)
  options(cores=5)
  input <- MCMCinput( run = 200, run.S = 1,
    rho.family = "rhoPowerExp", Y.family = "Poisson", ifkappa = 0,
    scales = c(0.5, 1.65^2+0.8, 0.8, 0.7, 0.15), 
    phi.bound = c(0.005, 1), 
    initials = list(c(1), 1.5, 0.2, 1) )
  res.prl <- runMCMC.multiChain(Y, L=0, loc=loc, X=loc, MCMCinput = input, n.chn = 5 )
}
}

\seealso{
   \code{\link{MCMCinput}}, \code{\link{runMCMC}}.
}

\keyword{MCMC}
