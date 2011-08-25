\name{repYeb}
\alias{repYeb}

\title{Generation of Replicated Data Sets}

\usage{
repYeb(N.sim, loc, L, X = NULL, rho.family = "rhoPowerExp", 
    res.m = NULL, est = "mode", beta = NULL, sigma = NULL, phi = NULL, 
    k = 1)
}

\arguments{
  \item{N.sim}{the number of replicated data sets to be simulated.}
  \item{loc}{\eqn{n \times 2} matrix that indicates the coordinates of locations.}
  \item{L}{the time duration for all the locations.}
  \item{X}{the covariate matrix.}
  \item{rho.family}{the correlation function to be used for latent variables.}
  \item{res.m}{the posterior samples to be used for parameter estimation.}
  \item{est}{the way to estimate.}
  \item{beta}{alternative estimation for \eqn{\beta}}
  \item{sigma}{alternative estimation for \eqn{\sigma}}
  \item{phi}{alternative estimation for \eqn{\phi}}
  \item{k}{alternative estimation for \eqn{\kappa}}
}

\description{
  This function generats replicated data sets based on estimated parameters.}

\details{
  ...
} 
\value{
 A matrix of replicated data sets.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
# Estimation from posterior samples
Yrep.eb <- repYeb(N.sim=2000, loc, L, res.m, est = "mode")
# Given estimated parameters
Yrep.eb <- repYeb(N.sim=2000, loc, L, beta = 5, sigma = 1, phi = 0.1, 
    k = 1)
}
}

\seealso{
\code{\link{repYpost}}
}

\keyword{Replicated Data}
