\name{predY}
\alias{predY}

\title{Predict for Given Locations}

\usage{
predY(res.m, loc, locp, X = NULL, Xp = NULL, Lp = 0, k = 1, 
      rho.family = "rhoPowerExp", Y.family = "Poisson")
}

\arguments{
  \item{res.m}{the posterior samples of latent variables and parameters at observed locations.}
  \item{loc}{the observed locations.}
  \item{locp}{the predicted locations.}
  \item{X}{the covariate matrix for observed locations.}
  \item{Xp}{the covariate matrix for predicted locations.}
  \item{Lp}{the time duration for predicted locations.}
  \item{k}{\eqn{\kappa} for correlation function.}
  \item{rho.family}{\code{"rhoPowerExp"} indicates powered exponential correlation function.}
  \item{Y.family}{the distribution for response variables.}
}

\description{
  This function generates posterior predictive samples of latent and reponse variables for predicted locations.}

\details{
  ...
}

\value{
 A list of posterior predictive samples.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
locp <- matrix(runif(200),,2)
Ypred <- predY(res.m, loc, locp, X=loc, Xp=locp, k=1, rho.family="rhoPowerExp")
Ypred.avg <- rowMeans(Ypred$Y); EYpred.avg <- rowMeans(exp(Ypred$Sp))
}
}

\seealso{
   \code{\link{runMCMC}}, \code{\link{runMCMC.multiChain}}.
}

\keyword{prediction}
