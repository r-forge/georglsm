\name{plotACF}
\alias{plotACF}

\title{Auto-correlation Plot for Latent Variables}

\usage{
plotACF(S.mcmc)
}

\arguments{
  \item{S.mcmc}{the posterior samples of latent variables.}
}

\description{
  This function plots auto-correlation for latent variables.}

\details{
  ...
} 
\value{
 No return value. A plot of auto-correlation.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
require(coda)
S.mcmc <- mcmc(t(res.m$S))
plotACF(S.mcmc)
}
}

\seealso{

}

\keyword{ACF}
