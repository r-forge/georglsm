\name{cdfU}
\alias{cdfU}

\title{CDF Value for Observed Data}

\usage{
cdfU(Y.obs, Y.rep, discrete = TRUE)
}

\arguments{
  \item{Y.obs}{the observed data set.}
  \item{Y.rep}{the replicated data sets.}
  \item{discrete}{if the random variable is discrete.}
}

\description{
  This function approximates CDF value for the observed data by using replicated data sets.}

\details{
  ...
} 
\value{
 A vector of CDF values.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Yrep <- repYeb(N.sim=2000, loc, L, beta = 5, sigma = 1, phi = 0.1)
cdfU(Y.obs, Y.rep, discrete = TRUE) }
}

\seealso{
\code{\link{tranR}}
}

\keyword{Model Checking}
