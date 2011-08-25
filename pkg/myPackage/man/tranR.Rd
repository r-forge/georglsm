\name{tranR}
\alias{tranR}

\title{Calculate Transformed Residuals for Observed Data}

\usage{
tranR(Y.obs, Y.rep, discrete = TRUE)
}

\arguments{
  \item{Y.obs}{the observed data set.}
  \item{Y.rep}{the replicated data sets.}
  \item{discrete}{if the random variable is discrete.}
}

\description{
  This function approximates transformed residuals for the observed data by using replicated data sets.}

\details{
  ...
} 
\value{
 A vector of transformed residuals.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Yrep <- repYeb(N.sim=2000, loc, L, beta = 5, sigma = 1, phi = 0.1)
tranR(Y.obs, Y.rep, discrete = TRUE) }
}

\seealso{
\code{\link{cdfU}}, \code{\link{plot_etran}}, \code{\link{e2dist}}
}

\keyword{Model Checking}
