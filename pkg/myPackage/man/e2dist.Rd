\name{e2dist}
\alias{e2dist}

\title{Calculate Distances between Transformed Residuals and Standard Normal}

\usage{
e2dist(e.tran)
}

\arguments{
  \item{e.tran}{the transformed residuals.}
}

\description{
  This function calcualtes different distances between the empirical distribution of transformed residuals and standard normal.}

\details{
  ...
} 
\value{
 A matrix of distances.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Yrep <- repYeb(N.sim=2000, loc, L, beta = 5, sigma = 1, phi = 0.1)
etran <- tranR(Y.obs, Y.rep, discrete = TRUE)
require(distrEx)
e2dist(etran)       }
}

\seealso{
\code{\link{tranR}}, \code{\link{plot_etran}}
}

\keyword{Model Checking}
