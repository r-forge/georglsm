\name{baseline.dist}
\alias{baseline.dist}

\title{Generate Baseline Samples for Standard Normal}

\usage{
baseline.dist(n, iter)
}

\arguments{
  \item{n}{the number of residuals.}
  \item{iter}{the number of baseline samples to generate.}
}

\description{
  This function generates baseline samples for standard normal distribution.}

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
# (time-consuming! only need to run once for date sets with same size)
d.base <- baseline.dist(length(etran), iter=1000) }
}

\seealso{
\code{\link{tranR}}, \code{\link{e2dist}}
}

\keyword{Model Checking}
