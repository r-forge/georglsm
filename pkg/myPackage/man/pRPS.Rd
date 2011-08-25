\name{pRPS}
\alias{pRPS}

\title{Calculate P-value and RPS}

\usage{
pRPS(T.obs, T.rep)
}

\arguments{
  \item{T.obs}{the observed statistic.}
  \item{T.rep}{the replicated statistics.}
}

\description{
  This function calculates p-value of RPS by comparing observed and replicated statistics.}

\details{
  ...
} 
\value{
 p-value and RPS.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
pRPS(2, rnorm(1000)) }
}

\seealso{
\code{\link{BMCT}}, \code{\link{plot_pRPS}}
}

\keyword{Model Checking}
