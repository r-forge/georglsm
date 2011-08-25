\name{plot_pRPS}
\alias{plot_pRPS}

\title{Plot for Observed and Replicated Statistics}

\usage{
plot_pRPS(T.obs, T.rep, nm = "x")
}

\arguments{
  \item{T.obs}{the observed statistic.}
  \item{T.rep}{the replicated statistics.}
  \item{nm}{the label of statistics.}
}

\description{
  This function plots the observed statistic and the empirical density of replicated statistics.}

\details{
  ...
} 
\value{
 A plot.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
plot_pRPS(2, rnorm(1000)) }
}

\seealso{
\code{\link{BMCT}}, \code{\link{pRPS}}
}

\keyword{Model Checking}
