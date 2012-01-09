\name{baseline.parallel}
\alias{baseline.parallel}

\title{Generate Distance Samples to Build Baseline Distribution (Parallel Version)}

\usage{
baseline.parallel(n, iter, n.cores = getOption("cores"))
}

\arguments{
  \item{n}{the number of residuals}
  \item{iter}{the number of distance samples to generate}
  \item{n.cores}{the number of CPUs that will be used for parallel computing}
}

\description{
  This function generates distance samples in parallel to build the baseline distribution for standard normal.}

\details{
  \code{HellingerDist} and \code{KolmogorovDist} functions in \code{{distrEx}} are used to compute the distances. See \code{?HellingerDist} and \code{?KolmogorovDist} for details about how the distances are computed.
  
  This function performs parallel computing with the help of \code{{multicore}} package. Be aware that \code{{multicore}} package currently is not available in Windows.
} 
\value{
 A \eqn{iter \times 3} matrix for three types of distance: "Discrete Hellinger", "Smooth Hellinger" and "Kolmogorov". 
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
# Time-consuming! Run once with large "iter" and 
# save the results for future use
require(multicore)
d.base <- baseline.parallel(50, iter=100, n.cores = 4) }
}

\seealso{
\code{\link{d.base}}, \code{\link{baseline.dist}}, \code{\link{plot_baseline}}, \code{\link{pOne}}.
}

\keyword{Transformed Residual Checking}
