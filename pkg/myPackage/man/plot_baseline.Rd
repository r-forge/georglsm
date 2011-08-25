\name{plot_baseline}
\alias{plot_baseline}

\title{Plot Baseline Samples}

\usage{
plot_baseline(res.in, dist.name)
}

\arguments{
  \item{res.in}{the baseline samples.}
  \item{dist.name}{the name of distance.}
}

\description{
  This function plots the baseline samples.}

\details{
  ...
} 
\value{

}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
# (time-consuming! only need to run once for date sets with same size)
d.base <- baseline.dist(length(etran), iter=1000) 
plot_baseline(d.base[,1], colnames(d.base)[1])
plot_baseline(d.base[,2], colnames(d.base)[2])
plot_baseline(d.base[,3], colnames(d.base)[3])      }
}

\seealso{
\code{\link{e2dist}}, \code{\link{baseline.dist}}
}

\keyword{Model Checking}
