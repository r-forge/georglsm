\name{plot_etran}
\alias{plot_etran}

\title{Calculate Transformed Residuals for Observed Data}

\usage{
plot_etran(e.tran, fig = 1:4)
}

\arguments{
  \item{e.tran}{the transformed residuals.}
  \item{fig}{indicates which plot(s) to plot.}
}

\description{
  This function plots transformed residuals in different type of plots.}

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
Yrep <- repYeb(N.sim=2000, loc, L, beta = 5, sigma = 1, phi = 0.1)
etran <- tranR(Y.obs, Y.rep, discrete = TRUE)
plot_etran(etran)   }
}

\seealso{
\code{\link{tranR}}
}

\keyword{Model Checking}
