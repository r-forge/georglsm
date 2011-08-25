\name{BMCT}
\alias{BMCT}

\title{Bayesian Model Checking}

\usage{
BMCT(Y.obs, Y.rep, funcT, ifplot = FALSE)
}

\arguments{
  \item{Y.obs}{the observed data set.}
  \item{Y.rep}{the replicated data sets.}
  \item{funcT}{the dignostic statistic.}
  \item{ifplot}{whether plot the distribution of dignostic statistic.}
}

\description{
  This function conducts Bayesian model checking by comparing observed and replicated data sets and plots the distribution of dignostic statistic.}

\details{
  ...
} 
\value{
 The checking result: p-value and RPS.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
Yrep.eb <- repYeb(N.sim=2000, loc, L, res.m, est = "mode")
funcT <- function(Y){ max(Y)-min(Y) }
BMCT(Y, Yrep.eb, funcT, ifplot=TRUE)  }
}

\seealso{
\code{\link{repYeb}}, \code{\link{repYpost}}, \code{\link{pRPS}}, \code{\link{plot_pRPS}}
}

\keyword{Model Checking}
