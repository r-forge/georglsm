\name{mixChain}
\alias{mixChain}

\title{Mix Multiple MCMC Chains}

\usage{
mixChain(res.m.prl)
}

\arguments{
  \item{res.m.prl}{multiple chains of the posterior samples.}
}

\description{
  This function mix multiple chains into one chain.}

\details{
  ...
} 
\value{
 A list of mixed chains.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
res.prl <- runMCMC.multiChain(Y, L=0, loc=loc, X=loc, 
    MCMCinput = input.rong, n.chn = 5)
res.m.prl <- lapply(res.prl, cutChain, chain.ind=1:4, burnin=200, thining=20)
res.mix <- mixChain(res.m.prl)
}
}

\seealso{
   \code{\link{runMCMC.multiChain}}, \code{\link{cutChain}}.
}

\keyword{MCMC}
