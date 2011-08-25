\name{cutChain}
\alias{cutChain}

\title{Modify MCMC Chains with Burn-in and Thining}

\usage{
cutChain(res, chain.ind, burnin, thining)
}

\arguments{
  \item{res}{the posterior samples of latent variables and parameters.}
  \item{chain.ind}{the index of chains.}
  \item{burnin}{the number of samples for burn-in.}
  \item{thining}{the number of samples for thining.}
}

\description{
  This function takes the results from running MCMC algorithms and modifies the chains with burn-in and thining.}

\details{
  ...
} 
\value{
 A list of modified chains.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
res <- runMCMC(Y, L=0, loc=loc, X=loc, MCMCinput = input )
res.m <- cutChain(res, chain.ind=1:4, burnin=100, thining=10)
}
}

\seealso{
   \code{\link{runMCMC}}, \code{\link{runMCMC.multiChain}}.
}

\keyword{MCMC}
