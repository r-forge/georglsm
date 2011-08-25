\name{repYpost}
\alias{repYpost}

\title{Generation of Replicated Data Sets}

\usage{
repYpost(res.m, L)
}

\arguments{
  \item{res.m}{the posterior samples of latern variables to be used for generation.}
  \item{L}{the time duration for all the locations.}
}

\description{
  This function generats replicated data sets based on posterior samples of latern variables.}

\details{
  ...
} 
\value{
 A matrix of replicated data sets.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
# Estimation from posterior samples
Yrep.post <- repYpost(res.m, L) }
}

\seealso{
\code{\link{repYeb}}
}

\keyword{Replicated Data}
