\name{locSquad}
\alias{locSquad}
\title{Simulate Squared Locations}

\usage{
locSquad(a, np)
}

\arguments{
  \item{a}{half length of the edge}
  \item{np}{the number of locations on each edge}
}

\description{
  This function simulates a given number of locations equally distributed on the square.
}

\details{
  The center of the square is (0, 0).
}

\value{
  A (4*np-4)*2 matrix indicates the x-y coordinates of the locations.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
  plot(locSquad(0.5, 4))
}
}

\seealso{
   \code{\link{locCircle}},
   \code{\link{simData}},
   \code{\link{plotData}}.
}

\keyword{location}
