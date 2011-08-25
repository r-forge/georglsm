\name{TexasCounty.center}
\alias{TexasCounty.center}

\title{Data Set of Texas County Centers}

\usage{
data(TexasCounty_center)
}

\arguments{

}

\description{
  This data set contains the center locations for all Texas countries.}

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

data(TexasCounty_center)
str(TexasCounty.center)
plotDataBD(TexasCounty.boundary)
points(TexasCounty.center[,2:3], col=2, pch=3)
}
}

\seealso{
\code{\link{TexasCounty.boundary}}, \code{\link{TexasCounty.population}}.
}

\keyword{Data}
