\name{TexasCounty.boundary}
\alias{TexasCounty.boundary}

\title{Data Set of Texas County Boundries}

\usage{
data(TexasCounty_boundary)
}

\arguments{

}

\description{
  This data set contains the boundary information for all Texas countries.}

\details{
  
} 
\value{
  A list with 254 elements each of which contains the boundary information for one county.
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
data(TexasCounty_boundary)
length(TexasCounty.boundary); names(TexasCounty.boundary)
plotDataBD(TexasCounty.boundary, xlab = "Longitude", ylab = "Latitude")
text(TexasCounty.center[,2:3], names(TexasCounty.boundary), cex=0.4)
}
}

\seealso{
\code{\link{plotDataBD}}, \code{\link{TexasCounty.center}}, \code{\link{TexasCounty.population}}.
}

\keyword{Data}
