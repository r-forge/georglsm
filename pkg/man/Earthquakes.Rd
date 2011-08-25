\name{Earthquakes}
\alias{Earthquakes}

\title{Data Set of Earthquakes}

\usage{
data(datEarthquake)
}

\arguments{

}

\description{
  This data set contains the informations of earthquakes.}

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

data(datEarthquake)
str(Earthquakes)
plotData(Earthquakes$Magnitude, Earthquakes[,c("Lat","Lon")])
}
}

\seealso{
\code{\link{plotData}}.
}

\keyword{Data}
