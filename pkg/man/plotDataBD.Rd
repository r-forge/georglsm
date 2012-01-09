\name{plotDataBD}
\alias{plotDataBD}

\title{Plot Geostatistical Data}

\usage{
plotDataBD(bdry, Y = NULL, loc = NULL, 
        Yp = NULL, locp = NULL, Yt = NULL, loct = NULL, 
        col = 1:2, colt = 3, pch = 1, size = c(0.3, 2.7), ...)
}

\arguments{
  \item{bdry}{a list containing the coordinates of boundaries}
  \item{Y, Yp, Yt}{the vector of response variables}
  \item{loc, locp, loct}{\eqn{n \times 2} matrix that indicates the coordinates of locations}
  \item{col, colt}{the colors used for different sets of response variables}
  \item{pch}{the shape}
  \item{size}{the minimum and maximum of the sizes}
  \item{...}{other parameters that control the plotting}
}

\description{
  This function plots the given boundaries and geostatistical data sets.
}

\details{
  
}

\value{
 
}

\author{
  Liang Jing \email{ljing918@gmail.com}
  }

\examples{
\dontrun{
data(TexasCounty_boundary)
plotDataBD(TexasCounty.boundary, xlab = "Longitude", ylab = "Latitude")
}
}

\seealso{
   \code{\link{simData}}, \code{\link{plotData}}.
}

\keyword{Data}
