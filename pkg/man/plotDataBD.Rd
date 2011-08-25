\name{plotDataBD}
\alias{plotDataBD}

\title{Plot Geostatistical Data}

\usage{
plotDataBD(bdry, Y = NULL, loc = NULL, 
        Yp = NULL, locp = NULL, Yt = NULL, loct = NULL, 
        col = 1:2, colt = 3, pch = 1, size = c(0.3, 2.7))
}

\arguments{
  \item{bdry}{a list of boundaries.}
  \item{Y, Yp, Yt}{the response variables at different locations.}
  \item{loc, locp, loct}{\eqn{n \times 2} matrix that indicates the coordinates of locations.}
  \item{col, colt}{the colors for three types of response variables.}
  \item{pch}{the shape.}
  \item{size}{the range of the size.}
}

\description{
  This function plots geostatistical data for up to three data sets and the given boundaries.
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
plotDataBD(TexasCounty.boundary)
}
}

\seealso{
   \code{\link{simData}}, \code{\link{plotData}}.
}

\keyword{Data}
