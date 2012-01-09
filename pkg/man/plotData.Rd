\name{plotData}
\alias{plotData}

\title{Plot Geostatistical Data}

\usage{
plotData(Y, loc, Yp = NULL, locp = NULL, Yt = NULL, loct = NULL, 
         col = 1:2, colt = 3, pch = 1, size = c(0.3, 2.7), ...)
}

\arguments{
  \item{Y, Yp, Yt}{the vector of response variables}
  \item{loc, locp, loct}{\eqn{n \times 2} matrix that indicates the coordinates of locations}
  \item{col, colt}{the colors used for different sets of response variables}
  \item{pch}{the shape}
  \item{size}{the minimum and maximum of the sizes}
  \item{...}{other parameters that control the plotting}
}

\description{
  This function plots geostatistical data for up to three data sets.
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
  loc <- rbind(locCircle(1, 60), 
               locCircle(0.667, 40), 
               locCircle(0.333, 20)  
               )
  dat <- simData(loc, cov.par = c(1, 0.1, 1))
  Y <- dat$data
  plotData(Y[1:60], loc[1:60, ], Y[61:100], loc[61:100, ], 
           Y[101:120], loc[101:120, ], pch = 16
           )
}
}

\seealso{
   \code{\link{locCircle}},
   \code{\link{locGrid}},
   \code{\link{locSquad}}, \code{\link{plotDataBD}}, \code{\link{simData}}.
}

\keyword{Data}
