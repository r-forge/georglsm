\name{plotData}
\alias{plotData}

\title{Plot Geostatistical Data}

\usage{
plotData(Y, loc, Yp = NULL, locp = NULL, Yt = NULL, loct = NULL, 
        col = 1:2, colt = 3, pch = 1, size = c(0.3, 2.7))
}

\arguments{
  \item{Y, Yp, Yt}{the response variables at different locations.}
  \item{loc, locp, loct}{\eqn{n \times 2} matrix that indicates the coordinates of locations.}
  \item{col, colt}{the colors for three types of response variables.}
  \item{pch}{the shape.}
  \item{size}{the range of the size.}
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
  loc <- rbind(locCircle(0.3, 10),
               locCircle(0.6, 30),
               locCircle(1.0, 50)
               )
  dat <- simData(loc, cov.par = c(1, 0.1, 1))
  plotData(dat$data, loc)
}
}

\seealso{
   \code{\link{simData}}.
}

\keyword{plot}
