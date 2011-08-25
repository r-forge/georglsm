\name{simData}
\alias{simData}
\title{Simulate Data Set from Generalized Linear Spatial Model on Given Locations}

\usage{
simData(loc, L = 0, X = NULL, beta = 0, cov.par, 
        rho.family = "rhoPowerExp", Y.family = "Poisson")
}

\arguments{
  \item{loc}{ \eqn{n \times 2} matrix that indicates the coordinates of given locations.}
  \item{L}{a vector of length n that indicates the time duration for each location.}
  \item{X}{a n*p covariate matrix.}
  \item{beta}{a vector of length (p+1) that indicates the coefficients}
  \item{cov.par}{a vector of length 3 that indicates \eqn{(\sigma,\phi,\kappa)} }
  \item{rho.family}{\code{"rhoPowerExp"} indicates powered exponential correlation function.}
  \item{Y.family}{\code{"Poisson"} indicates Poisson distribution for response variables}
}

\description{
  This function simulates a data set from Generalized Linear Spatial Model on given locations.
}

\details{
  ...
}

\value{
  A list with two elements:
  \item{data}{a vector indicates the response variables.}
  \item{latent.variables}{a vector indicates the latent variables}
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
   \code{\link{locCircle}},
   \code{\link{locSquad}},
   \code{\link{simData}},
   \code{\link{plotData}}.
}

\keyword{spatial}
