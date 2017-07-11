\name{coef.fgamma}

\alias{coef.fgamma}
\alias{vcov.fgamma}
\alias{terms.fgamma}
\alias{model.matrix.fgamma}

\title{Extractor Methods for fgamma Objects}

\description{
  Methods for extracting information from fitted \code{fgamma} objects.
}

\usage{
\S3method{coef}{fgamma}(object, model = c("full", "mu", "sigma"), \dots)
\S3method{vcov}{fgamma}(object, model = c("full", "mu", "sigma"), \dots)

\S3method{terms}{fgamma}(x, model = c("mu", "sigma", "full"), \dots)
\S3method{model.matrix}{fgamma}(object, model = c("mu", "sigma"), \dots)
}
\arguments{
  \item{object, x}{an object of class \code{"fgamma"}.}
  \item{model}{character indicating (sub)model for which information
    should be extracted.}
  \item{\dots}{currently not used.}
}

\details{
  In addition to the methods above, a set of standard extractor functions for 
  \code{"fgamma"} objects is available, see \code{\link{fgamma}} for an overview.
}

\seealso{\code{\link{fgamma}}}

\examples{
## full gamma regression model for averge hourly earnings
data("Wages", package = "fgamma")
f1 <- fgamma(wage ~ educ + exper + tenure  + sex + marriage  |
               tenure, data = Wages)
## extract coefficients
coef(f1)
coef(f1, model = "mu")
coef(f1, model = "sigma")

## corresponding model matrices
head(model.matrix(f1, model = "mu"))
head(model.matrix(f1, model = "sigma"))
}

\keyword{regression}
