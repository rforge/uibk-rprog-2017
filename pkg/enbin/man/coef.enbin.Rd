\name{coef.enbin}

\alias{coef.enbin}
\alias{vcov.enbin}
\alias{terms.enbin}
\alias{model.matrix.enbin}

\title{Extractor Methods for enbin Objects}

\description{
  Methods for extracting information from fitted \code{enbin} objects.
}

\usage{
\S3method{coef}{enbin}(object, model = c("full", "location", "scale"), \dots)
\S3method{vcov}{enbin}(object, model = c("full", "location", "scale"), \dots)

\S3method{terms}{enbin}(x, model = c("location", "scale", "full"), \dots)
\S3method{model.matrix}{enbin}(object, model = c("location", "scale"), \dots)
}
\arguments{
  \item{object, x}{an object of class \code{"enbin"}.}
  \item{model}{character indicating (sub)model for which information
    should be extracted.}
  \item{\dots}{currently not used.}
}

\details{
  In addition to the methods above, a set of standard extractor functions for 
  \code{"enbin"} objects is available, see \code{\link{enbin}} for an overview.
}

\seealso{\code{\link{enbin}}}

\examples{
## extended negative binomial model
data("RecreationDemand", package = "AER")
m <- enbin(trips ~ . - income | . - income, data = RecreationDemand)

## extract coefficients
coef(m)
coef(m, model = "location")
coef(m, model = "scale")

## corresponding model matrices
head(model.matrix(m, model = "location"))
head(model.matrix(m, model = "scale"))
}

\keyword{regression}
