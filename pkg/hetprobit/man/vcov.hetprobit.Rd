\name{vcov.hetprobit}
\alias{vcov.hetprobit}

\title{
Calculate Variance-Covariance Matrix for a Heteroscedastic Probit Model Object
}

\description{
Returns the variance-covariance matrix from a fitted heteroscedastic probit model object.
}

\usage{
\method{vcov}{hetprobit}(object, model = c("full", "mean", "scale"), \ldots)
}

\arguments{
  \item{object}{a fitted model object of class \code{\link[hetprobit]{hetprobit}}.}
  \item{model}{for which component of the model the variance-covariance matrix should be returned.} 
  \item{\dots}{further arguments to be passed.}
}

\details{
  The default returns the variance-covariance estimate of the full model, i.e. of the mean and scale 
  coefficient estimates. Alternatively, the variance-covariance estimate of either the \code{mean} or 
  \code{scale} estimates can be extracted. 
}

\value{
  A matrix of variance-covariance estimates.
}

\seealso{
  \code{\link[hetprobit]{hetprobit}}
}
\keyword{regression}
