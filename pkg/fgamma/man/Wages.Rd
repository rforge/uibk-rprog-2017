\name{Wages}
\alias{Wages}

\title{U.S. cross-section wage data}

\description{
  Cross-section wage data consisting of a random sample with 526 observations taken from the 
  U.S. Current Population Survey for the year 1976. 
}

\usage{data("Wages")}

\format{
  A data frame containing 526 observations on 8 variables.
  \describe{
    \item{wage}{Average hourly earnings.}
    \item{educ}{Years of education.}
    \item{exper}{Years potential experience.}
    \item{tenure}{Years with current employer.}
    \item{marriage}{Factor indicating \emph{married} and \emph{nonmarried}.}
    \item{numdep}{Number of dependents.}
    \item{ethnicity}{Factor indicating \emph{nonwhite} and \emph{white}.}
    \item{sex}{Factor indicating \emph{female} and \emph{male}.}
  }
}

\source{
  Jeffrey M. Wooldridge 
}

\references{
  Wooldridge, J.M.  (2000). \emph{Introductory Econometrics: A Modern Approach}, South-Western College Publishing.
}

\examples{
  data("Wages", package = "fgamma")
  
  ## Gamma model with only intercept for sigma parameter
  f1 <- fgamma(wage ~ ., data = Wages)
  summary(f1)
  
  ## Gamma model with all explanatory variables for sigma parameter
    f2 <- fgamma(wage ~ . | ., data = Wages)
  summary(f2)
  
  ## Gamma model with selected explanatory variables for sigma parameter
  f3 <- fgamma(wage ~ . | educ + exper + tenure, data = Wages)
  summary(f3)
}

\keyword{datasets}
