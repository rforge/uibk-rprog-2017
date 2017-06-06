\name{publications}
\alias{publications}

\title{Number of publications produced by Ph.D biochemists}

\description{
Data on the number of publications produced by Ph.D. biochemists, 
taken from Long (1990). 
}

\usage{data("publications")}

\format{
A data frame containing 915 observations on 6 variables.
\describe{
  \item{articles}{integer. articles in the last three years of Ph.D.}
  \item{gender}{factor. coded one for females.}
  \item{married}{factor. coded one if married.}
  \item{kids}{integer. number of children under age six.}
  \item{prestige}{numeric. prestige of Ph.D. program}
  \item{mentor}{integer. articles by mentor in last three years}
}
}

\source{
Supplementary material for Long (1990), in Stata format availlable at
\url{http://www.indiana.edu/~jslsoc/stata/socdata/couart2.dta}.
}

\references{
Long, J. Scott (1990). "The origins of sex differences in science", \emph{Social Forces},  68:1297--1315.
Long, J. Scott (1997). \emph{Regression Models for Categorical and Limited Dependent Variables}. Thousand Oaks, California: Sage.
}

\examples{
data("publications", package = "negbin1")

## Poisson model for the number of publications
poi <- glm(articles ~ gender + married + kids + prestige + mentor,
  family = "poisson", data = publications)
summary(poi)

## NB1 model for the number of publications
nb11 <- negbin1(articles ~ gender + married + kids + prestige + mentor,
   data = publications)
summary(nb11)

}

\keyword{datasets}
