\name{AlcoholTobacco}
\alias{AlcoholTobacco}

\title{Alcohol and Tobacco Expenditures in Belgian Households}

\description{
Data on budget shares of alcohol and tobacco for 2724 Belgian households, 
taken from the Belgian household budget survey of 1995/96. 
}

\usage{data("AlcoholTobacco")}

\format{
A data frame containing 2724 observations on 9 variables.
\describe{
  \item{alcohol}{Budget share of alcohol.}
  \item{tobacco}{Budget share of tobacco.}
  \item{expenditure}{Total expenditure.}
  \item{age}{10-year age intervals ranging from 0 (younger than 30) to 4 (60 or older).}
  \item{adults}{Number of adults.}
  \item{oldkids}{Number of old kids (2 years or older).}
  \item{youngkids}{Number of young kinds (younger than 2 years).}
  \item{occupation}{Factor indicating occupation of household head
    (\code{"blue"} collar worker, \code{"white"} collar worker, or \code{"other"}).}
  \item{region}{Factor indicating region of residence (Dutch-speaking \code{"flanders"},
    French-speaking \code{"wallonie"}, or bilingual \code{"brussels"}).}
}
}

\source{
Supplementary material for Verbeek (2004), originally provided
by the National Institute of Statistics (NIS), Belgium.
}

\references{
Verbeek, M. (2004). \emph{A Guide to Modern Econometrics}, 2nd ed. Chichester, UK: John Wiley.
}

\examples{
data("AlcoholTobacco", package = "htobit")

## homoscedastic tobit model for budget share of alcohol
ma <- htobit(alcohol ~ (age + adults) * log(expenditure) + oldkids + youngkids,
  data = AlcoholTobacco)
summary(ma)

## homoscedastic tobit model for budget share of tobacco
mt <- htobit(tobacco ~ (age + adults) * log(expenditure) + oldkids + youngkids,
  data = AlcoholTobacco)
summary(mt)
}

\keyword{datasets}
