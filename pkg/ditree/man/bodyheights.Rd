\name{bodyheights}
\alias{bodyheights}

\title{Heights of children and their parents}

\description{
This data set is extracted from the data set 'PearsonLee' which is available in the package \pkg{HistData}.
It contains the heights of children and their parents together with two variables indicating whether the height of the son or daughter and height of the mother or father was meassured.
}

\usage{data("bodyheights")}

\format{
A data frame containing 746 observations on 4 variables.
\describe{
  \item{child}{child height in inches, a numeric vector}
  \item{parent}{parent height in inches, a numeric vector}
  \item{par}{a factor with levels Father Mother}
  \item{chl}{a factor with levels Daughter Son}
}
}

\source{
Pearson, K. and Lee, A. (1896). Mathematical contributions to the theory of evolution. On telegony in man, etc. Proceedings of the Royal Society of London, 60 , 273-283.

Pearson, K. and Lee, A. (1903). On the laws of inheritance in man: I. Inheritance of physical characters. Biometika, 2(4), 357-462. (Tables XXII, p. 415; XXV, p. 417; XXVIII, p. 419 and XXXI, p. 421.)
}

\references{
Wachsmuth, A.W., Wilkinson L., Dallal G.E. (2003). Galton's bend: A previously undiscovered nonlinearity in Galton's family stature regression data. The American Statistician, 57, 190-192. \url{http://www.cs.uic.edu/~wilkinson/Publications/galton.pdf}
}

\examples{
library(ditree)
## load data
data("bodyheights", package = "ditree")

## distributional model for heights of daughters
hd <- bodyheights[levels(bodyheights$chl) == "Daughter", "child"]
m_hd <- difit(hd, family = NO())
print(m_hd)
plot(m_hd)
summary(m_hd)

## distributional model for heights of sons
hs <- bodyheights[levels(bodyheights$chl) == "Son", "child"]
m_hs <- difit(hs, family = NO())
print(m_hs)
plot(m_hs)
summary(m_hs)

## building tree models
th_chl <- ditree(child~chl, data = bodyheights)
plot(th_chl)

th <- ditree(child~parent, data = bodyheights)
plot(th)

th_all <- ditree(child~parent + chl, data = bodyheights, control = mob_control(minsize = 80))
plot(th_all)
}

\keyword{datasets}
