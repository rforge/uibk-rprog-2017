\name{VoterTurnout}
\alias{VoterTurnout}

\title{
Voter Turnout in U.S. Presidential Elections 1984
}
\description{
Census data from the U.S. Presidential Elections 1984 
}

\usage{data("VoterTurnout")}

\format{
  A data frame with 98857 observations on the following 6 variables.
  \describe{
    \item{vote}{factor. Did the respondent vote?}
    \item{education}{Years of education of the respondent}
    \item{age}{Age of the respondent}
    \item{south}{factor. Is the respondent from the South?}
    \item{govelection}{factor. Whether gubernatorial elections took place}
    \item{closing}{Number of days before the election that the registration closes}
  }
}

\details{
  Nagler (1991) first analyzed the data using the standard probit model. 
  In 1994 he introduced the skewed logit model (a.k.a scobit-model) and 
  Altman and McDonald (2003) replicated this study with focus on numerical accuracy.
}

\source{
  Supplementary material to Altman and McDonald (2003).  
  \doi{10.1093/pan/mpg016}
}

\references{
  Altman M., McDonald M.P. (2003) Replication with Attention to Numerical Accuracy. 
  \emph{Political Analysis}, \bold{11}, 302--307.

  Nagler J. (1991) The Effects of Registration Laws and Education on U.S. Voter Turnout. 
  \emph{The American Political Science Review}, \bold{85}(4), 1393--1405.

  Nagler J. (1994) Scobit: An Alternative Estimator to Logit and Probit. 
  \emph{Political Science}, \bold{38}(1), 230--55.
}

\seealso{
  \code{\link[hetprobit]{hetprobit}}
}


\examples{
data("VoterTurnout")

## homoscedastic probit model
## regressors are the same as in Table 3 in Nagler (1991)
m0 <- glm(vote ~ age + I(age^2) + south + govelection + 
                (education + I(education^2)) * closing, data = VoterTurnout,
                family = binomial(link = "probit"))
summary(m0)

## heteroscedastic probit model 
## main effects in the mean model and one regressor in the scale model 
m1 <- hetprobit(vote ~ education + age + south + govelection + closing | 
                education, data = VoterTurnout) 

summary(m1)
}

\keyword{datasets}
