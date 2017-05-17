\name{negbin1}
\alias{negbin1}
\alias{negbin1_fit}
\alias{negbin1_control}

\title{
Negative Binomial 1 Regression
}

\description{
Fit negative binomial (NB) regression models with NB1 response distribution.
}

\usage{
negbin1(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = TRUE,
  control = negbin1_control(...), ...)

negbin1_fit(y, x, control)

negbin1_control(maxit = 5000, start = NULL, ...)
}

\arguments{
  \item{formula}{FIXME}
  \item{data}{FIXME}
  \item{subset}{FIXME}
  \item{na.action}{FIXME}
  \item{model}{FIXME}
  \item{y}{FIXME}
  \item{x}{FIXME}
  \item{control}{FIXME}
  \item{\dots}{FIXME}
  \item{maxit, start}{control arguments passed to \code{\link[stats]{optim}}.}
}

\details{
FIXME
}

\value{
FIXME
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\references{
Winkelmann and Boes
}

\seealso{
\code{\link[MASS]{glm.nb}}
}

\examples{
formals(negbin1)
}

\keyword{regression}

