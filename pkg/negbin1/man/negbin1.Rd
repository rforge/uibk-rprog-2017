\name{negbin1}
\alias{negbin1}
\alias{negbin1_fit}
\alias{negbin1_control}

\title{
Negative Binomial 1 Regression
}

\description{
Fit negative binomial (NB) regression models with NB1 response distribution. The negative binomial distribution can arise as a gamma mixture of Poisson distributions, and can be used for modelling over-dispersed count data.
}

\usage{
negbin1(formula, data, subset, na.action,
  model = TRUE, y = TRUE, x = TRUE,
  control = negbin1_control(...), ...)

negbin1_fit(y, x, control)

negbin1_control(maxit = 5000, start = NULL, ...)
}

\arguments{
  \item{formula}{an object of class "\code{\link[Formula]{formula}}": a symbolic description of the model to be fitted.}
  \item{data}{a data frame containing the variables of the model. If not found in the data, variables are taken from \code{environment(formula)}}.
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{na.action}{a function which indicates what should happen when the data contain NAs.}
  \item{model, y, x}{logicals. If TRUE the corresponding components of the fit are returned.}
  \item{control}{a list of parameters for controlling the fitting process. For \code{negbin1_fit} this is passed to \code{negbin1_control}.}
  \item{\dots}{additional arguments to be passed.}
  \item{maxit, start}{control arguments passed to \code{\link[stats]{optim}}.}
}

\details{
A continous mixture of Poisson distribution with Gamma weights is called negative binomial distribution.
To employ the NB1 distribution in a regression, the natural mean equation is employed:

\eqn{log(\mu_{i}) = x_i^\top \beta}.

For the variance, the NB1 model assumes varying \eqn{\theta_i} such that \eqn{\alpha = \mu_i/\theta_i} and thus 

\eqn{Var(y_i | x_i) = (1 + \alpha) \cdot \mu_i}.

The workhorse function is \code{negbin1_fit}, which is normally not called directly, but when the model response and model matrix have already been calculated. Starting values in the optimization are by default taken from a Poisson glm.
}

\value{
\code{negbin} returns an object of \code{\link{class}} "\code{negbin}".
}

\references{
Cameron AC & Trivedi PK (1986).
Econometric Models Based on Count Data: Comparisons and Applications of Some Estimators and Tests,
\emph{Journal of Applied Econometrics}, \bold{1}, 29--53.

Cameron AC & Trivedi PK (2013).
\dQuote{Regression Analysis of Count Data},
Cambridge University Press.

Lawless JF (1987).
Negative Binomial and Mixed Poisson Regression,
\emph{The Canadian Journal of Statistics}, \bold{15}(3), 209--225.

Winkelmann R & Boes S (2009).
\dQuote{Analysis of Microdata},
Springer, Second Edition.
}

\seealso{
\code{\link[MASS]{glm.nb}}, \code{\link[VGAM]{vglm}}
}

\examples{
## packages
require("Formula") 
require("VGAM")
require("lmtest")
require("MASS")

## data generating process
dgp <- function(n = 1000, coef = c(2, 3, 0, 0.7)) {
  d <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
    )
  d$mu <- exp(coef[1] + coef[2] * d$x1 + coef[3] * d$x2)
  d$y <- rnbinom(n, mu = d$mu, size = coef[4] * d$mu)
  return(d)
}

## simulate data
set.seed(2007-05-15)
d <- dgp()

## model (with function negbin1)
m1 <- negbin1(y ~ x1 + x2, data = d, control = list(maxit = 5000))
    
## model (with vglm from VGAM package)
m2 <- vglm(y ~ x1 + x2, negbinomial(parallel = TRUE, zero = NULL), data = d, trace = TRUE)    
summary(m2)

## model (with glm.nb from MASS package)     
m3 <- glm.nb(y ~ x1 + x2, data = d) 

## comparison of coefficient vectors
cbind(m1$coefficients$location, coef(m2, matrix = TRUE)[1:3,1], coef(m3))   

## comparison of log likelihoods
cbind(m1$loglik, logLik(m2))

## model comparison
m0 <- negbin1(y ~ x1, data = d, control = list(maxit = 5000))
AIC(m0, m1)
BIC(m0, m1)
lrtest(m0, m1)

## comparison of theta
thetam1 <- m1$coefficients$theta
co <- coef(m2, matrix = TRUE)
thetam2 <- exp(co["(Intercept)", "loge(size)"] - co["(Intercept)", "loge(mu)"]) 
    
## benefit of formula/terms interface
update(m1, subset = x2 > 0)
head(model.frame(m1))
head(model.matrix(m1))

## predictions
newd <- data.frame(x1 = c(2, 0, 1, 4), x2 = c(1, 0, -1, 1))
predict(m1, newd, type = "location")
predict(m1, newd, type = "response")
predict(m1, newd, type = "quantile")
predict(m1, newd, type = "probability", at = 0)    

}

\keyword{regression}

