\name{GmmEst}
\alias{GmmEst}

\title{Generalized Methods of Moments Estimation}

\description{
Fitting a model by GMM
}
\usage{
GmmEst(func=NULL, initial_vals=NULL, data=NULL, 
        esttype=c("2step","1step","iter"), 
		initial_W=NULL,
        maxiter=100, ...)
}
\arguments{
  \item{func}{a user supplied function with the first input being a parameter vector. The function is assumed to return a (nobs x kmoms) matrix of moments. Necessary input.}
  \item{initial_vals}{a vector of initial values of the parameters. Necessary input.}
  \item{data}{a dataframe or matrix containing the necessary data to compute the moments. Necessary input.}
  \item{esttype}{One of "2step", "1step" or "iter". Default value is "2step".}
  \item{initial_W}{An initial weighting matrix of class "matrix" that can be supplied by the user. If not supplied, the identity matrix is used.}
  \item{maxiter}{If esttype=="iter", maxiter is used to determine the maximal number of iterations. Iterations stops if the change in the norm of the estimated parameter vector is smaller than 1e-03 or the number specified in maxiter is exceeded. Default value is 100.}
  \item{...}{Further arguments that can be passed to R's \code{"optim"} function.}
}
\details{Will be available soon}
\value{An object of class \code{"GmmEst"}.}
\references{Will be available soon}
\seealso{Will be available soon}

\examples{
# Description
# A simple example of a linear regression performed using the package

## Data
set.seed(23)
n = 500
x = rnorm(n, mean = 5, sd = 3)
eps = rnorm(n, mean = 0, sd = 3)
y = 10 + 5*x + eps
df = data.frame(y,x)

## User specified function
func1 = function(params, data){
  b0 = params[1]
  b1 = params[2]
  
  y = data$y
  x = data$x
  
  yfit = b0 + b1*x
  eps = y - yfit
  
  gt1 = eps
  gt2 = eps*x
  
  gt = cbind(gt1, gt2)
  return(gt)
}

## Model fitting
mod = GmmEst(func1, initial_vals = c(5,20), data=df)

}
\keyword{regression}
