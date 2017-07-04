\name{GmmEst}
\alias{GmmEst}
\alias{GmmEst_control}

\alias{coef.GmmEst}
\alias{nobs.GmmEst}
\alias{vcov.GmmEst}
\alias{print.GmmEst}
\alias{summary.GmmEst}
\alias{print.summary.GmmEst}

\alias{estfun.GmmEst}
\alias{bread.GmmEst}

\title{Generalized Methods of Moments Estimation}

\description{
Function to estimate a parameter vector by one-step, two-step or iterative GMM}
\usage{
GmmEst(func, theta0, data,
    est_type=c("2step","1step","iter"), 
    func_jac=NULL, initial_W=NULL, 
    crit=10e-7, itermax=100, 
    optim_method=c("BFGS","Nelder-Mead","L-BFGS-B"),
    control = GmmEst_control(\dots), \dots)

GmmEst_control(maxit = 5000, ...)

}
\arguments{
  \item{func}{a user supplied \code{function} of the form \emph{f(theta,data)}, with the first input being a vector of the parameters, which are to be estimated. The function is assumed to return a (nobs x kmoms) matrix of moments.}
  \item{theta0}{a vector of initial values of the parameters.}
  \item{data}{a dataframe or matrix containing the necessary data to compute the moments in \code{func}.}
  \item{est_type}{\emph{Optional:} one of \emph{2step}, \emph{1step} or \emph{iter}. Default value is \emph{2step}.}
  \item{func_jac}{\emph{Optional:} a user supplied \code{function} of the form \emph{j(theta,data)}, with the first input being a vector of the parameters which are to be estimated. The function is assumed to return the expected / average jacobian of the moment conditions. Should be matrix of size (kmoms x kparam).}
  \item{initial_W}{\emph{Optional:} an initial weighting matrix of class \emph{matrix} that can be supplied by the user. If not supplied, the identity matrix is used.}
  \item{crit}{\emph{Optional:} the stopping rule for the iterative GMM. It can be reduce to increase the precision.}
  \item{itermax}{\emph{Optional:} the maximal number of iterations in the iterated GMM procedure.}
  \item{optim_method}{\emph{Optional:} the possible implemented methods in \code{optim} that are used to minimize the objective function. Default is \emph{BFGS}.}
  \item{control, maxit}{\emph{Optional:} a list of control parameters passed to \code{\link[stats]{optim}}. \code{maxit} defines the maximal number of iterations in \code{\link[stats]{optim}}.}
  \item{\dots}{further arguments that can be passed to \code{\link[stats]{optim}}.}
}
\details{Please have a look into the vignette for }
\value{An object of class \code{"GmmEst"}.

The object of class "gmm" is a list containing at least:

\item{coefficients}{(kparams by 1) vector of coefficients.}
\item{jstat}{a list with elements \code{value} and {pval} returning the value of the J-Statistic (overidentification test). If the model is overidentified, the corresponding p-value is included, otherwise it is set to \code{NA}.}
\item{nobs}{number of observations.}
\item{vcov}{variance-covariance matrix of the estimated parameters.}
\item{kparams}{number of estimated parameters.}
\item{kmoms}{number of moment conditions.}
\item{est_type}{string of the applied estimation type.}
\item{S}{the estimated S matrix}
\item{W}{the inverse of the S matrix.}
\item{dmat}{the d-matrix.}
\item{gt}{an (nobs by kmoms) matrix of moment condition, based on the estimated parameters.}
\item{gt_mean}{the (kmoms by 1) vector of moment condition, based on the estimated parameters.}
\item{niter}{If \code{est_type=='iter'}, the number of performed GMM iterations.}
\item{vcov_gt}{The variance-covariance matrix of the moment condition matrix}
\item{identification}{String that shows if the model is over-, just-, or underidentified.}
}

\references{
  Cochrane, J.H. (2005),
  Asset Pricing. 
  Revised edition. Princeton University Press

  Greene, W.H. (2011),
  Econometric Analysis. 
  7. edition. Prentice Hall. 

  Hall, A.R. (2005),
  Generalized Method of Moments.
  1. edition. Oxford University Press. 

  Hansen, L.P. (1982),
  Large Sample Properties of Generalized Method of Moments Estimators.
  \emph{Econometrica}, \bold{50},
  1029-1054,
}

\examples{
# ====================================
# A minimal example: GMM can do OLS.

# Load data
data(mtcars)

## User specified function
func = function(params, data){
  b0 = params[1]
  b1 = params[2]
  
  y = data$mpg
  x = data$hp
  
  yfit = b0 + b1*x
  eps = y - yfit
  
  gt1 = eps
  gt2 = eps*x
  
  gt = cbind(gt1, gt2)
  return(gt)
}

# Estimation
mod = GmmEst(func, c(20,0), mtcars)

# Print result
summary(mod)


# ====================================
# Second example: Consumption based Asset pricing model

# Load dataset from within the package
data("data_consumption",package = "GmmEst")

# Define user specific function FUNC 
# (moment conditions using rmrf and hml portfolio returns as instruments)
mom_cond = function(params, data){
  gamma = params[1]
  
  beta = 1
  dc = data$dc
  rmrf = data$rmrf
  hml = data$hml
  rf = data$rf
  
  gt1 = beta*dc^(-gamma)*rmrf
  gt2 = beta*dc^(-gamma)*hml
  gt = cbind(gt1, gt2)
  return(gt)
}

# Define FUNC_JAC (in this case the gradient of the moment conditions)
mom_cond_grad = function(params,data){
  gamma = params[1]
  
  beta = 1
  dc = data$dc
  rmrf = data$rmrf
  hml = data$hml
  
  d1 = -dc^(-gamma)*log(dc)*rmrf*beta
  d2 = -dc^(-gamma)*log(dc)*hml*beta
  d = c(mean(d1),mean(d2))
  return(d)
}

# Estimate the gamma parameter
mod2 = GmmEst(mom_cond, 100, data, est_type = "2step", 
             optim_method = 'BFGS', func_jac = mom_cond_grad)
summary(mod2)

}
\keyword{GMM, regression}
