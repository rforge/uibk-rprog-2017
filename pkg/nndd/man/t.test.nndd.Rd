\name{S3 Methods in nndd}
\alias{t.test.nndd}
\alias{nndd_ttest}
\alias{predict.nndd}
\alias{predict_fun_glm}
\alias{print.nndd}
\alias{plot.nndd}
\alias{vcov.nndd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ttest}
\description{
Performs a ttest on a \code{"nndd"} object. The ttest is applied on only on the time span on which the NN are specified.}
\usage{
\method{t.test}{nndd}(x, \dots)

nndd_ttest(ttest_data, indexes = c("year", "firm_id", "tg"), nn_time, \dots)

\method{predict}{nndd}(object, prediction = c("nn", "dd"), newdata, \dots)

predict_fun_glm(pformula, link = "logit", data)

\method{print}{nndd}(x, digits = max(3, getOption("digits") - 3), \dots)

\method{plot}{nndd}(x, data, which = c(1L:7L, 8L), ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...)

\method{vcov}{nndd}(object, \dots)

}
\arguments{
  \item{x, object}{an object of class \code{"nndd"}}
  
  \item{ttest_data}{a data.frame containing the data on which the ttest should be performed on.}

\item{indexes}{a list containing the name of the time, group, treatment identifier, and the outcome variable #Fixme: last two not really necessary}
 
\item{nn_time}{a list containing the start timing and end timing for the NN estimation. The GLM is estimated at the end timing. However, it contains all regressor variables lagged in the given interval (start , end). 
}

\item{pformula}{
    a formula expression of the form \code{tg | outcome ~ x | z} where
    \code{tg} is the response and \code{z} regressor variable of the GLM. \code{outcome} is the response and \code{x} the regressor variable of the DD model.}
\item{prediction}{a list containing specifying which model should be applied for the prediction. If the list contains "nn" the pscores of the GLM model are predicted. In the case of containing "dd" the prediction of the linear difference in differences model is returned. Both predictions are returned if the list is specified as "nn", "dd".}

\item{newdata}{a data.frame containing the data on which the prediction should be based on.}

\item{link}{family for the GLM either binomial or logistic.}

\item{digits}{see \code{\link{print}}}    

\item{data}{a data frame containing the variables occurring in the
    formulas such as time and group identifiers. The data has to be a panel.}
\item{\dots}{space for objects to be passed on}

\item{which}{if a subset of the plots is required, specify a subset of the numbers 1:8.}

\item{ask}{logical; if TRUE, the user is asked before each plot, see \link{par}(ask=.).}

}
\details{...}
\value{...}
\references{...}
\examples{x <- "Fix: me"
print(x)}
\keyword{ttest}