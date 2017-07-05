\name{S3 Methods of nndd}
\alias{t.test.nndd}
\alias{nndd_ttest}
\alias{predict.nndd}
\alias{predict_fun_glm}
\alias{print.nndd}
\alias{plot.nndd}
\alias{vcov.nndd}
\alias{getSummary.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{S3 Methods of nndd}
\description{
Methods for testing, predicting, printing, plotting for a \code{nndd} object.}
\usage{
\method{t.test}{nndd}(x, \dots)

nndd_ttest(ttest_data, indexes = c("year", "firm_id", "tg"), nn_time, \dots)

\method{predict}{nndd}(object, prediction = c("nn", "dd"), newdata, \dots)

predict_fun_glm(pformula, link = "logit", data)

\method{print}{nndd}(x, digits = max(3, getOption("digits") - 3), \dots)

\method{plot}{nndd}(x, data, which = c(1L:7L, 8L), 
					ask = prod(par("mfcol")) < length(which) 
					&& dev.interactive()
					, ...)

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

\item{which}{if a subset of the plots is required, specify a subset of the numbers 1:8.}

\item{ask}{logical; if TRUE, the user is asked before each plot, see \link{par}(ask=.).}
\item{\dots}{space for objects to be passed on}

}
\details{
\code{t.test} is applied only on the time span on which the NN of the \code{nndd} object are specified.
\code{plot.nndd} Creates eight plots in one window. The first tow plots are the distribution of the treated and control pscores pre and post the NN selection. If data is not specified the pre NN plot is omitted. In the following the outcome variable of the DD estimation is plotted over time. For a description of the last four plots see \code{\link{plot.lm}}.

\code{\link[memisc]{getSummary}} (from the \pkg{memisc} package, enabling \code{\link[memisc]{mtable}}).
  }


\value{\item{predict_fun_glm}{A list containing pscores.}
 \item{predict.nndd}{A data.frame containing the specified predictions.}
 \item{print.nndd}{Returns the imputed object.}
 \item{plot.nndd}{#Fixme: Ns are not well chosen at the Density plots (the N is of the control) Name of the plot is returned.}
 }
\references{...}
\examples{
data("IRCA", package = "nndd")
library(Formula)
formula <- Formula(treated | v_crime ~ officers_pc + income  | 
                     officers_pc + income  + crack_index + unemprate) 

ma <- nndd(formula = formula, data = IRCA, 
           indexes = c("year", "county", "treated", "v_crime"), 
           nn_time = c("1985","1985" ), t_time = "1986" )

ma
t.test(ma)
}
\keyword{S3 Methods}