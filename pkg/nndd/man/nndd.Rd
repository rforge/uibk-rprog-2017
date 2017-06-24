\name{nndd}
\alias{nndd}
\alias{nndd_reshape}
\alias{nndd_reshape_other}

\title{Nearest Neighbour Matching (NN) followed by a Linear Model with Difference in Differences (DD)}
\description{Constructs NN by estimating pscores with a GLM and estimates a linear model in a DD setting}
\usage{


nndd(formula, data, indexes = c("year", "firm_id", "tg", "outcome"),  t_time , nn_time,
                     time_ids = c("year", ""),
                     link = "logit",
                     subset , na.action, clustervariables, 
                     model = TRUE, y = TRUE, x = FALSE,
                      ...)


nndd_reshape(object_nndd)

nndd_reshape_other(call, data)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    a formula expression of the form \code{tg | outcome ~ x | z} where
    \code{tg} is the response and \code{z} regressor variable of the GLM. \code{outcome} is the response and \code{x} the regressor variable of the DD model.}

  \item{data}{a data frame containing the variables occurring in the
    formulas such as time and group identifiers. The data has to be a panel.}
  \item{indexes}{a list containing the name of the time, group, treatment identifier, and the outcome variable #Fixme: last two not really necessary}

  \item{nn_time}{a list containing the start timing and end timing for the NN estimation. The GLM is estimated at the end timing. However, it contains all regressor variables lagged in the given interval (start , end). 
}
  \item{t_time}{a string containing the timing of the treatment.}
  \item{time_ids}{#Fixme: delete this}
  \item{link}{family for the GLM either binomial or logistic. 
}
\item{clustervariables}{List of variables on which the standard errors should be clustered.}
 \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
 \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
\item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{#Fixme: do I need this?}
  \item{\dots}{arguments to be passed on to the estimation of the DD (\link{lm.fit})}
\item{object_nndd}{an object of class \code{"nndd"}}
\item{call}{a call of a \code{"nndd"} object}  

  }
\details{
%%  ~~ If necessary, more details than the description above ~~

\code{plot.nndd} Creates eight plots in one window. The first tow plots are the distribution of the treated and control pscores pre and post the Nn selection. If data is not specified the pre Nn plot is omitted. In the following the outcome variable of the Dd estimation is plotted over time. For a description of the last four plots see \code{\link{plot.lm}}.
}
\value{
 \item{nndd}{An object of class \code{"nndd"}.}
 \item{nndd_reshape}{A data.frame containing a reshaped model matrix of the \code{"NdDd"} object.}
 \item{nndd_reshape_other}{A data.frame containing reshaped data.}
 \item{predict_fun_glm}{A list containing pscores.}
 \item{predict.nndd}{A data.frame containing the specified predictions.}
 \item{print.nndd}{Returns the imputed object.}
 \item{plot.nndd}{#Fixme: Ns are not well chosen at the Density plots (the N is of the control) Name of the plot is returned.}
 
}
\references{
Fixme
}
\author{
[aut]
}
\note{
Fixme
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Fixme:
}
\examples{

data("IRCA", package = "nndd")
library(Formula)
formula <- Formula(treated | v_crime ~ officers_pc + income + p_crime 
						   | officers_pc + income + p_crime + crack_index + unemprate) 
nndd <- nndd(formula = formula, data = IRCA, indexes = c("year", "county", "treated", "v_crime"), t_time = "1986" )
names(nndd)
print(nndd)
plot(nndd, IRCA)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ kwd1 }
\keyword{ kwd2 }% __ONLY ONE__ keyword per line
