\name{nndd}
\alias{nndd}
\alias{nndd_reshape}
\alias{nndd_reshape_other}

\title{Nearest Neighbour Matching (NN) followed by a Linear Model with Difference in Differences (DD)}
\description{Constructs NN by estimating pscores with a GLM and estimates a linear model in a DD setting}
\usage{


nndd(formula, data, index = c("firm_id", "year"),  t_time , nn_time,
                     time_ids = c("year", ""),
                     link = "logit",
                     subset , na.action, clustervariables, 
                     model = TRUE, y = TRUE, x = FALSE, displ_coefs, 
                      ...)


nndd_reshape(object_nndd)

nndd_reshape_other(call, data)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    a formula expression of the form \code{tg | outcome ~ x | z} where
    \code{tg} is the response and the \code{z} regressor variable of the GLM. \code{outcome} is the response and \code{x} the regressor variable of the DD model.}

  \item{data}{a data frame containing the variables occurring in the
    formulas such as time and group identifiers. The data has to be a panel.}
  \item{index}{a list containing the name of the group, time identifier.}

  \item{nn_time}{a list containing the start timing and end timing for the NN estimation. The GLM is estimated at the end timing. However, it contains all regressor variables lagged in the given interval (start , end). In the case of omitting \code{nn_time}, the GLM is estimated at the given \code{t_time}. 
}
  \item{t_time}{a string containing the timing of the treatment.}
  \item{time_ids}{#Fixme: might be implemented in future. }
  \item{link}{family for the GLM either binomial or logistic. 
}
\item{clustervariables}{list of variables on which the standard errors should be clustered.}
 \item{subset}{an optional vector specifying a subset of observations to be
    used for fitting.}
 \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.}
\item{model}{logical. If \code{TRUE} \emph{model frame} is
    included as a component of the returned value.}
  \item{x, y}{#Fixme: do I need this?}
 \item{displ_coefs}{a list of coefficients being displayed of the DD if a summary function is applied}   
\item{object_nndd}{an object of class \code{nndd}}
 \item{\dots}{arguments to be passed on to the estimation of the DD (\link{lm.fit})}
\item{call}{a call of a \code{nndd} object}  
  }

\value{
 \item{nndd}{An object of class \code{nndd}.}
 \item{nndd_reshape}{A data.frame containing a reshaped model matrix of the \code{nndd} object.}
 \item{nndd_reshape_other}{A data.frame containing reshaped data.}
}
\references{
Angrist JD, Pischke JS (2008). Mostly Harmless Econometrics: An Empiricist's Companion.
Princeton University Press.

Caliendo M, Kopeinig S (2008). "Some Practical Guidance for the Implementation of
Propensity Score Matching." Journal of Economic Surveys, 22(1), 31--72.

Rubin DB (1973). "Matching to Remove Bias in Observational Studies." Biometrics, 29(1),
159--183. 

}



\examples{

data("IRCA", package = "nndd")
library(Formula)
formula <- Formula(treated | v_crime ~ officers_pc + income + p_crime 
						   | officers_pc + income + p_crime + crack_index 
						   + unemprate) 
						   
nndd <- nndd(formula = formula, data = IRCA, 
			index = c("county", "year"), 
			t_time = "1986" )
names(nndd)
print(nndd)
plot(nndd, IRCA)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{matchin }
\keyword{difference in difference }% __ONLY ONE__ keyword per line
\keyword{nndd}