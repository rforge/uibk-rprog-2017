\name{IRCA}
\alias{IRCA}

\title{Data on the Evaluation of the Immigration Reform and Control Act (IRCA)}

\description{
The data is adapted data of Baker(2015). The author used the original data to examine the effects of the Immigration Reform and Control Act on crime. The IRCA was implemented in 1986 and forbid to hire or recruit undocumented immigrants. However the IRCA also implemented a near-universal legalization of immigrants in the United States.
The theory behind a positive impact of the IRCA on crime is that an increased labour market opportunity due to IRCA increases legal work and decreases crime. The labour market opportunity is expected to increase because legal (documented) immigrants have a higher salary and lower change to be fired. 
}

\usage{data("IRCA")}

\format{
A data frame containing 31206 observations on 21 variables. In detail it is a balanced data panel of 1486 US counties over 21 years (the time span is 1980 till 2000).
\describe{
  \item{county}{Observation county.}
  \item{year}{Observation year.}
  \item{pop}{Log population.}
  \item{immigrants_pc}{CumWeightedAllImmPerCap.}
  \item{unemprate}{Unemployment rate.}
  \item{povrate}{Poverty rate.}
  \item{immigrant_all}{AllTimeWeightedImm.}
  \item{crack_index}{Crack Index}
  \item{abortions}{Number of abortions lag 13 per capita.}
  \item{crimes}{The logarithm of the total amount of all crimes per capita.}
  \item{v_crime}{The logarithm of the total amount of violent crimes per capita.}
  \item{p_crime}{The logarithm of the total amount of property crimes per capita.}
  \item{officers_pc}{Number of officers per capita.}
  \item{income}{The logarithm of the average income.}  
  \item{empl}{The logarithm of the employment.}
  \item{sim_CumWeightedAllImmPerCap}{sim_CumWeightedAllImmPerCap. }
  \item{sim_diff}{sim_diff.}
  \item{treated}{Treated states. Defined as stated which  immigrants_pc was larger then 0 after the reform.}
  \item{post}{Identifier for the post treatment period.}
  \item{StateFIPS}{StateFIPS}
  \item{CountyFIPS}{CountyFIPS}
}
}

\source{
Adapted supplementary material of Baker(2015).
}

\references{
Baker, S. R. (2015). \emph{Effects of Immigrant Legalization on Crime}, American Economic Review, 105(5): 210-213.

}

\examples{
data("IRCA", package = "nndd")
library(Formula)
formula <- Formula(treated | v_crime ~ officers_pc + income + p_crime | 
					officers_pc + income + p_crime + crack_index + unemprate) 
					
ma <- nndd(formula = formula, data = IRCA, 
			indexes = c("year", "county", "treated", "v_crime"), 
			nn_time = c("1985","1985" ), t_time = "1986" )

ma
#plot(ma, IRCA)
}

\keyword{datasets}
