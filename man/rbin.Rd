\name{rbin}
\alias{rbin}
\title{
Simulating Correlated Binary Responses under Marginal Model Specification
}
\description{
Simulates correlated binary responses assuming that a prespecified model holds for the marginal probabilities.
}
\usage{
rbin(clsize, lin.pred, cor.matrix, intercepts, link = "probit")
}
\arguments{
  \item{clsize}{
integer indicating the common cluster size.
}
  \item{lin.pred}{
matrix with \code{clsize} columns. See details for more info.
}
  \item{cor.matrix}{
matrix indicating the correlation matrix of the multivariate normal distribution in the NORTA method.
}
  \item{intercepts}{
numerical (or numeric vector of length \code{clsize}) containing the intercept(s) of the marginal model. See details for more info.
}
  \item{link}{
character string indicating the link function of the marginal model. Options include \code{"probit"}, \code{"logit"}, \code{"cloglog"} or \code{"cauchit"}.
}
}
\details{
The assumed marginal model is  
\deqn{Pr(Y_{it} = 1 |x_{it})=F(\beta_{t0} +\beta^{'}_{t} x_{it})}
where \eqn{F} is the cumulative distribution function determined by \code{link}. For subject \eqn{i}, \eqn{Y_{it}} is the \eqn{t}-th multinomial response and \eqn{x_{it}} is the associated covariates vector. Finally, \eqn{\beta_{t0}} and \eqn{\beta_{t}} is the intercept and the regression parameter vector at the \eqn{t}-th measurement occasion.

The binary variable \eqn{Y_{it}} takes the value 1 if and only if the corresponding latent random variable is less than or equal to \eqn{\beta_{0t}+\beta^{'}_{t}x_{it}}. Otherwise \eqn{Y_{it}} is set equal to 0. See vignete for more info.

The sample is equal to the number of rows in \code{lin.pred}. Row \eqn{i} corresponds to subject \eqn{i}. The (\eqn{i,t})-th element of \code{lin.pred} equals \eqn{\beta^{'}_{t} x_{it}}. 

When \eqn{\beta_{t0}=\beta_{0}} for all \eqn{t}, then \code{intercepts} should be provided as a single number. Otherwise, \code{intercepts} must be provided as a numerical vector such that the \eqn{t}-th element is the intercept at measurement occasion \eqn{t}.
}
\value{
\item{Ysim}{the correlated nominal multinomial responses. Element (\eqn{i},\eqn{t}) represents the realization of \eqn{Y_{it}}.}
\item{correlation}{the correlation matrix of the multivariate normal distribution in the NORTA method.}
\item{rlatent}{the latent random vectors.}
}
\references{
Cario, M. C. and Nelson, B. L. (1997)  \emph{Modeling and generating random vectors with arbitrary marginal distributions and correlation matrix}. Technical Report, Department of Industrial Engineering and Management Sciences, Northwestern University, Evanston, Illinois.

Li, S. T. and Hammond, J. L. (1975)  Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557--561. 
}
\author{
Anestis Touloumis
}
\seealso{
For simulating correlated nominal multinomial responses use \code{rmult.bcl} and for simulating correlated ordinal multinomial responses use \code{rmult.clm} or \code{rmult.crm}.
}
\examples{
## See vignette for the sampling scheme.
set.seed(1)
N <- 500
clsize <- 4
intercepts <- 1
cor.matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
lin.pred <- matrix(rnorm(N), N, clsize)
Y <- rbin(clsize, lin.pred, cor.matrix, intercepts, "probit")   
head(Y$Ysim)               
}
