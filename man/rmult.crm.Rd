\name{rmult.crm}
\alias{rmult.crm}
\title{
Simulating Correlated Ordinal Multinomial Responses under a Marginal Continuation Ratio Model
}
\description{
Simulates correlated ordinal multinomial responses assuming that a continuation ratio model holds for the marginal probabilities.
}
\usage{
rmult.crm(clsize, lin.pred, cor.matrix, intercepts, link = "probit")
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
numerical vector or matrix containing the intercepts of the marginal continuation ratio model. See details for more info.
}
  \item{link}{
character string indicating the link function of the continuation ratio model. Options include \code{"probit"}, \code{"logit"}, \code{"cloglog"} or \code{"cauchit"}.
}
}
\details{
The assumed marginal cumulative link model is  
\deqn{Pr(Y_{it}=j |Y_{it} \ge j,x_{it})=F(\beta_{t0j} +\beta^{'}_{t} x_{it})}
where \eqn{F} is the cumulative distribution function determined by \code{link}. For subject \eqn{i}, \eqn{Y_{it}} is the \eqn{t}-th multinomial response and \eqn{x_{it}} is the associated covariates vector. Finally, \eqn{\beta_{t0j}} is the \eqn{j}-th category specific intercept and \eqn{\beta^{'}_{t}} is the regression parameter vector corresponding to the \eqn{t}-th measurement occasion.

The multinomial responses \eqn{Y_{it}} are determined by extending the latent variable threshold approach of Tutz (1991) to the correlated latent continuous variables. See vignette for more info.

The sample is equal to the number of rows in \code{lin.pred}. Row \eqn{i} corresponds to subject \eqn{i}. The (\eqn{i,t})-th element of \code{lin.pred} equals \eqn{\beta^{'}_{t} x_{it}}. 

When \eqn{\beta_{0tj}=\beta_{0j}} for all \eqn{t}, then \code{intercepts} should be provided as a numerical vector. Otherwise, \code{intercepts} must be a numerical matrix such that \eqn{t} contains the category specific intercepts at the \eqn{t}-th measurement occasion.
}
\value{
\item{Ysim}{the correlated ordinal multinomial responses. Element (\eqn{i},\eqn{t}) represents the realization of \eqn{Y_{it}}.}
\item{correlation}{the correlation matrix of the multivariate normal distribution in the NORTA method.}
\item{rlatent}{the latent random vectors.}
}
\references{
Cario, M. C. and Nelson, B. L. (1997)  \emph{Modeling and generating random vectors with arbitrary marginal distributions and correlation matrix}. Technical Report, Department of Industrial Engineering and Management Sciences, Northwestern University, Evanston, Illinois.

Li, S. T. and Hammond, J. L. (1975)  Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557--561.

Tutz, G. (1991) Sequential models in categorical regression. \emph{Computational Statistics & Data Analysis}, \bold{11}, 275--295.
}
\author{
Anestis Touloumis
}
\seealso{
For simulating correlated ordinal multinomial responses see also \code{rmult.clm}.

For simulating correlated nominal multinomial responses use \code{rmult.bcl}, and for simulating correlated binary responses use \code{rbin}. 
}
\examples{
## See vignette for more details about the sampling scheme.
set.seed(1)
N <- 500
clsize <- 4
intercepts <- c(-1.5, -0.5, 0.5, 1.5)
cor.matrix <- diag(1, 16)
x <- rnorm(N)
lin.pred <- matrix(rep(x,clsize), N, clsize, byrow=TRUE)
Y <- rmult.crm(clsize, lin.pred, cor.matrix, intercepts, link="probit")
}

