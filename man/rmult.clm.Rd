\name{rmult.clm}
\alias{rmult.clm}
\title{
Simulating Correlated Ordinal Multinomial Responses under a Marginal Cumulative Link Model 
}
\description{
Simulates correlated ordinal multinomial responses assuming that a cumulative link model holds for the marginal probabilities.
}
\usage{
rmult.clm(clsize, lin.pred, cor.matrix, intercepts, link = "probit")
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
numerical vector or matrix containing the intercepts of the marginal cumulative link model. See details for more info.
}
  \item{link}{
character string indicating the link function of the cumulative link model. Options include \code{"probit"}, \code{"logit"}, \code{"cloglog"} or \code{"cauchit"}.
}
}
\details{
The assumed marginal cumulative link model is  
\deqn{Pr(Y_{it}\le j |x_{it})=F(\beta_{t0j} +\beta^{'}_{t} x_{it})}
where \eqn{F} is the cumulative distribution function determined by \code{link}. For subject \eqn{i}, \eqn{Y_{it}} is the \eqn{t}-th multinomial response and \eqn{x_{it}} is the associated covariates vector. Finally, \eqn{\beta_{t0j}} is the \eqn{j}-th category specific intercept and \eqn{\beta^{'}_{t}} is the regression parameter vector corresponding to the \eqn{t}-th measurement occasion.

The multinomial response \eqn{Y_{it}} takes the value \eqn{j} if and only if the corresponding latent random variable lies on \eqn{(\beta_{t0(j-1)}-\beta^{'}_{t}x_{it},\beta_{t0j}-\beta^{'}_{t}x_{it}]}. See vignette for more info.

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

McCullagh, P. (1980) Regression models for ordinal data. \emph{Journal of the Royal Statistical Society B}, \bold{42}, 109--142.

Touloumis, A., Agresti, A. and Kateri, M. (2013) GEE for multinomial responses using a local odds ratios parameterization. \emph{Biometrics}, \bold{69}, 633--640.
}
\author{
Anestis Touloumis
}
\seealso{
For simulating correlated ordinal multinomial responses see also \code{rmult.crm}.

For simulating correlated nominal multinomial responses use \code{rmult.bcl} and for simulating correlated binary responses use \code{rbin}.
}
\examples{
## See vignette for the sampling scheme.
set.seed(12345)
N <- 500
clsize <- 4
intercepts <- c(-1.5, -0.5, 0.5, 1.5)
cor.matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
lin.pred <- rsmvnorm(N, toeplitz(c(1, rep(0.85, clsize-1))))
Y <- rmult.clm(clsize, lin.pred, cor.matrix, intercepts, "probit")$Ysim
data <- cbind(c(t(Y)), c(t(lin.pred)))
data <- data.frame(data)
data$id <- rep(1:N, each=clsize)
data$time <- rep(1:clsize, N)
colnames(data) <- c("y", "x", "id", "time")

## Fit a GEE model (Touloumis et al., 2013) to estimate the regression coefficients.
library(multgee)
fitmod <- ordLORgee(y~x, data=data, id=id, repeated=time, link="probit", IM="cholesky")
coef(fitmod)
}
