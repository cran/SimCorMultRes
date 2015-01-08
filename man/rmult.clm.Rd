\name{rmult.clm}
\alias{rmult.clm}
\title{
Simulating Correlated Ordinal Multinomial Responses Under A Marginal Cumulative Link Model 
}
\description{
This function simulates correlated ordinal multinomial responses assuming that a cumulative link model holds for the marginal probabilities.
}
\usage{
rmult.clm(clsize, lin.pred, corr, cuts, link = "probit")
}
\arguments{
  \item{clsize}{
integer indicating the common cluster size.
}
  \item{lin.pred}{
matrix with \code{clsize} columns. See details for more info.
}
  \item{corr}{
matrix or constant describing the latent correlation structure. See details for more info.
}
  \item{cuts}{
vector containing the intercepts of the cumulative link model, and \code{-Inf} and \code{Inf} as first and last elements respectively.
}
  \item{link}{
character string indicating the link function of the cumulative link model. Options include \code{"probit"}, \code{"logit"}, \code{"cloglog"} or \code{"cauchit"}.
}
}
\details{
The assumed marginal cumulative link model is  
\deqn{Pr(Y_{it}\le j |x_{it})=F(\beta_{0j} +\beta^{'} x_{it})}
where \eqn{Y_{it}} is the \eqn{t}-th multinomial response for subject \eqn{i}, \eqn{x_{it}} is the associated covariates vector, \eqn{F} is the cumulative distribution function determined by \code{link}, \eqn{\beta_{0j}} is the \eqn{j}-th response category specific intercept and \eqn{\beta} is the marginal regression parameter vector excluding intercepts.

The multinomial response \eqn{Y_{it}} takes the value \eqn{j} if and only if the corresponding latent random variable lies on \eqn{(\beta_{0(j-1)}-\beta^{'}x_{it},\beta_{0j}-\beta^{'}x_{it}]}.

The row dimension of \code{lin.pred} equals the sample size. Row \eqn{i} corresponds to subject \eqn{i}. The (\eqn{i,t})-th element of \code{lin.pred} equals \eqn{\beta^{'} x_{it}}. 

More details about the generative process can be found in the vignette.

If \code{corr} is a matrix, then the NORTA method is employed and the correlation matrix corresponds to that of the multivariate normal distribution. Otherwise, it must be a constant on [0,1). For the \code{cloglog} link, it generates latent random vectors from a multivariate Gumbel distribution with correlation parameter \code{corr}. For the \code{logit} link, it generates the latent random vectors as the difference of two indepedent random vectors that follow a multivariate Gumbel distribution with correlation parameter \code{corr}.
}
\value{
\item{Ysim}{the correlated ordinal multinomial responses. Element (\eqn{i},\eqn{t}) represents the realization of \eqn{Y_{it}}.}
\item{correlation}{the latent correlation matrix.}
\item{rlatent}{the underlying latent random vectors.}
}
\references{
Balakrishnan, N. (1992). \emph{Handbook of the logistic distribution.} CRC Press.

Gumbel, E.J. (1961). \emph{Multivariate extremal distributions}. Bulletin de l'Institut International de Statistiques.

Li, S.T. and Hammond, J.L. (1975). Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557-561.  

Touloumis, A., Agresti, A. and Kateri, M. (2013). GEE for multinomial responses using a local odds ratios parameterization. \emph{Biometrics}, \bold{69}, 633-640.
}
\author{
Anestis Touloumis
}
\seealso{
For simulating correlated nominal responses use \code{\link{rmult.bcl}}. For simulating correlated ordinal responses see also \code{\link{rmult.crm}}.
}
\examples{
## Simulating 500 correlated ordinal multinomial responses assuming a cumulative probit
## link model for the marginal univariate probabilities. The sampling scheme is 
## described in Touloumis et al. (2013).

set.seed(12345)
N <- 500
clustersize <- 4
intercepts <- c(-Inf,-1.5,-0.5,0.5,1.5,Inf)
cormat <- toeplitz(c(1,0.85,0.5,0.15))
linpred <- rsmvnorm(N,toeplitz(c(1,rep(0.85,clustersize-1))))
Y <- rmult.clm(clsize=clustersize,lin.pred=linpred,corr=cormat,
               cuts=intercepts,link="probit")$Ysim
data <- cbind(c(t(Y)),c(t(linpred)))
data <- data.frame(data)
data$id <- rep(1:N,each=clustersize)
data$time <- rep(1:clustersize,N)
colnames(data) <- c("y","x","id","time")

## Fit a GEE model to estimate the marginal regression coefficients.
library(multgee)
fitmod <- ordLORgee(y~x,data=data,id=id,repeated=time,link="probit",IM="cholesky")
coef(fitmod)
}
