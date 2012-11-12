\name{rmult.clm}
\alias{rmult.clm}
\title{
Simulating Correlated Ordinal Multinomial Responses
}
\description{
This function generates correlated ordinal multinomial responses assuming a cumulative link model for the univariate marginal probabilities.
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
character string indicating the link function of the cumulative link model. Options include \code{"probit"}, \code{"logit"}, \code{"clogclog"} or \code{"cauchit"}.
}
}
\details{
\enumerate{
\item Each multinomial response takes the value \eqn{j} if and only if the corresponding latent random variable is on \eqn{(cuts_{j-1},cuts_{j}]}.
\item The row dimension of \code{lin.pred} equals the sample size. 
\item \code{lin.pred} contains the linear predictor of the cumulative link model after subtracting the intercepts.
\item If \code{corr} is a matrix, then the NORTA method is employed and the correlation matrix corresponds to that of the multivariate normal distribution. Otherwise, it must be a constant on [0,1). For the \code{cloglog} link, it generates latent random vectors from a multivariate Gumbel distribution with correlation parameter \code{corr}. For the \code{logit} link, it generates the latent random vectors as the difference of two indepedent random vectors that follow a multivariate Gumbel distribution with correlation parameter \code{corr}.
}
}
\value{
\item{Ysim}{the correlated ordinal multinomial responses.}
\item{correlation}{the latent correlation matrix.}
\item{rlatent}{the underlying latent random vectors.}
}
\references{
Li, S.T. and Hammond, J.L. (1975). Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. IEEE Transacions on Systems, Man and Cybernetics (5), pp 557-561.  

Balakrishnan, N. (1992). Handbook of the logistic distribution.

Gumbel, E.J. (1961). Multivariate extremal distributions. Bulletin de l'Institut International de Statistiques.

Touloumis, A., Agresti, A. and Kateri, M. (2012). GEE for multinomial responses using a local odds ratios parameterization. \emph{Submitted}.
}
\author{
Anestis Touloumis
}
\seealso{
\code{\link{rnorta}}, \code{\link{rmult.bcl}}
}
\examples{
## Simulating 500 correlated ordinal multinomial responses assuming a probit cumulative
## link model for the univariate marginal probabilities. The model is described in Touloumis,
## Agresti and Kateri (2012) with a toeplitz underlying correlation structure and a 
## time-varying covariates design. 

set.seed(12345)
N <- 500
clustersize <- 4
intercepts <- c(-Inf,-1.5,-0.5,0.5,1.5,Inf)
cormat <- toeplitz(c(1,0.85,0.5,0.15))
linpred <- rmvnorm(N,sigma=toeplitz(c(1,rep(0.85,clustersize-1))))

Y <- rmult.clm(clsize=clustersize,lin.pred=linpred,corr=cormat,
               cuts=intercepts,link="probit")$Ysim
data <- cbind(c(t(Y)),c(t(linpred)))
data <- data.frame(data)
data$id <- rep(1:N,each=clustersize)
data$time <- rep(1:clustersize,N)
colnames(data) <- c("y","x","id","time")

## We fit a GEE model to illustrate that the regression coefficients can be estimated
library(multgee)
fitmod <- ordLORgee(y~x,data=data,id="id",repeated="time",link="probit",IM="cholesky")

## The GEE estimates of the regression coefficients
coef(fitmod)

}
