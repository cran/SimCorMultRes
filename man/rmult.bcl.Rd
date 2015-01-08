 \name{rmult.bcl}
\alias{rmult.bcl}
\title{
Simulating Correlated Nominal Multinomial Responses Under A Marginal Baseline Category Logit Model
}
\description{
This function simulates correlated nominal multinomial responses assuming that a baseline category logit model holds for the marginal probabilities.
}
\usage{
rmult.bcl(clsize, ncategories, lin.pred, cor.matrix)
}
\arguments{
  \item{clsize}{
integer indicating the common cluster size.
}
  \item{ncategories}{
integer indicating the number of nominal response categories. 
}
  \item{lin.pred}{
matrix with \code{clsize} x \code{ncategories} columns. See details for more info.
}
  \item{cor.matrix}{
matrix indicating the correlation matrix of the multivariate normal distribution in the NORTA method.
}
}
\details{
The assumed marginal baseline category logit model is 
\deqn{log \frac{Pr(Y_{it}=j |x_{it})}{Pr(Y_{it}=I |x_{it})}=(\beta_{0j}-\beta_{0I}) + (\beta^{'}_j-\beta^{'}_I) x_{it}=\beta^{*}_{0j}+ \beta^{*'}_j x_{it}}
where \eqn{Y_{it}} is the \eqn{t}-th multinomial response for subject \eqn{i}, \eqn{x_{it}} is the associated covariates vector, \eqn{\beta_{0j}} is the \eqn{j}-th response category specific intercept and \eqn{\beta_{j}} is the \eqn{j}-th response category specific parameter vector.

The multinomial response \eqn{Y_{it}} is obtained by applying the principle of maximum random utility to the corresponding underlying latent variables.

The row dimension of \code{lin.pred} equals the sample size. Row \eqn{i} corresponds to subject \eqn{i}. Columns \eqn{(t-1)*\code{ncategories}+1,...,t*\code{ncategories}} correspond to the linear predictors \eqn{\beta_{01}+ \beta^{'}_1 x_{it},...,\beta_{0I}+ \beta^{'}_I x_{it}}, respectively, for \eqn{t=1,\ldots,\code{clsize}}. 

More details about the generative process can be found in the vignette.

The algorithm forces \code{cor.matrix} to respect the local independence assumption between the alternatives. 
}
\value{
\item{Ysim}{the correlated nominal multinomial responses. Element (\eqn{i},\eqn{t}) represents the realization of \eqn{Y_{it}}.}
\item{correlation}{the correlation matrix of the multivariate normal distribution in the NORTA method.}
\item{rlatent}{the latent random vectors after applying the NORTA method.}
}
\references{
Li, S.T. and Hammond, J.L. (1975). Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557-561.

McFadden, D. (1973). \emph{Conditional logit analysis of qualitative choice behavior}. Institute of Urban and Regional Development, University of California.
}
\author{
Anestis Touloumis
}
\seealso{
For simulating correlated ordinal responses use \code{\link{rmult.clm}} or \code{\link{rmult.crm}}.
}
\examples{
## Simulating 500 correlated nominal multinomial responses assuming a baseline category 
## logit model for the marginal univariate probabilities. In particular, assume that the
## first (I-1) category-specific intercepts are all equal to 1 and the first (I-1) 
## category-specific coefficients are all equal to 2. For identifiability reasons, set 
## equal to 0 all the reference-specific parameters. Further, sample the cluster specific
## covariates from the standard normal distribution.

set.seed(1)
N <- 500
ncategories <- 4
clustersize <- 3
Xmat <- matrix(rnorm(N),N,ncategories)
betas <- c(1,2,1,2,1,2)
linpred <- matrix(c(betas[c(2,4,6)],0),N,4,byrow=TRUE)*Xmat+
           matrix(c(betas[c(1,3,5)],0),N,4,byrow=TRUE)
linpred <- matrix(linpred,N,ncategories*clustersize)
cormat <- toeplitz(c(1,rep(0,3),rep(c(0.9,0,0,0),2)))


## Fit a GEE model to estimate the regression coefficients of the marginal model,
## i.e., the beta stars.

library(multgee)
Y <- rmult.bcl(clsize=3,ncategories=4,lin.pred=linpred,cor.matrix=cormat)$Ysim
data <- cbind(c(t(Y)),c(t(Xmat[,-ncategories])))
data <- data.frame(data)
data$id <- rep(1:N,each=clustersize)
data$time <- rep(1:clustersize,N)
colnames(data) <- c("y","x","id","time")
fitmod <- nomLORgee(y~x,id=id,repeated=time,data=data,add=0.01)
coef(fitmod)
}
