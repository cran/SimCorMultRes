\name{rmult.bcl}
\alias{rmult.bcl}
\title{
Simulating Correlated Nominal Multinomial Responses
}
\description{
This function generates correlated nominal multinomial responses assuming a baseline category logit model for the univariate marginal probabilities.
}
\usage{
rmult.bcl(clsize, ncategories, lin.pred, cor.matrix)
}
\arguments{
  \item{clsize}{
integer indicating the common cluster size.
}
  \item{ncategories}{
integer indicating the number of response categories. 
}
  \item{lin.pred}{
matrix with \code{clsize} x \code{ncategories} columns. See details for more info.
}
  \item{cor.matrix}{
matrix indicating the correlation matrix of the multivariate normal distribution in the NORTA method.
}
}
\details{
\enumerate{
\item Each multinomial response is obtained by utilizing the principle of maximum random utility.
\item The row dimension of \code{lin.pred} equals the sample size. Rows correspond to subjects. Columns \eqn{t*(\code{ncategories})+1,...,2*\code{ncategories}} correspond to the \eqn{1-,...,\code{ncategories}}-th category specific linear predictor at occasion \eqn{t} respectively.
\item The algorithm forces \code{cor.matrix} to respect the local independence assumption between the alternatives. 
}
}
\value{
\item{Ysim}{the correlated nominal multinomial responses.}
\item{correlation}{the correlation matrix of the multivariate normal distribution in the NORTA method.}
\item{rlatent}{the latent random vectors after applying the NORTA method.}
}
\references{
Li, S.T. and Hammond, J.L. (1975). Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. IEEE Transacions on Systems, Man and Cybernetics (5), pp 557-561.

McFadden, D. (1973). Conditional logit analysis of qualitative choice behavior. Institute of Urban and Regional Development, University of California.
}
\author{
Anestis Touloumis
}
\seealso{
\code{\link{rnorta}}, \code{\link{rmult.clm}}
}
\examples{
## Simulating 500 correlated nominal multinomial responses assuming a baseline category 
## logit model for the univariate marginal probabilities.
## In particular we assume that all the category-specific intercepts are equal to 1
## and all the category-specific coefficients are equal to 2.
## Further we sample a cluster specific covariate from a standard normal distribution.

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


## We fit a GEE model to illustrate that the regression coefficients can be estimated
library(multgee)
Y <- rmult.bcl(clsize=3,ncategories=4,lin.pred=linpred,cor.matrix=cormat)$Ysim
data <- cbind(c(t(Y)),c(t(Xmat[,-ncategories])))
data <- data.frame(data)
data$id <- rep(1:N,each=clustersize)
data$time <- rep(1:clustersize,N)
colnames(data) <- c("y","x","id","time")
fitmod <- nomLORgee(y~x,id="id",repeated="time",data=data,add=0.01)

## The GEE estimates of the regression coefficients
coef(fitmod)

}

