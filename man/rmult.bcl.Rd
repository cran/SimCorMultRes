 \name{rmult.bcl}
\alias{rmult.bcl}
\title{
Simulating Correlated Nominal Multinomial Responses under a Marginal Baseline Category Logit Model
}
\description{
Simulates correlated nominal multinomial responses assuming that a baseline category logit model holds for the marginal probabilities.
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
\deqn{log \frac{Pr(Y_{it}=j |x_{it})}{Pr(Y_{it}=I |x_{it})}=(\beta_{t0j}-\beta_{t0I}) + (\beta^{'}_{tj}-\beta^{'}_{tI}) x_{it}=\beta^{*}_{t0j}+ \beta^{*'}_{tj} x_{it}}
For subject \eqn{i}, \eqn{Y_{it}} is the \eqn{t}-th multinomial response and \eqn{x_{it}} is the associated covariates vector. Finally, \eqn{\beta_{t0j}} and \eqn{\beta^{'}_{tj}} is the \eqn{j}-th category specific intercept and regression parameter vector corresponding to the \eqn{t}-th measurement occasion respectively.

The multinomial response \eqn{Y_{it}} is obtained by extending the principle of maximum random utility (McFadden, 1974) to correlated latent continuous variables. See vignette for more info.

The sample size is equal to the number of the rows in \code{lin.pred}. Row \eqn{i} corresponds to subject \eqn{i}. Columns \eqn{(t-1)*\code{ncategories}+1,...,t*\code{ncategories}} correspond to the linear predictors \eqn{\beta_{t01}+ \beta^{'}_{t1} x_{it},...,\beta_{t0I}+ \beta^{'}_{tI} x_{it}}, respectively, for \eqn{t=1,\ldots,\code{clsize}}. 

The algorithm forces \code{cor.matrix} to respect the local independence assumption between the alternatives. 
}
\value{
\item{Ysim}{the correlated nominal multinomial responses. Element (\eqn{i},\eqn{t}) represents the realization of \eqn{Y_{it}}.}
\item{correlation}{the correlation matrix of the multivariate normal distribution in the NORTA method.}
\item{rlatent}{the latent random vectors.}
}
\references{
Cario, M. C. and Nelson, B. L. (1997)  \emph{Modeling and generating random vectors with arbitrary marginal distributions and correlation matrix}. Technical Report, Department of Industrial Engineering and Management Sciences, Northwestern University, Evanston, Illinois.

Li, S. T. and Hammond, J. L. (1975)  Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557--561.   

McFadden, D. (1974) \emph{Conditional logit analysis of qualitative choice behavior}. New York: Academic Press, 105--142.

Touloumis, A., Agresti, A. and Kateri, M. (2013) GEE for multinomial responses using a local odds ratios parameterization. \emph{Biometrics}, \bold{69}, 633--640.
}
\author{
Anestis Touloumis
}
\seealso{
For simulating correlated ordinal multinomial responses use \code{rmult.clm} or \code{rmult.crm}, and for simulating correlated binary responses use \code{rbin}.
}
\examples{
## See vignette for the sampling scheme.
set.seed(1)
N <- 500
ncategories <- 4
clsize <- 3
Xmat <- matrix(rnorm(N), N, ncategories)
betas <- c(1, 2, 3, 4, 5, 6)
lin.pred <- matrix(c(betas[c(2, 4, 6)], 0), N, 4, byrow=TRUE) * Xmat +
           matrix(c(betas[c(1, 3, 5)], 0), N, 4, byrow=TRUE)
lin.pred <- matrix(lin.pred, N, ncategories*clsize)
cor.matrix <- diag(1, 12)

## Fit a GEE model (Touloumis et al., 2013) to estimate the regression coefficients,
## i.e., the beta stars. 
library(multgee)
Y <- rmult.bcl(clsize, ncategories, lin.pred, cor.matrix)$Ysim
data <- cbind(c(t(Y)), c(t(Xmat[,-ncategories])))
data <- data.frame(data)
data$id <- rep(1:N, each=clsize)
data$time <- rep(1:clsize, N)
colnames(data) <- c("y", "x", "id", "time")
fitmod <- nomLORgee(y~x, id=id, repeated=time, data=data, add=0.01)
coef(fitmod)
}
