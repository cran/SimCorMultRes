\name{rmult.crm}
\alias{rmult.crm}
\title{
Simulating Correlated Ordinal Multinomial Responses Under A Marginal Continuation Ratio Model
}
\description{
This function simulates correlated ordinal multinomial responses assuming that a continuation ratio model holds for the marginal probabilities.
}
\usage{
rmult.crm(clsize, lin.pred, cor.matrix, cuts, link = "probit")
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
  \item{cuts}{
vector containing the intercepts of the continuation ratio model, and \code{-Inf} and \code{Inf} as first and last elements respectively.
}
  \item{link}{
character string indicating the link function of the continuation ratio model. Options include \code{"probit"}, \code{"logit"}, \code{"cloglog"} or \code{"cauchit"}.
}
}
\details{
The assumed marginal cumulative link model is  
\deqn{Pr(Y_{it}=j |Y_{it} \ge j,x_{it})=F(\beta_{0j} +\beta^{'} x_{it})}
where \eqn{Y_{it}} is the \eqn{t}-th multinomial response for subject \eqn{i}, \eqn{x_{it}} is the associated covariates vector , \eqn{F} is the cumulative distribution function determined by \code{link}, \eqn{\beta_{0j}} is the \eqn{j}-th response category specific intercept and \eqn{\beta} is the marginal regression parameter vector excluding intercepts.

The multinomial responses \eqn{Y_{it}} are determined by applying the latent variable threshold approach of Tutz (1991) to the corresponding latent variables.

The row dimension of \code{lin.pred} equals the sample size. Row \eqn{i} corresponds to subject \eqn{i}. The (\eqn{i,t})-th element of \code{lin.pred} equals \eqn{\beta^{'} x_{it}}. 

More details about the generative process can be found in the vignette.
}
\value{
\item{Ysim}{the correlated ordinal multinomial responses. Element (\eqn{i},\eqn{t}) represents the realization of \eqn{Y_{it}}.}
\item{correlation}{the correlation matrix of the multivariate normal distribution in the NORTA method.}
\item{rlatent}{the underlying latent random vectors.}
}
\references{
Li, S.T. and Hammond, J.L. (1975). Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557-561.  

Tutz, G. (1991). Sequential models in categorical regression, \emph{Computational Statistics & Data Analysis}, \bold{11}, 275-295.
}
\author{
Anestis Touloumis
}
\seealso{
For simulating correlated nominal responses use \code{\link{rmult.bcl}}. For simulating correlated ordinal responses see also \code{\link{rmult.clm}}.
}

