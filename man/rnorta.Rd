\name{rnorta}
\alias{rnorta}
\title{
Simulating Continuous Random Vectors Using The NORTA Method
}
\description{
This function simulates continuous random vectors with prescribed univariate marginal distributions using the NORTA method.
}
\usage{
rnorta(R = R, cor.matrix = cor.matrix, distr = "normal")
}
\arguments{
  \item{R}{
integer indicating the sample size.
}
  \item{cor.matrix}{
matrix indicating the correlation matrix of the multivariate normal distribution in the NORTA method.
}
  \item{distr}{
character string indicating the desired univariate marginal distributions. Options include \code{"normal"}, \code{"logistic"}, \code{"extreme"} or \code{"cauchy"}.
}
}
\details{
Checks are made to ensure that \code{cor.matrix} is a positive definite correlation matrix. This is determined by the eigenvalues.
}
\value{
Returns \code{R} random vectors of size equal to the row dimension of \code{cor.matrix} and univariate marginal distributions \code{distr}.  
}
\references{
Li, S.T. and Hammond, J.L. (1975). Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557-561.
}
\author{
Anestis Touloumis
}
\examples{
## An example with marginal standard logistic distributions.
set.seed(1)
## The latent correlation matrix
latent.cor <- toeplitz(c(1,rep(0.8,3)))
logit.sim <- 
rnorta(R = 1000, cor.matrix = latent.cor, distr = "logistic")

## The following illustrates how the NORTA method works.
set.seed(1)
norm.sim <- rsmvnorm(1000,latent.cor)
all(logit.sim==qlogis(pnorm(norm.sim)))
}
