\name{rnorta}
\alias{rnorta}
\title{
Simulating Continuous Random Vectors using the NORTA Method  
}
\description{
Utility function to simulate continuous random vectors with predefined marginal distributions via the NORTA method.  
}
\usage{
rnorta(R = R, cor.matrix = cor.matrix, distr = "normal")
}
\arguments{
  \item{R}{
integer indicating the sample size.
}
  \item{cor.matrix}{
matrix indicating the correlation matrix of the multivariate normal distribution employed in the NORTA method.
}
  \item{distr}{
character string indicating the desired univariate marginal distributions. Options include \code{"normal"}, \code{"logistic"}, \code{"extreme"} or \code{"cauchy"}.
}
}
\details{
Checks are made to ensure that \code{cor.matrix} is a positive definite correlation matrix. This is determined by the eigenvalues.
}
\value{
Returns \code{R} random vectors of size equal to the row dimension of \code{cor.matrix} and marginal distributions \code{distr}.  
}
\references{
Cario, M. C. and Nelson, B. L. (1997)  \emph{Modeling and generating random vectors with arbitrary marginal distributions and correlation matrix}. Technical Report, Department of Industrial Engineering and Management Sciences, Northwestern University, Evanston, Illinois.

Li, S. T. and Hammond, J. L. (1975)  Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557--561.
}
\author{
Anestis Touloumis
}
\examples{
## An example with marginal distributions the standard logistic.
set.seed(1)
## The latent correlation matrix
latent.cor <- toeplitz(c(1, rep(0.8, 3)))
logit.sim <- rnorta(R = 1000, cor.matrix = latent.cor, distr = "logistic")

## The following code illustrates the NORTA method.
set.seed(1)
norm.sim <- rsmvnorm(1000, latent.cor)
all(logit.sim==qlogis(pnorm(norm.sim)))
}
