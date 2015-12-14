\name{rsmvnorm}
\alias{rsmvnorm}
\title{
Simulating Continuous Random Vectors from a Multivariate Normal Distribution
}
\description{
Utility function to simulate continuous random vectors from a multivariate normal distribution such that all the marginal distributions are univariate standard normal.
}
\usage{
rsmvnorm(R = R, cor.matrix = cor.matrix)
}
\arguments{
  \item{R}{
integer indicating the sample size.
}
  \item{cor.matrix}{
matrix indicating the correlation matrix of the multivariate normal distribution. 
}
}
\details{
Checks are made to ensure that \code{cor.matrix} is a positive definite correlation matrix. This is determined by the eigenvalues.
}
\value{
Returns \code{R} random vectors of size equal to the row dimension of \code{cor.matrix}.  
}
\author{
Anestis Touloumis
}

