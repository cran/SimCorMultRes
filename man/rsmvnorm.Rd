\name{rsmvnorm}
\alias{rsmvnorm}
\title{
Simulating Continuous Random Vectors From A Multivariate Normal Distribution
}
\description{
This function simulates continuous random vectors from a multivariate normal distribution such that the marginal distributions are standard normal.
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

