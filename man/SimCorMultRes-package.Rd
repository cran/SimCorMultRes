\name{SimCorMultRes-package}
\alias{SimCorMultRes-package}
\alias{SimCorMultRes}
\docType{package}
\title{
Simulating Correlated Multinomial Responses
}
\description{
Functions to simulate correlated nominal or ordinal multinomial responses (three or more response categories) and correlated binary responses subject to  a marginal model specification.
}
\details{
The simulated correlated binary or multinomial responses are drawn as realizations of a latent regression model for continuous random vectors with the correlation structure expressed in terms of the latent random vectors and not of the correlated multinomial responses. 

For an ordinal response scale, the multinomial variables are simulated either under a marginal cumulative link model (\code{rmult.clm}) or under a marginal continuation ratio model (\code{rmult.crm}).

For a nominal response scale, the multinomial responses are simulated under a marginal baseline category logit model (\code{rmult.bcl}).

Correlated binary responses are simulated using the function \code{rbin}. 

See vignette for more details about the threshold approaches and the implemented marginal models.
}
\author{
Anestis Touloumis <A.Touloumis@brighton.ac.uk>

Maintainer: Anestis Touloumis <A.Touloumis@brighton.ac.uk>
}
\references{
Cario, M. C. and Nelson, B. L. (1997)  \emph{Modeling and generating random vectors with arbitrary marginal distributions and correlation matrix}. Technical Report, Department of Industrial Engineering and Management Sciences, Northwestern University, Evanston, Illinois.

Li, S. T. and Hammond, J. L. (1975)  Generation of pseudorandom numbers with specified univariate distributions and correlation coefficients. \emph{IEEE Transacions on Systems, Man and Cybernetics}, \bold{5}, 557--561.

McCullagh, P. (1980) Regression models for ordinal data. \emph{Journal of the Royal Statistical Society B}, \bold{42}, 109--142.

McFadden, D. (1974) \emph{Conditional logit analysis of qualitative choice behavior}. New York: Academic Press, 105--142.

Touloumis, A., Agresti, A. and Kateri, M. (2013) GEE for multinomial responses using a local odds ratios parameterization. \emph{Biometrics}, \bold{69}, 633--640.

Tutz, G. (1991) Sequential models in categorical regression. \emph{Computational Statistics & Data Analysis}, \bold{11}, 275--295.
}
\keyword{ package }
