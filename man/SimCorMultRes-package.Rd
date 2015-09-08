\name{SimCorMultRes-package}
\alias{SimCorMultRes-package}
\alias{SimCorMultRes}
\docType{package}
\title{
Simulating Correlated Multinomial Responses
}
\description{
Functions to simulate correlated nominal or ordinal multinomial responses (three or more response categories) under a given marginal model specification for the univariate probabilities.
}
\details{
The simulated correlated multinomial responses are drawn as realizations of a latent underlying regression model for continuous random vectors with the correlation structure expressed in terms of the latent random vectors and not of the correlated multinomial responses. 

For an ordinal response scale, the multinomial responses are simulated either under a cumulative link model by chopping off the corresponding latent variables (McCullagh, 1980), or under a marginal continuation ratio model by utilizing the latent variable representation described in Tutz (1991).

For a nominal response scale, the multinomial responses are simulated under a marginal baseline category logit model by utilizing the principle of maximum random utility (McFadden, 1973).

The package aims to ease the simulation of correlated multinomial responses when the interest lies on estimating the marginal regression coefficient parameters, e.g., by fitting a local odds ratios GEE model (Touloumis et al., 2013).
}
\author{
Anestis Touloumis <Anestis.Touloumis@brighton.ac.uk>

Maintainer: Anestis Touloumis <Anestis.Touloumis@brighton.ac.uk>
}
\references{
McCullagh, P. (1980). Regression models for ordinal data. \emph{Journal of the Royal Statistical Society B}, \bold{42}, 109-142.

McFadden, D. (1973). \emph{Conditional logit analysis of qualitative choice behavior}. Institute of Urban and Regional Development, University of California.

Touloumis, A., Agresti, A. and Kateri, M. (2013). GEE for multinomial responses using a local odds ratios parameterization. \emph{Biometrics}, \bold{69}, 633-640.

Tutz, G. (1991). Sequential models in categorical regression, \emph{Computational Statistics & Data Analysis}, \bold{11}, 275-295.
}
\keyword{ package }
