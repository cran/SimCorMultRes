\name{SimCorMultRes-package}
\alias{SimCorMultRes-package}
\alias{SimCorMultRes}
\docType{package}
\title{
Simulating Correlated Multinomial Responses
}
\description{
Functions to simulate correlated nominal or ordinal multinomial responses assuming a cumulative link model or a baseline category logit model for the univariate marginal probabilities.
}
\details{
The simulated correlated multinomial responses are drawn as realizations of an underlying regression model for continuous random vectors. The correlation structure is expressed in terms of the latent random vectors. 
For ordinal response categories, the multinomial responses are produced by chopping off the corresponding latent variables (McCullagh, 1980), and for nominal response categories by utilizing the principle of maximum random utility (McFadden, 1973).
The functions are suitable to generate correlated ordinal or nominal multinomial responses when a cumulative link model or a baseline category logit model, respectively, holds for the univariate marginal probabilities.

The package can ease the simulation of correlated multinomial responses when the interest lies on estimating the marginal regression coefficient parameters, e.g. using a GEE model. An example is given in Touloumis, Agresti and Kateri (2012).
}
\author{
Anestis Touloumis <anestis@ebi.ac.uk>

Maintainer: Anestis Touloumis <anestis@ebi.ac.uk>
}
\references{
McCullagh, P. (1980). Regression models for ordinal data. Journal of the Royal Statistical Society B, pp 109-142.

McFadden, D. (1973). Conditional logit analysis of qualitative choice behavior. Institute of Urban and Regional Development, University of California.

Touloumis, A., Agresti, A. and Kateri, M. (2012). GEE for multinomial responses using a local odds ratios parameterization. \emph{Submitted}.
}
\keyword{ package }
