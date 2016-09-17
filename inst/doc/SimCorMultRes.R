### R code from vignette source 'SimCorMultRes.Rnw'

###################################################
### code chunk number 1: SimCorMultRes.Rnw:46-58
###################################################
library(SimCorMultRes)
set.seed(1)
ncategories <- 4
N <- 500
clsize <- 3
betas <- c(1, 2, 1.5, 3, 4, 3.5, 5, 6, 5.5, 0, 0, 0)
x1 <- rep(rnorm(N), each = clsize)
x2 <- rnorm(N * clsize)
xdata <- data.frame(x1, x2)
cor.matrix <- diag(1, 12)
CorNorRes <- rmult.bcl(clsize = clsize, ncategories = ncategories, betas = betas, 
    xformula = ~x1 + x2, xdata = xdata, cor.matrix = cor.matrix)


###################################################
### code chunk number 2: SimCorMultRes.Rnw:61-62
###################################################
head(CorNorRes$Ysim)


###################################################
### code chunk number 3: SimCorMultRes.Rnw:65-70
###################################################
library(evd)
rlatent <- rmvevd(n = N, dep = 1, model = "log", d = clsize * ncategories)
CorNorRes <- rmult.bcl(clsize = clsize, ncategories = ncategories, betas = betas, 
    xformula = ~x1 + x2, xdata = xdata, rlatent = rlatent)
head(CorNorRes$Ysim)


###################################################
### code chunk number 4: SimCorMultRes.Rnw:107-116
###################################################
set.seed(12345)
N <- 500
clsize <- 4
intercepts <- c(-1.5, -0.5, 0.5, 1.5)
betas <- matrix(c(1, 2, 3, 4), 4, 1)
x <- rep(rnorm(N), each = clsize)
cor.matrix <- toeplitz(c(1, 0.85, 0.5, 0.15))
CorOrdRes <- rmult.clm(clsize = clsize, intercepts = intercepts, betas = betas, 
    xformula = ~x, cor.matrix = cor.matrix, link = "probit")


###################################################
### code chunk number 5: SimCorMultRes.Rnw:119-120
###################################################
head(CorOrdRes$Ysim)


###################################################
### code chunk number 6: SimCorMultRes.Rnw:147-156
###################################################
set.seed(1)
N <- 500
clsize <- 4
intercepts <- c(-1.5, -0.5, 0.5, 1.5)
cor.matrix <- diag(1, 16)
x <- rnorm(N * clsize)
CorOrdRes <- rmult.crm(clsize = clsize, intercepts = intercepts, betas = 1, 
    xformula = ~x, cor.matrix = cor.matrix, link = "probit")



###################################################
### code chunk number 7: SimCorMultRes.Rnw:159-160
###################################################
head(CorOrdRes$Ysim)


###################################################
### code chunk number 8: SimCorMultRes.Rnw:194-203
###################################################
set.seed(123)
N <- 5000
clsize <- 4
intercepts <- 0
betas <- 0.2
cor.matrix <- toeplitz(c(1, 0.9, 0.9, 0.9))
x <- rep(rnorm(N), each = clsize)
CorBinRes <- rbin(clsize = clsize, intercepts = intercepts, betas = betas, 
    xformula = ~x, cor.matrix = cor.matrix, link = "probit")


###################################################
### code chunk number 9: SimCorMultRes.Rnw:206-209
###################################################
library(gee)
binGEEmod <- gee(y ~ x, family = binomial("probit"), id = id, data = CorBinRes$simdata)
summary(binGEEmod)


###################################################
### code chunk number 10: SimCorMultRes.Rnw:212-221
###################################################
set.seed(123)
library(evd)
rlatent1 <- rmvevd(N, dep = sqrt(1 - 0.9), model = "log", d = clsize)
rlatent2 <- rmvevd(N, dep = sqrt(1 - 0.9), model = "log", d = clsize)
rlatent <- rlatent1 - rlatent2
CorBinRes <- rbin(clsize = clsize, intercepts = intercepts, betas = betas, 
    xformula = ~x, rlatent = rlatent)
binGEEmod <- gee(y ~ x, family = binomial("logit"), id = id, data = CorBinRes$simdata)
summary(binGEEmod)


###################################################
### code chunk number 11: SimCorMultRes.Rnw:224-225
###################################################
citation("SimCorMultRes")


