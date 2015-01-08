### R code from vignette source 'SCMRvin.Rnw'

###################################################
### code chunk number 1: SCMRvin.Rnw:40-52
###################################################
library("SimCorMultRes")
set.seed(1)
N <- 500
ncategories <- 4
clustersize <- 3
Xmat <- matrix(rnorm(N),N,ncategories)
betas <- c(1,2,3,4,5,6)
linpred <- matrix(c(betas[c(2,4,6)],0),N,4,byrow=TRUE)*Xmat+
           matrix(c(betas[c(1,3,5)],0),N,4,byrow=TRUE)
linpred <- matrix(linpred,N,ncategories*clustersize)
cormat <- diag(1,12)
Y <- rmult.bcl(clsize=3,ncategories=4,lin.pred=linpred,cor.matrix=cormat)


###################################################
### code chunk number 2: SCMRvin.Rnw:55-56
###################################################
head(Y$Ysim)


###################################################
### code chunk number 3: SCMRvin.Rnw:87-96
###################################################
set.seed(1)
N <- 500
clustersize <- 4
intercepts <- c(-Inf,-1.5,-0.5,0.5,1.5,Inf)
cormat <- toeplitz(c(1,0.85,0.5,0.15))
x <- rnorm(N)
linpred <- matrix(rep(x,clustersize),N,clustersize,byrow=TRUE)
Y <- rmult.clm(clsize=clustersize,lin.pred=linpred,corr=cormat,
               cuts=intercepts,link="probit")


###################################################
### code chunk number 4: SCMRvin.Rnw:99-100
###################################################
head(Y$Ysim)


###################################################
### code chunk number 5: SCMRvin.Rnw:124-133
###################################################
set.seed(1)
N <- 500
clustersize <- 4
intercepts <- c(-Inf,-1.5,-0.5,0.5,1.5,Inf)
cormat <- diag(1,16)
x <- rnorm(N)
linpred <- matrix(rep(x,clustersize),N,clustersize,byrow=TRUE)
Y <- rmult.crm(clsize=clustersize,lin.pred=linpred,cor.matrix=cormat,
               cuts=intercepts,link="probit")


###################################################
### code chunk number 6: SCMRvin.Rnw:136-137
###################################################
head(Y$Ysim)


