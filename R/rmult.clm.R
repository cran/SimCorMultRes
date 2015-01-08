rmult.clm <-
function(clsize,lin.pred,corr,cuts,link="probit")
{
 if(!is.numeric(clsize) | clsize < 2)
     stop("'clsize' must be greater than or equal to two")
 clsize <- as.integer(clsize) 
 lin.pred <- as.matrix(lin.pred)
 if(!is.numeric(lin.pred))
     stop("'lin.pred' must be a numeric")
 if(ncol(lin.pred)!=clsize) 
     stop("the matrix 'lin.pred' must have ",clsize,"columns")
 R <- nrow(lin.pred)
 if(!is.vector(cuts))
    stop("'cuts' must be a list or a vector")
 if(!is.numeric(cuts))
    stop("'cuts' must be numeric")
 if(cuts[1]!=-Inf)
   stop("'-Inf' must be the first cutpoint")
 if(cuts[length(cuts)]!=Inf)
   stop("'Inf' must be the last cutpoint")
 if(length(cuts) < 4)
   stop("'cuts' must have at least four elements")
 if(any(diff(cuts)<=0)) 
    stop("'cuts' must be increasing") 
 links <- c("probit","logit","cloglog","cauchit")
 if(!is.element(link,links)) 
   stop("'link' must be either 'probit','logit','cloglog' or'cauchit'") 
 distr <- switch(link,"probit"="normal","logit"="logistic",
                       "cloglog"="extreme","cauchit"="cauchit")
 if(!is.numeric(corr)) 
    stop("'corr' must be numeric")
 if(!is.matrix(corr) & !is.vector(corr))
    stop("'corr' must be matrix or a vector")
 if(is.matrix(corr))
   {
  if(ncol(corr)!=clsize | nrow(corr)!=clsize) 
    stop("'corr' must be a ",clsize,"x",clsize," matrix")
  if(!isSymmetric(corr)) 
    stop("'corr' must be a symmetric matrix") 
  if(any(diag(corr)!=1)) 
    stop("the diagonal elements of 'corr' must be one")
  if(any(corr>1) | any(corr< -1))
    stop("all the elements of 'corr' must be on [-1,1]")
  if(any(eigen(corr,symmetric=TRUE,only.values=TRUE)$values<=0))
    stop("'corr' must be positive definite")
  err <- rnorta(R=R,cor.matrix=corr,distr=distr)
   } else {
 corr <- corr[1]
 if(distr=="normal" | distr=="cauchit") 
    stop("'corr' must be a matrix when 'link'='probit' or 'cauchit'")
 if(any(corr>1) | any(corr<=0))
     stop("'corr' must be on [0,1)") 
 if(distr=="extreme") 
    err <- rmvevd(n=R,dep=sqrt(1-corr),d=clsize) else {
    err <- rmvevd(n=R,dep=sqrt(1-corr),d=clsize)-
           rmvevd(n=R,dep=sqrt(1-corr),d=clsize)
         }
 corr <- toeplitz(c(1,rep(corr,clsize-1)))
 }
 U <- if(distr=="extreme") lin.pred+err else -lin.pred+err
 Ysim <- matrix(cut(U,cuts,labels=FALSE),R,clsize)
 if(distr=="extreme") Ysim <- max(Ysim)-Ysim+1
 list(Ysim=Ysim,correlation=corr,rlatent=err)
 }