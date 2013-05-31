rmult.bcl <-
function(clsize,ncategories,lin.pred,cor.matrix)
{
 if(!is.numeric(clsize) | clsize < 2)
     stop("'clsize' must be greater than or equal to two")
 clsize <- as.integer(clsize) 
 if(!is.numeric(ncategories) | ncategories < 3 )
     stop("'ncategories' must be greater than or equal to three")
 ncategories <- as.integer(ncategories) 
 dims.lp <- clsize*ncategories
 if(!is.numeric(lin.pred))
     stop("'lin.pred' must be a numeric matrix")
 lin.pred <- as.matrix(lin.pred)
 if(ncol(lin.pred)!= dims.lp) 
     stop("'lin.pred' must have ",dims.lp," columns")
 R <- nrow(lin.pred)

 if(!is.numeric(cor.matrix)) 
     stop("'cor.matrix' must be numeric")
 cor.matrix <- as.matrix(cor.matrix)
 if(ncol(cor.matrix)!=dims.lp) 
    stop("'cor.matrix' must be a ",dims.lp,"x",dims.lp," matrix")
 if(!isSymmetric(cor.matrix))
     stop("'cor.matrix' must be a symmetric matrix")
 for(i in 1:clsize){
 diag.index <- 1:ncategories+(i-1)*ncategories
 cor.matrix[diag.index,diag.index]<- diag(1,ncategories)
  }
 if(any(cor.matrix>1) | any(cor.matrix< -1))
     stop("all the elements of 'cor.matrix' must be on [-1,1]")
 if(any(eigen(cor.matrix,symmetric=TRUE,only.values=TRUE)$values<=0))
     stop("'cor.matrix' must respect the local independence of the alternatives and it must be a positive definite matrix")
 err <- rnorta(R,cor.matrix=cor.matrix,distr="extreme")
 U <- lin.pred + err
 U <- matrix(as.vector(t(U)),nrow=clsize*R,ncol=ncategories,TRUE)
 Ysim <- apply(U,1,which.max)
 Ysim <- matrix(Ysim,ncol=clsize,byrow=TRUE)
 if(!is.matrix(cor.matrix)) 
    cor.matrix <- toeplitz(
                  c(1,rep(0,ncategories-1),
                  rep(c(cor.matrix,rep(0,ncategories-1)),clsize-1)))
 list(Ysim=Ysim,correlation=cor.matrix,rlatent=err)
 }