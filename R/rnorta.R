rnorta <-
function(R=R,cor.matrix=cor.matrix,distr="normal")
{
 if(!is.numeric(R) | R<1) 
    stop("'R' must be greater than or equal to one")
 R <- as.integer(R)
 if(!is.numeric(cor.matrix))  
    stop("'cor.matrix' must be numeric")
 cor.matrix <- as.matrix(cor.matrix)
 if(!isSymmetric(cor.matrix)) 
    stop("'cor.matrix' must be symmetric") 
 if(any(diag(cor.matrix)!=1)) 
    stop("the diagonal elements of 'cor.matrix' must be one")
 if(any(cor.matrix> 1) | any(cor.matrix< -1))
    stop("all the elements of 'cor.matrix' must be on [-1,1]")
 if(any(eigen(cor.matrix,symmetric=TRUE,only.values=TRUE)$values<0))
    stop("'cor.matrix' must be semi-positive definite")
 distrs <- c("normal","logistic","extreme","cauchit")
 if(!is.element(distr,distrs)) 
    stop("'distr' must be either 'normal','logistic','extreme' or 'cauchit'") 
 ans <- rmvnorm(R,sigma=cor.matrix)
 if(distr=="logistic")  ans <- qlogis(pnorm(ans)) 
 if(distr=="extreme")   ans <- qgumbel(pnorm(ans))  
 if(distr=="cauchit")   ans <- qcauchy(pnorm(ans))
 ans
}

