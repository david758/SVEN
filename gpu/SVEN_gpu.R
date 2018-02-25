#' Train an elastic net using a reduction to a binary SVM as described in Zhou et al.
#' 
#' @param A nxd matrix containing the training data
#' @param A vector of length n with labels 
#' 
SVEN <- function(B,y,t,lambda,usegpu=F) {
  n <- dim(B)[1]
  d <- dim(B)[2]
  
  #Globals used by svm functions
  X <<- cbind(B - (y/t) %*% rep(1,d),B + (y/t) %*% rep(1,d))
  Y <- matrix(c(rep(1,d),rep(-1,d)),nrow=2*d,ncol=1)
  out <- matrix(1,nrow=2*d,ncol=1)
  w <- matrix(0,nrow=n,ncol=1)
  
  if (usegpu) {
    if (!require(gpuR)) {
      stop("gpuR package required.")
    }
    out <- vclMatrix(out)
    w <- vclMatrix(w)
    Y <- vclMatrix(Y)
    X <<- vclMatrix(X)
  }
  
  if (dim(X)[2] < dim(X)[1]) {
    print("here")
    K <<- t(X) %*% X
    beta <- svm(F,Y,lambda,w,out,usegpu)
    beta <- beta * Y
    beta <- t * (beta[1:d,] - beta[(d+1):(2*d),]) / sum(beta)
    return(beta)
  } 
  else {
    print("here2")
    alpha <- svm(T,Y,lambda,w,out,usegpu)
    beta <- alpha / sum(alpha) * t
    return(beta[1:d,] - beta[(d+1):(2*d),])
  }
}


  


  