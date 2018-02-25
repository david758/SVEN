#' Trains an elastic net model using the method described in Zhou et al. 2014
#' 
#' @param B A nxd matrix of training data
#' @param y A n length vector of labels
#' @param t L1 regularization coefficient 
#' @param lambda L2 regularization coefficient
#' 
#' @return Elastic net weight vector
#' 
#' @example B <- matrix(rnorm(100,100),100,100)
#' @example y <- rnorm(100)
#' @example SVEN(B,y,0.1,0.5)
#' 
SVEN <- function(B,y,t,lambda) {
  source("svm.R")
  
  n <- dim(B)[1]
  d <- dim(B)[2]
  
  X <<- cbind(B - replicate(d,y/t),B + replicate(d,y/t))
  Y <- c(rep(1,d),rep(-1,d))
  out <- rep(1,2*d)
  w <- rep(0,n)
  
  if (dim(X)[2] < dim(X)[1]) {
    K <<- t(X) %*% X
    beta <- svm(F,Y,lambda,w,out)
    beta <- beta * Y
    beta <- t * (beta[1:d] - beta[(d+1):(2*d)]) / sum(beta)
    return(beta)
  } 
  else {
    alpha <- svm(T,Y,lambda,w,out)
    beta <- alpha / sum(alpha) * t
    return(beta[1:d] - beta[(d+1):(2*d)])
  }
}


  


  