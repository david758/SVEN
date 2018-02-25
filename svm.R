svm <- function(linear,Y,lambda,w,out) {

  opt.iter_max_Newton <- 50
  opt.lin_cg <- FALSE
  opt.prec <- 1e-6
  opt.cg_it <- 5

  norm <- function(x) {sqrt(sum(x^2))}
  
  obj_fun_linear <- function(w,Y,lambda,out) {
    out <- pmax(out,0)
    obj <- sum(out^2) / 2 + lambda * (w %*% w)[1] / 2
    grad <- lambda * w - as.vector(X %*% (out * Y))
    sv <- which(out > 0)
    return(list(obj=obj,grad=grad,sv=sv))
  }
  
  line_search_linear <- function(w,d,out,Y,lambda,tol) {
    t = 0
    Xd <- as.vector(d %*% X)
    wd <- lambda * (d %*% w)[1]
    dd <- lambda * (d %*% d)[1]
    while (TRUE) {
      out2 <- out - t * (Y * Xd)
      sv <- which(out2 > 0)
      g <- wd + t * dd - ((out2 * Y)[sv] %*% Xd[sv])[1]
      h <- dd + (Xd[sv] %*% Xd[sv])[1]
      t <- t - g/h
      if (g^2/h < tol) {
        break
      }
    }
    return(list(t=t,out=out2))
  }
  
  #Requires: X, Y
  primal_svm_linear <- function(Y,lambda,w,out) {
    d <- dim(X)[1]
    n <- dim(X)[2]
    
    iter <- 0
    diagD <- diag(d)
    
    while(TRUE) {
      iter = iter + 1
      if (iter > opt.iter_max_Newton) {
        print("Maximum number of Newton steps reached. Try larger lambda.")
        break
      }
      
      objFunOut <- obj_fun_linear(w,Y,lambda,out)
      obj <- objFunOut[["obj"]]
      grad <- objFunOut[["grad"]]
      sv <- objFunOut[["sv"]]

      if (length(sv) == n) {
        Xsv <- X
      }
      else {
        Xsv <- X[,sv]
      }
        
      hess <- lambda * diagD + Xsv %*% t(Xsv)
      step <- as.vector(-solve(hess) %*% grad)
      
      lsOut <- line_search_linear(w,step,out,Y,lambda,1e-12)
      out <- lsOut[["out"]]
      t <- lsOut[["t"]]
      
      w <- w + t * step
      
      if (-(step %*% grad)[1] < opt.prec * obj) {
        break
      }
    }
    
    out_l <- out / norm(out)
    w_l <- w / norm(out)
    alpha <- pmax(out,0)
    
    return(list(out_l=out_l,obj=obj,w_l=w_l,alpha=alpha))
  }
  
  primal_svm_linear_cg <- function(Y,lambda,w,out) {
    d <- dim(X)[1]
    n <- dim(X)[2]
    
    stop <- 0
    iter <- 0
    go <- as.vector(X %*% Y)
    
    s <- go
    while (TRUE) {
      iter <- iter + 1
      if (iter > opt.cg_it) {
        break
      }
      
      lsOut <- line_search_linear(w,s,out,Y,lambda,1e-10)
      out <- lsOut[["out"]]
      t <- lsOut[["t"]]
      w = w + t * s
      
      objFunOut <- obj_fun_linear(w,Y,lambda,out)
      obj <- objFunOut[["obj"]]
      gn <- -objFunOut[["grad"]]
      sv <- objFunOut[["sv"]]
      
      if (t * (s %*% go)[1] < opt.prec * obj) {
        stop <- 1
        break
      }
      
      be <- ((gn %*% gn)[1] - 0 * (gn %*% go)[1]) / (go %*% go)[1]
      s <- be * s + gn
      go <- gn
    }
    
    out_l <- out / norm(out)
    w_l <- w / norm(out)
    alpha <- pmax(out,0)
    
    return(list(w_l=w_l,obj=obj,alpha=alpha,out_l=out_l,stop=stop))
  }
  
  line_search_nonlinear <- function(step,Kb,b,Y,lambda,fullstep) {
    training <- which(Y != 0)
    act <- which(step != 0)
    Ks <- as.vector(K[training,training[act]] %*% step[act])
    
    Kss <- (step[act] %*% Ks[act])[1]
    Kbs <- (step[act] %*% Kb[act])[1]
    
    t <- 0
    Y <- Y[training]
    out <- 1 - Y * (Kb+Ks)
    sv <- which(out > 0)
    obj1 <- (lambda * (2 * Kbs + Kss) + sum(out[sv]^2)) / 2
    
    while (TRUE) {
      out <- 1 - Y * (Kb + t * Ks)
      sv <- which(out > 0)
      obj <- (lambda * (2 * t * Kbs + t^2 * Kss) + sum(out[sv]^2)) / 2
      g <- lambda * (Kbs + t * Kss) - (Ks[sv] %*% (Y[sv] * out[sv]))[1]
      
      if (fullstep & t == 0 & obj-obj1 > -0.2*g) {
        t <- 1
        break
      }
      
      h <- lambda * Kss + norm(Ks[sv])^2
      t <- t - g/h
      
      if (g^2 / h < 1e-10) {
        break
      }
    }
    
    return(list(t=t,Kb=Kb + t * Ks))
  }

  primal_svm_nonlinear <- function(Y,lambda) {
      training <- which(Y != 0)
      n <- length(training)
      
      if (n >= 1100) {
        perm <- sample(1:n,n)
        ind <- training[perm[1:floor(0.75*n)]]
        Y2 <- Y
        Y2[ind] <- 0
        beta <- primal_svm_nonlinear(Y2,lambda)
        sv <- which(beta != 0)
        Kb <- K[training,sv] %*% beta[sv]
      }
      else {
        sv <- training
        beta <- rep(0,length(Y))
        Kb <- rep(0,n)
      }
      
      iter <- 0
      
      while (TRUE) {
        old_sv <- sv
        out <- 1 - Y[training] * Kb
        sv <- training[out > 0]
        obj <- lambda * (beta[training] %*% Kb)[1] + sum(pmax(out,0)^2) / 2

        iter <- iter + 1
        if (iter > 1 & length(unique(c(setdiff(sv,old_sv), setdiff(old_sv,sv)))) == 0) {
          break
        } 
        if (iter > opt.iter_max_Newton) {
          print("Maximum number of Newton steps reached. Try larger lambda")
          break
        }
        
        H <- K[sv,sv] + lambda * diag(length(sv))
        beta_new <- rep(0,length(Y))
        beta_new[sv] <- as.vector(solve(H) %*% Y[sv])

        step <- beta_new - beta

        lsOut <- line_search_nonlinear(step[training],Kb,beta[length(beta)],Y,lambda,TRUE)
        t <- lsOut[["t"]]
        Kb <- lsOut[["Kb"]]
        beta <- beta + t * step
      }   
      return(beta)
  }
  
  if (linear) {
    pslcgOut <- primal_svm_linear_cg(Y,lambda,w,out)
    w <- pslcgOut[["w_l"]]
    out <- pslcgOut[["out_l"]]
    pslOut <- primal_svm_linear(Y,lambda,w,out)
    return(pslOut[["alpha"]])
  } 
  else {
    sol <- primal_svm_nonlinear(Y,lambda)
    return(sol)
  }
}



