pglasso_mod <- function(y, X=NULL, fe_mat_indices_list, fe_mat_indices, D, 
                        origin_mat_in, dest_mat_in, n_origins, n_dest,
                        mu.init=NULL, beta.init=NULL, 
                        offset=NULL, lambda=0, gamma=0, delta=0, gmaxsteps=20, 
                        maxsteps=5e3, tol=1e-11, gconv.eps=1e-4, conv.eps=1e-4, 
                        gverbose=FALSE, verbose=FALSE){
  require(genlasso)
  require(Matrix)
  require(parallel)
  require(tidyverse)
  n <- length(y)
  if(!is.null(X)){
    p <- ncol(X)
    coeff <- rep(0, p+1)
    X1 <- cbind(1,X)
  }else{
    p <- 0
    coeff <- 0
    X1 <- matrix(1, nrow=n)
  }
  tX1 <- t(X1)
  num.edge <- nrow(D)
  
  mu <- rep(0, n_origins + n_dest)
  if(!is.null(mu.init)){
    mu.prev <- mu <- mu.init
  } else {
    mu.prev <- mu
  }
  if(!is.null(beta.init)){
    if(length(beta.init)==(p)){
      beta.init <- c(0, beta.init)
    }
    coeff <- beta.init
  } 
  #
  if(!is.null(offset)){
    loff <- log(offset)
  }else{
    loff <- 0
  }  
  erisk.prev <- 1e10
  
  risk <- c(erisk.prev)
  coeff.prev <- coeff
  for(gcount in 1:gmaxsteps){ 
    if(gcount > gmaxsteps) break
    cat(gcount, "th pglasso learning\n")
    ################
    # 1. beta step #
    ################
    coeff <- beta.sol(y, n, p, 
                      (loff + mclapply(fe_mat_indices_list, \(x) {
                        mu[x[1]] + mu[x[2] + n_origins]
                      }, mc.cores = 8) %>% unlist()),
                      X1, tX1, delta, maxsteps, verbose, conv.eps)
    if(gverbose){
      cat("beta=", coeff, "\n")
    }
    hx <- loff + (X1 %*% coeff)
    fx <- hx + map(fe_mat_indices_list, \(x) {
      mu[x[1]] + mu[x[2] + n_origins]
    }) %>% unlist()
    erisk.b <- mean(-y*fx + exp(fx)) + 
      (lambda*sum(abs(D%*%(mu))) + 
         lambda*gamma*sum(abs(mu)) + 
         delta*sum(abs(coeff[-1])))/n
    
    if(erisk.b>erisk.prev+1e-15 | erisk.b>=erisk.prev){
      coeff <- coeff.prev
      break
    }else if(abs(erisk.b-erisk.prev)/(abs(erisk.prev)) <gconv.eps){
      if(gverbose){
        cat("pglasso: beta convergence\n")
        cat(erisk.prev, "=> beta update=> ", erisk.b, "\n")
      }
      break
    }else{
      coeff.prev <- coeff
      erisk.prev <- erisk.b
    }
    
    ###############
    # 2. mu step  #
    ###############
    mout <- mu.sol_mod(n, hx, mu, fe_mat_indices_list, fe_mat_indices, 
                       n_origins, n_dest, y, fe_selection_mat, D, lambda, gamma, 
                       maxsteps, verbose, 10, conv.eps)
    
    if(mout$feasibility==FALSE) {
      return(list(feasibility=FALSE))
    }
    mu <- mout$mu
    print(summary(mu))
    #cat("median", median(mu), "\n")
    fx <- (hx + mclapply(fe_mat_indices_list, \(x) {
      mu[x[1]] + mu[x[2] + n_origins]
    }, mc.cores = 8) %>% unlist())
    eta <- exp(fx)
    erisk <- mean(-y*fx + eta) + (lambda*sum(abs(D%*%mu))) / n
    
    if(erisk > erisk.prev+1e-15 | erisk>=erisk.prev){
      if (gverbose) {
        cat("pglasso didn't help")
      }
      break
    } else
      if(abs(erisk-erisk.prev)/(abs(erisk.prev)) <gconv.eps){
        if(gverbose){
          cat("pglasso: mu convergence\n")
          cat(round(erisk.prev,18), "=> mu update=> ", round(erisk,18), "\n")         
        }
        risk <- c(risk, erisk)
        break
      }
    mu.prev <- mu
    erisk.prev <- erisk
    risk <- c(risk, erisk)
    cat("at step", gcount)
  }
  #cat("median", median(mu), "\n")
  ################
  #cval <- median(mu)
  cval <- 0
  if(cval!=0){
    mu2 <- mu-cval
    if(p==0){
      coeff[1] <- log(sum(y)/sum(exp(loff + mclapply(fe_mat_indices_list, \(x) {
        mu2[x[1]] + mu2[x[2] + n_origins]
      }, mc.cores = 8) %>% unlist())))
      fx <- (loff + mclapply(fe_mat_indices_list, \(x) {
        mu2[x[1]] + mu2[x[2] + n_origins]
      }, mc.cores = 8) %>% unlist() + coeff[1])      
    }else{
      coeff[1] <- log(sum(y) / 
                        sum(exp(loff + (X%*%coeff[-1]) + mclapply(fe_mat_indices_list, \(x) {
                          mu2[x[1]] + mu2[x[2] + n_origins]
                        }, mc.cores = 8) %>% unlist())))
      fx <- loff + mclapply(fe_mat_indices_list, \(x) {
        mu2[x[1]] + mu[x[2] + n_origins]
      }, mc.cores = 8) %>% unlist() + (X1 %*% coeff)
    } 
    eta <- exp(fx)
    if(gverbose){
      cat("Adjusted mu\n")
      erisk <- mean(-y*fx + eta) + (lambda*sum(abs(D%*%mu2)) + 
                                      lambda*gamma*sum(abs(mu2)) + 
                                      delta*sum(abs(coeff[-1])))/n
    }
    #mu <- mu2
    risk <- c(risk, erisk)
  }
  dfm <- length(unique(mu[mu!=0]))
  dfb <- sum(coeff!=0)
  
  clust.infor <- clusters.pglasso(mu) 
  obj <- list(feasibility=TRUE, coefficients=coeff, mu=mu, 
              risk=risk, df=dfb+dfm, dfb=dfb, dfm=dfm, clust.infor=clust.infor)
  obj
}

cal.risk.mu <- function(n, y, fx, eta, lambda, gamma, mu, D){
  mean(-y*fx + eta) + (lambda*(sum(abs(D%*%mu))+gamma*lambda*sum(abs(mu))))/n
}
mu.sol_mod <- function(n, hx, mu, fe_mat_indices_list, fe_mat_indices, 
                       n_origins, n_dest, y, fe_selection_mat, D, lambda, 
                       gamma, maxsteps, verbose, bound = 10, conv.eps=1e-9){
  if(verbose){
    cat("2. mu step\n")
  }  
  df <- 0
  mu_calc <- mclapply(fe_mat_indices_list, \(x) {
    mu[x[1]] + mu[x[2] + n_origins]
  }, mc.cores = 8) %>% unlist()
  fx <- (hx + mu_calc)
  eta <- exp(fx)
  sigma.prev <- sigma <- 2 * (map(1:(n_origins + n_dest), \(x) {
    if (x <= n_origins) {
      eta[which(fe_mat_indices[ ,1] == x)] %>% sum()
    } else {
      eta[which(fe_mat_indices[ ,2] == x - n_origins)] %>% sum()
    }
  }) %>% unlist() %>% max())
  wy <- sqrt(sigma) * (mu + (
    mclapply(1:(n_origins + n_dest), \(x) {
      if (x <= n_origins) {
        ((y - eta)[which(fe_mat_indices[,1] == x)]) %>% 
          sum()
      } else {
        ((y - eta)[which(fe_mat_indices[,2] == x - n_origins)]) %>%
          sum()
      }
    }, mc.cores = 8) %>% unlist()
  ) / sigma)
  
  risk <- c()
  erisk.prev <- cal.risk.mu(n, y, fx, eta, lambda, gamma, mu, D)
  mu.prev <- mu
  
  conv <- FALSE
  count <- 0 
  unusual <- FALSE
  feasibility <- TRUE
  
  if (verbose) {
    print("mu Initialization complete")
  }
  
  while(!conv){
    count <- count+1
    out <- try(fusedlasso(y=wy, X=diag(sqrt(sigma), n_origins+n_dest), 
                          D = D, gamma=0, approx=FALSE, 
                          maxsteps=maxsteps, minlam=lambda, 
                          verbose = verbose), 
               silent=FALSE)
    #print(out$lambda)
    if(class(out)[1]=="try-error"){
      return(list(feasibility=FALSE))
    }   
    if(min(out$lambda)<lambda){
      mu <- (coef(out, lambda=lambda)$beta)
    }else{
      mu <- (coef(out, lambda=min(out$lambda))$beta)    
    }
    print(summary(mu))
    mu[which(abs(mu) > bound)] <- mu[which(abs(mu) > bound)] / 
      abs(mu[which(abs(mu) > bound)]) * bound
    print(summary(mu))
    
    mu_calc <- mclapply(fe_mat_indices_list, \(x) {
      mu[x[1]] + mu[x[2] + n_origins]
    }, mc.cores = 8) %>% unlist()
    fx <- (hx + mu_calc) 
    eta <- exp(fx)
    sigma <- 2 * (map(1:(n_origins + n_dest), \(x) {
      if (x <= n_origins) {
        eta[which(fe_mat_indices[ ,1] == x)] %>% sum()
      } else {
        eta[which(fe_mat_indices[ ,2] == x - n_origins)] %>% sum()
      }
    }) %>% unlist() %>% max())
    
    if(abs(sigma-sigma.prev)/sigma.prev<conv.eps){
      cat(count, "th sigma=", sigma, sigma.prev, "\n")
      break  
    }
    erisk <- cal.risk.mu(n, y, fx, eta, lambda, gamma, mu, D)
    cat(count, "th risk=", erisk, "\n")
    if(verbose){
      cat("mu.max.sol>>>", count, "th objective", 
          erisk.prev, erisk, abs(erisk-erisk.prev)/(abs(erisk.prev)), "\n")
    }
    if(erisk > erisk.prev+1e-10){
      conv <- TRUE
      unusual <- TRUE
      print("unusual!")
    }else if(abs(erisk-erisk.prev)/(abs(erisk.prev)) < conv.eps){
      conv <- TRUE
      if(verbose){
        cat("mu step convergence\n")
      }
    }else if(count > maxsteps){
      conv <- TRUE
    }     
    if(!unusual){
      sigma.prev <- sigma
      mu.prev <- mu
      erisk.prev <- erisk
      risk <- c(risk, erisk)
    }else{
      #mu <- mu.prev
      break
    }  
  } # while
  df <- length(unique(mu[mu!=0]))
  list(mu=mu, df=df, feasibility=TRUE)
}

beta.sol <- function(y, n, p, hx, X1, tX1, delta, maxsteps, 
                     verbose, conv.eps=1e-9){
  erisk.prev <- 1e6
  risk <- c()
  
  if(verbose){
    cat("1. beta step\n")
  }
  if(p==0){
    coefficients <- log(sum(y)/sum(exp(hx)))
    return(coefficients)
  }  
  coeff.prev <- coeff <- rep(0, p+1)
  
  for(count in 1:maxsteps){
    fx <- hx + (X1 %*% coeff)
    eta <- exp(fx)        
    wCov <- tX1 %*% (as.vector(eta) * X1)
    sigma <- max(eigen(wCov)$values)
    wy <- (coeff + tX1%*%(y-eta)/sigma)
    
    coeff[1] <- wy[1]
    coeff[-1] <- sign(wy[-1])*pmax(abs(wy[-1])-delta/sigma,0)
    
    fx <- hx + (X1 %*% coeff)
    erisk <- mean(-y*fx + exp(fx)) + delta/n*sum(abs(coeff[-1]))
    
    if(verbose){
      cat(count, "th objective", erisk, abs(erisk-erisk.prev)/(abs(erisk.prev)), "\n")
    }
    if(erisk > erisk.prev+1e-10 | erisk>=erisk.prev){
      break
    }else if(abs(erisk-erisk.prev)/(abs(erisk.prev)) < conv.eps){
      if(verbose){
        cat("beta step convergence\n")
      }
      break
    }     
    coeff.prev <- coeff  
    erisk.prev <- erisk
    risk <- c(risk, erisk)
  } # while
  coefficients <- coeff
  return(coefficients)
}
clusters.pglasso <- function(mu){
  # mu=0 => 0 group
  n <- length(mu)
  if(sum(mu==0)==n){
    out <- list()
    out$membership <- rep(0,n)
    out$csize <- n
    out$upar <- 0
    return(out)
  }
  umu <- sort(unique(mu),decreasing=FALSE)
  ng <- length(umu)
  
  out <- list()
  out$membership <- 1:n
  out$csize <- rep(0,ng)
  out$upar <- rep(0,ng)
  
  umu0 <- setdiff(umu,0)
  
  if(length(umu0) == ng){
    for(i in 1:ng){
      flag <- mu==umu[i]
      out$membership[flag] <- i
      out$csize[i] <- sum(flag)
      out$upar[i] <- umu[i]
    }
    out$no <- ng
  }else{
    for(i in 1:(ng-1)){
      flag <- mu==umu0[i]
      out$membership[flag] <- i      
      out$csize[i+1] <- sum(flag)
      out$upar[i+1] <- umu0[i]      
    }
    flag <- mu==0
    out$membership[flag] <- 0
    out$csize[1] <- sum(flag)
    out$no <- ng-1
  }
  return(out)
}

mse_calc <- function(mu, epsilon, X_ij, delta, fe_mat_indices_list, 
                     origin_count) {
  mu_calc <- map(fe_mat_indices_list, \(x) {
    mu[x[1]] + mu[x[2] + length(origin_list)]
  }) %>% unlist()
  (X_ij - (exp(mu_calc) * delta^(-epsilon)))^2 %>% mean()
}