#moving window
library(glmnet)
library(doParallel)
library(MCMCpack)
library(readr)
library(bayesreg)
library(ggplot2)
library(gridExtra)
library(forecast)
library(hqreg)
library(stats)
library(forecastML)

lag.func <- function(x, k,pq) {
  t = length(x)
  n = t-k #collumn length
  y=c(rep(NA,n))
  for (i in k:t){
    y[i-k] = x[i-pq-1]
  }
  return(y)
}


adshrink123 <- function(x, h=1, long.ar.select, maxP, maxQ, updateMA=FALSE, ic = c("BIC", "AIC", "AICc"), eta=2, Method= c("ADLASSO", "ADRIDGE", "ADENET"), w.Method=c("LASSO", "LS", "RIDGE")) {
  # x : the time series
  # Information criteria, choose one of the three options
  # Methods to srink the coefficients "ADLASSO" adaptive lasso, "ADRIDGE" for adaptive Ridge & "ADENET" adaptive elastic net
  # Weight methods, "LASSO", "LS" & RIDGE
  
  require(glmnet)
  
  Method <- match.arg(Method)
  w.Method <-match.arg(w.Method)
  ic <- match.arg(ic)
  
  
  Nt <- length(x)

  #max.ar.order <- ceiling(10 * log10(Nt)) # returns the smallest integer 
  init.mod.est <- ar(x, aic=TRUE, order.max=Nt-1, demean=TRUE) 
  init.mod.error <- residuals(init.mod.est)
  init.mod.order <- length(which(is.na(init.mod.error))) 
  
  m <- init.mod.est$order + max(maxP,maxQ) + 1
  
  dataP <- foreach(p = 1:maxP, .combine = cbind) %do% {
    lag.func(x, m,p)  #ypirxe ena +H
  }
  dataQ <- foreach(q = 1:maxQ, .combine = cbind) %do% {
    lag.func(init.mod.error, m,q) #ypirxe ena +H
  }
  
  #first.modX <- as.matrix(cbind(dataP, dataQ))[-(1:(init.mod.order + max(maxP, maxQ) + h - 1)), ]
  first.modX <- as.matrix(cbind(dataP, dataQ))
  #first.y <- x[-(1:(init.mod.order + max(maxP, maxQ) + h - 1))]
  first.y <- x[m:(Nt-1)]
  if (all(diff(first.y)==0)){
    first.y <- first.y + rnorm(length(first.y),mean=0,sd=0.1) #important when we have constant vectore, cant stadarize by the glmnet
  }
 

  
  if (w.Method == "LASSO") {
    
    first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 1))
    first.mod.RSS <- colSums((first.y - predict(first.mod.est, newx = first.modX))^2)
    
    if (ic == "BIC") {
      first.bic.out <- log(length(first.y)) * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.bic.out)]
      
    } else if (ic == "AICc"){

      k = length(coef(first.mod.est)!=0)
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y)) + (2*k*(k+1))/(Nt - k - 1)
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
    }
    else {
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
      
    }
    
    if (Method == "ADENET"){
      first.cv.out = c(1, first.mod.lambda, min(first.aic.out))
      alpha=seq(0, 0.9, 0.1)
      n.alpha <- length(alpha)
    }
    
    first.mod.coef <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[-1]
    first.mod.mu <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[1]
    weights <- abs(first.mod.coef + 1 / length(first.y)) ^ (-eta) 
  
    
  }else if( w.Method == "RIDGE"){

    
    first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 0))
    first.mod.RSS <- colSums((first.y - predict(first.mod.est, newx = first.modX))^2)
    
    if (ic == "BIC") {
      first.bic.out <- log(length(first.y)) * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.bic.out)]
    } else if (ic == "AICc"){
      k = length(coef(first.mod.est)!=0)
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y)) + (2*k*(k+1))/(Nt - k - 1)
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
    }else{
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
      
    }
    if (Method == "ADENET"){
      first.cv.out = c(0, first.mod.lambda, min(first.aic.out)) #alpha is zero here
      alpha=seq(0.1, 1, 0.1)
      n.alpha <- length(alpha)
    }
    
    first.mod.coef <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[-1]
    first.mod.mu <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[1]
    weights <- abs(first.mod.coef + 1 / length(first.y)) ^ (-eta) 

    
  }else if (w.Method == "LS"){
    
    first.mod.est <- lm(first.y ~ first.modX)
    weights <- abs(coef(first.mod.est)[-1] +1/length(first.y))
    weights[is.na(weights)] <- 1/length(first.y)^2
    weights <- weights^(-eta)
    
    if (Method == "ADENET"){
      alpha=seq(0.0, 1, 0.1)
      n.alpha <-length(alpha)
    }
    
    
  }
   
    

  if (updateMA) {
      update.mod.predict <- rep(NA, length(x))
      update.mod.error <- rep(0, length(x))
      
      for (v in (h + max(maxP, maxQ)):Nt) {
        update.mod.predict[v] <- first.mod.mu + x[(v - h):(v - maxP - h + 1)] %*% first.mod.coef[1:maxP] + 
        update.mod.error[(v - h):(v - maxQ - h + 1)] %*% first.mod.coef[-(1:maxP)]
        update.mod.error[v] <- x[v] - update.mod.predict[v]
      }
      
      update.dataQ <- foreach(q = 1:maxQ, .combine = cbind) %do% { lag.func(update.mod.error, k = (q + h - 1)) }
      second.modX <- as.matrix(cbind(dataP, update.dataQ))[-(1:(max(maxP, maxQ) + h - 1)), ]
      second.y <- x[-(1:(max(maxP, maxQ) + h - 1))]
    }else {
      second.modX <- first.modX
      second.y <- first.y
      
  }
  if (Method == "ADLASSO"){
    
    
      second.mod.est <- glmnet(y = second.y, x = second.modX, standardize = TRUE, alpha = 1, thresh = 1e-16, penalty.factor = abs(weights))
      second.mod.RSS <- colSums((second.y - predict(second.mod.est, newx = second.modX))^2)
      
      if (ic == "BIC") {
        second.bic.out <- log(length(second.y)) * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
        second.mod.lambda <- second.mod.est$lambda[which.min(second.bic.out)]
      } else if (ic == "AICc"){
        k = length(coef(second.mod.est)!=0)
        second.aic.out <- 2 * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y)) + (2*k*(k+1))/(Nt - k - 1)
        second.mod.lambda <- first.mod.est$lambda[which.min(second.aic.out)]
      }else {
        second.aic.out <- 2 * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
        second.mod.lambda <- second.mod.est$lambda[which.min(second.aic.out)]
      }
      
      final.mod.coef <- as.numeric(coef(second.mod.est, s = second.mod.lambda, method = "lambda"))[-1]
      nonzero.select <- which(final.mod.coef != 0)
      final.mod.int <- as.numeric(coef(second.mod.est, s = second.mod.lambda, method = "lambda"))[1]
      final.mod.s2 <- sum((second.y - predict(second.mod.est, newx = second.modX,
                                              s = second.mod.lambda, method = "lambda"))^2) / (length(second.y) - sum(final.mod.coef[nonzero.select] != 0) - 1)
      
      
    }else if(Method == "ADRIDGE"){
      
      
      second.mod.est <- glmnet(y = second.y, x = second.modX, standardize = TRUE, alpha = 0, thresh = 1e-16, penalty.factor = abs(weights))
      second.mod.RSS <- colSums((second.y - predict(second.mod.est, newx = second.modX))^2)
      
      if (ic == "BIC") {
        second.bic.out <- log(length(second.y)) * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
        second.mod.lambda <- second.mod.est$lambda[which.min(second.bic.out)]
      } else if (ic == "AICc"){
        k = sum(coef(second.mod.est)!=0)
        second.aic.out <- 2 * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y)) + (2*k*(k+1))/(Nt - k - 1)
        second.mod.lambda <- first.mod.est$lambda[which.min(second.aic.out)]
      }else {
        second.aic.out <- 2 * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
        second.mod.lambda <- second.mod.est$lambda[which.min(second.aic.out)]
      }
      
      final.mod.coef <- as.numeric(coef(second.mod.est, s = second.mod.lambda, method = "lambda"))[-1]
      nonzero.select <- which(final.mod.coef != 0)
      final.mod.int <- as.numeric(coef(second.mod.est, s = second.mod.lambda, method = "lambda"))[1]
      final.mod.s2 <- sum((second.y - predict(second.mod.est, newx = second.modX,
                                              s = second.mod.lambda, method = "lambda"))^2) / (length(second.y) - sum(final.mod.coef[nonzero.select] != 0) - 1)
      
      
    }else if (Method == "ADENET"){
      
    
      
        # Final Elastic Net Estimates (Search Through All Lambdas)
      second.cv.out = foreach(a = 1:n.alpha, .combine = rbind) %do% {
         
        second.mod.est = glmnet(y = second.y, x = second.modX, standardize = TRUE,
                                    alpha = alpha[a], penalty.factor = weights)
        second.mod.RSS = colSums((second.y - predict(second.mod.est, newx = second.modX))^2)
            
        if (ic == "BIC") {
          second.bic.out = log(length(second.y)) * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
          second.mod.lambda = second.mod.est$lambda[which.min(second.bic.out)]
          result = c(alpha[a], second.mod.lambda, min(second.bic.out))
        } else if (ic == "AICc"){
          k = length(coef(second.mod.est)!=0)
          second.aic.out <- 2 * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y)) + (2*k*(k+1))/(Nt - k - 1)
          second.mod.lambda <- first.mod.est$lambda[which.min(second.aic.out)]
        }else {
          second.aic.out = 2 * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
          second.mod.lambda = second.mod.est$lambda[which.min(second.aic.out)]
          result = c(alpha[a], second.mod.lambda, min(second.aic.out))
          }
        result 
      }
      
      if (w.Method == "LS"){
        second.mod.alpha = alpha[which.min(second.cv.out[,3])]
        second.mod.lambda = second.cv.out[which.min(second.cv.out[,3]),2]
      }else {
        if (first.cv.out[3] > min(second.cv.out[,3])){
          second.mod.alpha = alpha[which.min(second.cv.out[,3])]
          second.mod.lambda = second.cv.out[which.min(second.cv.out[,3]),2]
        }else{
            second.mod.alpha = first.cv.out[1]
            second.mod.lambda = first.cv.out[2]
        }
      }
      second.mod.est = glmnet(y = second.y, x = second.modX, standardize = TRUE, alpha = second.mod.alpha, lambda = second.mod.lambda, penalty.factor = weights)
      second.mod.coef <- as.numeric(coef(second.mod.est))[-1]
      second.mod.mu <- as.numeric(coef(second.mod.est))[1]
      final.mod.coef <- second.mod.coef
      nonzero.select <- which(final.mod.coef != 0)
      final.mod.int <- second.mod.mu
      final.mod.s2 <- sum((second.y - predict(second.mod.est, newx = second.modX))^2) / (length(second.y) - sum(final.mod.coef[nonzero.select] != 0) - 1)
    
      
    }
    
    
    out <- list(final.mod.coef = final.mod.coef, final.mod.int = final.mod.int, final.mod.s2 = final.mod.s2, nonzero.select = nonzero.select)
    
    return(out)
  } 
