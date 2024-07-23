library(astsa)
library(glmnet)
library(doParallel)
library(MCMCpack)
library(readr)
library(bayesreg)
library(forecast)
library(rqPen)
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

chen_chan <- function(x, h = 1, long.ar.select, maxP, maxQ, updateMA = FALSE, ic = c("BIC", "AIC", "AICc"),eta = 2, alpha = seq(0, 1, 0.1), Method = c("LADADLASSO", "ADLASSO","ADENET","LADADENET","RIDGE"), w.Method = c("LASSO", "RIDGE", "LS")) {
  require(glmnet)
  require(quantreg)
  Method <- match.arg(Method)
  w.Method <- match.arg(w.Method)
  ic <- match.arg(ic)
  
  Nt <- length(x)
  
  
  init.mod.est <- ar(x, aic = TRUE, order.max = Nt - 1, demean = TRUE)
  init.mod.error <- residuals(init.mod.est)
  init.mod.order <- length(which(is.na(init.mod.error))) 
  
  m <- init.mod.est$order + max(maxP, maxQ) + 1
  
  dataP <- foreach(p = 1:maxP, .combine = cbind) %do% {
    lag.func(x, m, p)
  }
  dataQ <- foreach(q = 1:maxQ, .combine = cbind) %do% {
    lag.func(init.mod.error, m, q)
  }
  
  first.modX <- as.matrix(cbind(dataP, dataQ))
  first.y <- x[m:(Nt - 1)]
  
  if (all(diff(first.y)==0)){
    first.y <- first.y + rnorm(length(first.y),mean=0,sd=0.1) #important when we have constant vector, cant stadarize by the glmnet
  }
  
  n.alpha <- length(alpha)
  
  if (Method == "RIDGE") {
    first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 0))
    first.mod.RSS <- colSums((first.y - predict(first.mod.est, newx = first.modX))^2)
    
    if (ic=='BIC') {
      first.bic.out <- log(length(first.y)) * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.bic.out)]
    } else if (ic=='AIC') {
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
    }
    
    final.mod.coef <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[-1]
    final.mod.int <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[1]
    final.mod.s2 <- sum((first.y - predict(first.mod.est, newx = first.modX, s = first.mod.lambda, method = "lambda"))^2) / (length(first.y) - sum(final.mod.coef != 0) - 1)
    
    out <- list(final.mod.coef = final.mod.coef, final.mod.int = final.mod.int, final.mod.s2 = final.mod.s2)
    return(out)
  }
  else if (Method == "ADLASSO") {
    if (w.Method == "RIDGE") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 0))
    } else if (w.Method == "LASSO") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 1))
    } else if (w.Method == "LS") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE,alpha = 0, lambda = 0))
    }
    first.mod.RSS <- colSums((first.y - predict(first.mod.est, newx = first.modX))^2)
    
    if (ic=='BIC') {
      first.bic.out <- log(length(first.y)) * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.bic.out)]
    } else if (ic=='AIC') {
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
    }
    
    first.mod.coef <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[-1]
    first.mod.mu <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[1]
    weights <- abs(first.mod.coef + 1 / length(first.y))^(-eta)
    weights[is.infinite(weights)]<-10e10
    
    
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
    } else {
      second.modX <- first.modX
      second.y <- first.y
    }
    
    second.mod.est <- glmnet(y = second.y, x = second.modX, standardize = TRUE, alpha = 1, thresh = 1e-16, penalty.factor = weights)
    second.mod.RSS <- colSums((second.y - predict(second.mod.est, newx = second.modX))^2)
    
    if(ic=='BIC'){
      second.bic.out <- log(length(second.y)) * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
      second.mod.lambda <- second.mod.est$lambda[which.min(second.bic.out)]
    } else if (ic=='AIC') {
      second.aic.out <- 2 * second.mod.est$df + length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
      second.mod.lambda <- second.mod.est$lambda[which.min(second.aic.out)]
    }
    
    final.mod.coef <- as.numeric(coef(second.mod.est, s = second.mod.lambda, method = "lambda"))[-1]
    nonzero.select <- which(final.mod.coef != 0)
    final.mod.int <- as.numeric(coef(second.mod.est, s = second.mod.lambda, method = "lambda"))[1]
    final.mod.s2 <- sum((second.y - predict(second.mod.est, newx = second.modX, s = second.mod.lambda, method = "lambda"))^2) / (length(second.y) - sum(final.mod.coef[nonzero.select] != 0) - 1)
    
    out <- list(final.mod.coef = final.mod.coef, final.mod.int = final.mod.int, nonzero.select = nonzero.select, final.mod.s2 = final.mod.s2)
    return(out)
  }else if (Method == "ADENET") {
    # Initial LASSO Weights
    if (w.Method == "RIDGE") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 0))
    } else if (w.Method == "LASSO") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 1))
    } else if (w.Method == "LS") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE,alpha = 0, lambda = 0))
    }
    first.mod.RSS <- colSums((first.y - predict(first.mod.est, newx = first.modX))^2)
    
    if (ic=='BIC') {
      first.bic.out <- log(length(first.y)) * first.mod.est$df +length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.bic.out)]
      first.cv.out <- c(1, first.mod.lambda, min(first.bic.out))
    } else if (ic=='AIC'){
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.RSS) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
      first.cv.out <- c(1, first.mod.lambda, min(first.aic.out))
    }
    
    first.mod.alpha <- 1
    first.mod.lambda <- first.cv.out[2]
    first.mod.est <- glmnet(y = first.y, x = first.modX, standardize = TRUE,alpha = first.mod.alpha, lambda = first.mod.lambda)
    first.mod.coef <- as.numeric(coef(first.mod.est))[-1]
    first.mod.mu <- as.numeric(coef(first.mod.est))[1]
    weights <- abs(first.mod.coef + 1/length(first.y)) ^ (-eta)
    weights[is.infinite(weights)]<-10e10
    
    # Update Model Matrix of AR and MA terms Based off Initial Estimation
    if (updateMA) {
      update.mod.predict <- rep(NA, length(x))
      update.mod.error <- rep(0, length(x))
      
      for (v in (h + max(maxP, maxQ)):Nt) {
        update.mod.predict[v] <- first.mod.mu +
          x[(v - h):(v - maxP - h + 1)] %*% first.mod.coef[1:maxP] +
          update.mod.error[(v - h):(v - maxQ - h + 1)] %*% first.mod.coef[-(1:maxP)]
        update.mod.error[v] <- x[v] - update.mod.predict[v]
      }
      
      update.dataQ <- foreach(q = 1:maxQ, .combine = cbind) %do% {
        lag.func(update.mod.error, k = (q + h - 1))
      }
      
      second.modX <- as.matrix(cbind(dataP, update.dataQ))[-(1:(max(maxP, maxQ) + h - 1)), ]
      second.y <- x[-(1:(max(maxP, maxQ) + h - 1))]
    } else {
      second.modX <- first.modX
      second.y <- first.y
    }
    
    # search Through All Lambdas
    second.cv.out <- foreach(a = 1:n.alpha, .combine = rbind) %do% {
      second.mod.est <- glmnet(y = second.y, x = second.modX, standardize = TRUE,alpha = alpha[a], penalty.factor = weights)
      second.mod.RSS <- colSums((second.y - predict(second.mod.est, newx = second.modX))^2)
      
      if (ic=='BIC') {
        second.bic.out <- log(length(second.y)) * second.mod.est$df +length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
        second.mod.lambda <- second.mod.est$lambda[which.min(second.bic.out)]
        result <- c(alpha[a], second.mod.lambda, min(second.bic.out))
      } else if (ic=='AIC'){
        second.aic.out <- 2 * second.mod.est$df +length(second.y) * log(as.vector(second.mod.RSS) / length(second.y))
        second.mod.lambda <- second.mod.est$lambda[which.min(second.aic.out)]
        result <- c(alpha[a], second.mod.lambda, min(second.aic.out))
      }
      result
    }
    
    second.mod.alpha <- alpha[which.min(second.cv.out[, 3])]
    second.mod.lambda <- second.cv.out[which.min(second.cv.out[, 3]), 2]
    second.mod.est <- glmnet(y = second.y, x = second.modX,standardize = TRUE, alpha = second.mod.alpha,lambda = second.mod.lambda, penalty.factor = weights)
    second.mod.coef <- as.numeric(coef(second.mod.est))[-1]
    second.mod.mu <- as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef <- second.mod.coef
    nonzero.select <- which(final.mod.coef != 0)
    final.mod.int <- second.mod.mu
    final.mod.s2 <- sum((second.y - predict(second.mod.est, newx = second.modX))^2) /(length(second.y) - sum(final.mod.coef[nonzero.select] != 0) - 1)
    
    out <- list(final.mod.coef = final.mod.coef,final.mod.int = final.mod.int, nonzero.select = nonzero.select,final.mod.s2 = final.mod.s2)
    return(out)
  }else if (Method == "LADADLASSO") {
    
    lambda_seq <- 10^seq(3, -2, by = -0.01)
    if (w.Method == "RIDGE") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 0))
    } else if (w.Method == "LASSO") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 1))
    } else if (w.Method == "LS") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE,alpha = 0, lambda = 0))
    }
    first.mod.SAD <- colSums(abs(first.y - predict(first.mod.est, newx = first.modX)))
    
    if (ic=='BIC') {
      first.bic.out <- log(length(first.y)) * first.mod.est$df + length(first.y) * log(as.vector(first.mod.SAD) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.bic.out)]
    } else if (ic=='AIC') {
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.SAD) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
    }
    
    first.mod.coef <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[-1]
    first.mod.mu <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[1]
    weights <- abs(first.mod.coef + 1 / length(first.y))^(-eta) ###acquired from LASSO
    weights[is.infinite(weights)]<-10e10
    
    rqpen_model <- rq.pen(as.matrix(first.modX), first.y, lambda = lambda_seq, penalty = "aLASSO",tau=0.5,a=1,penalty.factor = weights)###Ridge weights are the default in rqpen
    select_on_aic<-qic.select(rqpen_model,method=ic)###change to BIC when selecting on BIC
    aic_values<- select_on_aic$ic ### store the aic or bic values
    final.mod.int <- as.numeric(select_on_aic$coefficients)[1]
    final.mod.coef <- as.numeric(select_on_aic$coefficients)[-1] 
    nonzero.select <- which(final.mod.coef != 0)
    df <- sum(final.mod.coef != 0)  
    
    #store the results
    out <- list(final.mod.coef = final.mod.coef, final.mod.int = final.mod.int, nonzero.select = nonzero.select)
    return(out)
  }else if (Method == "LADADENET") {
    lambda_seq <- 10^seq(3, -2, by = -0.01)
    if (w.Method == "RIDGE") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 0))
    } else if (w.Method == "LASSO") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE, alpha = 1))
    } else if (w.Method == "LS") {
      first.mod.est <- suppressWarnings(glmnet(y = first.y, x = first.modX, standardize = TRUE,alpha = 0, lambda = 0))
    }    
    first.mod.SAD <- colSums(abs(first.y - predict(first.mod.est, newx = first.modX)))
    
    if (ic=='BIC') {
      first.bic.out <- log(length(first.y)) * first.mod.est$df + length(first.y) * log(as.vector(first.mod.SAD) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.bic.out)]
    } else if (ic=='AIC') {
      first.aic.out <- 2 * first.mod.est$df + length(first.y) * log(as.vector(first.mod.SAD) / length(first.y))
      first.mod.lambda <- first.mod.est$lambda[which.min(first.aic.out)]
    }
    
    first.mod.coef <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[-1]
    first.mod.mu <- as.numeric(coef(first.mod.est, s = first.mod.lambda, method = "lambda"))[1]
    weights <- abs(first.mod.coef + 1 / length(first.y))^(-eta) ###acquired from LASSO
    weights[is.infinite(weights)]<-10e10
    
    rqpen_model <- rq.pen(as.matrix(first.modX), first.y, lambda = lambda_seq, penalty = "ENet",a=0.5,tau=0.5,penalty.factor = weights)
    select_on_aic<-qic.select(rqpen_model,method=ic)###change to BIC when selecting on BIC
    aic_values<- select_on_aic$ic ### store the aic or bic values
    final.mod.int <- as.numeric(select_on_aic$coefficients)[1]
    final.mod.coef <- as.numeric(select_on_aic$coefficients)[-1] 
    nonzero.select <- which(final.mod.coef != 0)
    df <- sum(final.mod.coef != 0)  
    
    #store the results
    out <- list(final.mod.coef = final.mod.coef, final.mod.int = final.mod.int, nonzero.select = nonzero.select)
    return(out)
  }
  
}



############Simulations Models - Contaminated Data###################################################


model.one <- function(n,p,d){
  bernoulli.samples <- rbinom(n,1,p)
  w <- rnorm(n)
  for (u in 1:n){
    if (bernoulli.samples[u]==1){
      w[u] = rnorm(1,mean = 0,sd=1*d )
    }
  }
  
  
  y <- numeric(n)
  
  y[1] = w[1]
  for (t in 2:n) {
    if (t < 8) {
      y[t] = w[t] + 0.8 * y[t-1]
    } else {
      y[t] = w[t] + 0.8 * y[t-1] + 0.7 * y[t-6] - 0.56 * y[t-7]
    }
  }
  
  true.coef <- c(0.8,0,0,0,0,0.7,-0.56,0,0,0,0,0,0,0)
  number.coef <- 3
  out <- list(y, true.coef, number.coef)
  #return(y)
  return(out)
}


model.two <- function(n,p,d){
  bernoulli.samples <- rbinom(n,1,p)
  w <- rnorm(n)
  for (u in 1:n){
    if (bernoulli.samples[u]==1){
      w[u] = rnorm(1,mean = 0,sd=1*d )
    }
  }
  
  y <- numeric(n)
  y[1] = w[1]  
  for (t in 2:n) {
    if (t < 8) {
      y[t] = w[t] + 0.8 * y[t-1] + 0.8 * w[t-1]
    } else {
      y[t] = w[t] + 0.8 * y[t-1] + 0.7 * y[t-6] - 0.56 * y[t-7] + 0.8 * w[t-1] + 0.7 * w[t-6] + 0.56 * w[t-7]
    }
  }

  true.coef <- c(0.8, 0,0,0,0,0.7,-0.56,0.8,0,0,0,0,0.7,0.56)
  number.coef <- 6
  out<-list(y, true.coef, number.coef)
  return(out)
}

model.three <- function(n,p,d){
  bernoulli.samples <- rbinom(n,1,p)
  w <- rnorm(n)
  for (u in 1:n){
    if (bernoulli.samples[u]==1){
      w[u] = rnorm(1,mean = 0,sd=1*d )
    }
  }
  
  y <- numeric(n)
  for (t in 1:n) {
    if (t < 8) {
      y[t] = w[t]
    } else {
      y[t] = w[t] + 0.8 * w[t-1] + 0.7 * w[t-6] + 0.56 * w[t-7]
    }
  }

  true.coef <- c(0, 0,0,0,0,0,0,0.8,0,0,0,0,0.7,0.56)
  number.coef <- 3
  out<-list(y, true.coef, number.coef)
  return(out)
}

model.four <- function(n, p, d) {
  bernoulli.samples <- rbinom(n, 1, p)
  w <- rnorm(n)
  for (u in 1:n){
    if (bernoulli.samples[u]==1){
      w[u] = rnorm(1,mean = 0,sd=1*d )
    }
  }
  
  y <- numeric(n)
  
  for (t in 1:n) {
    if (t < 13) {
      y[t] = w[t]
    } else {
      y[t] = w[t] - 0.6 * w[t-1] - 0.8 * w[t-12]
    }
  }
  

  
  true.coef <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8)
  number.coef <- 2
  out <- list(y, true.coef, number.coef)
  return(out)
}


# Set parameters
n <- 360  # Number of samples
p <- 0.1  # Probability for Bernoulli samples
d <- 1    # Standard deviation for noise


################## Evaluations of the Chen & Chan method plus the Adaptive LAD and Elastic net method ##########################

n <- c(120, 240, 360) # simulation sample size
MC <- 1000 # Monte Carlo number of experiments


ic <- c("BIC","AIC") # Information critiria 
w.method <- c("LASSO","RIDGE","LS") # Method of calculating the Weights 
models <- c("one","two","three","four") #simulatuons Methods

max.p <- 14 # maximum AR order
max.q <- 14 # maximum MA order

for (model.c in models) {
  cat('Model:', model.c, "\n")
  for (icc in ic) {
    cat('Information Criterion:', icc, "\n")
    for (w in w.method) {

      cat('Weighting Method:', w, "\n")

      for (i in n) {
        corect.model <- 0
        significant.variables <- 0
        true.positive.rate <- numeric()
        false.positive.rate <- numeric()
        false.negative.rate <- numeric()
        for (mc in 1:MC) {
          
          model_function <- get(paste0("model.", model.c))
          temp <- model_function(i, p = 0.1, d = 1)  
          y <- temp[[1]]
          true.coef <- temp[[2]]
          number.coef <- temp[[3]]
          model <- chen_chan(x = y, long.ar.select = TRUE, maxP = max.p, maxQ = max.q, ic = icc, Method = "LADADENET", w.Method = w)
          
          false.positive.rate <- c(false.positive.rate, sum(model$final.mod.coef != 0 & true.coef == 0) / (sum(model$final.mod.coef != 0 & true.coef == 0) + sum(model$final.mod.coef != 0 & true.coef != 0)))
          false.negative.rate <- c(false.negative.rate, sum(model$final.mod.coef == 0 & true.coef != 0) / (sum(model$final.mod.coef != 0 & true.coef != 0) + sum(model$final.mod.coef == 0 & true.coef == 0)))
          true.positive.rate <- c(true.positive.rate, sum(model$final.mod.coef != 0 & true.coef != 0) / (sum(model$final.mod.coef != 0 & true.coef != 0) + sum(model$final.mod.coef != 0 & true.coef == 0)))
          

          
          # Initialize the order variables
          ar.order <- FALSE
          ma.order <- FALSE
          
          # Check AR order 
          if (any(true.coef[1:max.p] != 0) && any(model$final.mod.coef[1:max.p] != 0)) {
            ar.order <- tail(which(true.coef[1:max.p] != 0), 1) == tail(which(model$final.mod.coef[1:max.p] != 0), 1)
          } else if (!any(true.coef[1:max.p] != 0) && !any(model$final.mod.coef[1:max.p] != 0)) {
            ar.order <- TRUE
          }
          
          # Check MA order 
          if (any(true.coef[(max.p+1):(max.p*2)] != 0) && any(model$final.mod.coef[(max.p+1):(max.p*2)] != 0)) {
            ma.order <- tail(which(true.coef[8:14] != 0), 1) == tail(which(model$final.mod.coef[(max.p+1):(max.p*2)] != 0), 1)
          } else if (!any(true.coef[(max.p+1):(max.p*2)] != 0) && !any(model$final.mod.coef[(max.p+1):(max.p*2)] != 0)) {
            ma.order <- TRUE
          }
          
          if (ar.order && ma.order) {
            corect.model <- corect.model + 1
          }
        }
        
        corect.model.rate <- corect.model / MC
        cat("n =", i, "A:", median(true.positive.rate, na.rm = TRUE), "T:", corect.model.rate, '+', median(false.positive.rate, na.rm = TRUE), '-', median(false.negative.rate, na.rm = TRUE), '\n')
      }
    }
  }
}













