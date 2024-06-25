
n <- c(120, 240, 360, 1600)


MC <- 1000
model.one <- function(n){
  w <- rnorm(n)
  
  y = rep(NA, n)
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
  return(out)
}

model.two <- function(n){
  w <- rnorm(n)
  y = rep(NA, n)
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

model.three <- function(n){
  w <- rnorm(n)
  y = rep(NA, n)
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

model.four <- function(n){
  w <- rnorm(n)
  y <- w
  for (t in 13:n){
    y[t] = w[t] -0.6*w[t-1] - 0.8*w[t-12]
  }
  true.coef <- c(0,0,0,0,0,0,0,0.6,0,0,0,0.8,0,0)
  number.coef <- 2
  out<-list(y, true.coef, number.coef)
  return(out)
}




ic <-c("AIC","BIC")

models <- c("one", "two", "three", "four") 
models<-c("four")
#w.method <- c("LASSO")
for (model.c in models){
  cat('model', model.c,"\n")
  for (icc in ic){
    cat(icc, "\n")
    if (icc=="BIC"){
      w.method <- c("LASSO", "RIDGE", "LS")
    }else {
      w.method <- c("LASSO")
    }
    for (w in w.method){
      cat(w,"\n")
      for (i in n){
        corect.model <- 0 
        true.positive.rate <- numeric()
        false.positive.rate <- numeric()
        false.negative.rate <- numeric()
        for (mc in 1:MC){
  
          model_function <- get(paste0("model.", model.c))
          temp <- model_function(i)
          y <- as.double(unlist(temp[1]))
          true.coef <-  as.double(unlist(temp[2]))
          number.coef <-  as.double(unlist(temp[3]))
          window.model <-adshrink123(x = y , long.ar.select = TRUE , maxP = 7, maxQ = 7,ic = icc, Method = "ADENET", w.Method=w)
          false.positive.rate <- c(false.positive.rate, sum(window.model$final.mod.coef !=0 & true.coef==0)/( sum(window.model$final.mod.coef !=0 & true.coef==0) + sum(window.model$final.mod.coef !=0 & true.coef!=0) ))
          false.negative.rate <- c(false.negative.rate, sum(window.model$final.mod.coef == 0 & true.coef !=0)/(sum(window.model$final.mod.coef !=0 & true.coef !=0) + sum(window.model$final.mod.coef ==0 & true.coef==0) ))
          #A
          true.positive.rate <- c(true.positive.rate, sum(window.model$final.mod.coef != 0 & true.coef != 0 )/(sum(window.model$final.mod.coef != 0 & true.coef != 0 ) + sum(window.model$final.mod.coef != 0 & true.coef == 0 )) )
          
          if (all(true.coef[1:7] ==0)==FALSE && all(window.model$final.mod.coef[1:7] ==0)==FALSE){
            ar.order <- tail(which(true.coef[1:7] != 0),1)==tail(which(window.model$final.mod.coef[1:7] !=0),1)
            
          } else if (all(true.coef[1:7] ==0)==TRUE && all(window.model$final.mod.coef[1:7] ==0)==TRUE) {ar.order=FALSE}
          else if (all(true.coef[1:7] ==0)==TRUE && all(window.model$final.mod.coef[1:7] ==0)==FALSE){ar.order=FALSE}
          
          if (all(true.coef[8:14] ==0)==FALSE && all(window.model$final.mod.coef[8:14] ==0)==FALSE){
            ma.order <- tail(which(true.coef[8:14] != 0),1)==tail(which(window.model$final.mod.coef[8:14] !=0),1)
            
          } else if (all(true.coef[8:14] ==0)==TRUE && all(window.model$final.mod.coef[8:14] ==0)==TRUE) {ma.order=FALSE}
          else if (all(true.coef[8:14] ==0)==TRUE && all(window.model$final.mod.coef[8:14] ==0)==FALSE){ma.order=FALSE}
          
         if (ma.order == TRUE && ar.order == TRUE){
           corect.model<- corect.model +1 
         }
          cat("AR order:", ar.order,"\n")
          cat("MA order:", ma.order,"\n")
        }  
        corect.model.rate <- corect.model/MC
        false.negative.rate[is.nan(false.negative.rate)]<-0
        false.positive.rate[is.nan(false.positive.rate)]<-0
        true.positive.rate[is.nan(true.positive.rate)]<-0
        cat("n",i," A", mean(true.positive.rate)," T", corect.model.rate, '-', mean(false.positive.rate), '+', mean(false.negative.rate), '\n')
        
      }  
    }
     
  }
}
