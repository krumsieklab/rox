# deparse data to data matrix
f_deparse <- function(ys){
  if(is.vector(ys)){
    if(!is.factor(ys)) cbind(ys) 
    else model.matrix(~ys)[,-1,drop=F]
  }else{
    model.matrix(~., data.frame(ys))[,-1,drop=F]
  }
}

get_fit <- function(S, ym, w, func = c("concreg", "coxph")){
  func = match.arg(func)
  ws = c(w, 1)[S[,2]+1]
  sm <-
    if(func == "concreg"){
      # for concreg it makes difference
      S1r = f_tied_censored(S)
      fit = concreg2( Surv(time,event = event)~., data.frame(time = S1r[,1], event = S1r[,2], ym)[ws>0,], w_p_LOD = ws[ws>0])
      cbind(B = -fit$coefficients, se = sqrt(diag(fit$var)))
      # # returning these summary stats are possible
      # # this is a note for future improvements 
      # list( loglik = fit$loglik, Wald = fit$Wald[1] )
      
    }else{
      fit = coxph( S~., data.frame(S=S, ym)[ws>0,], weights = ws[ws>0])
      cbind(B = fit$coefficients, se = sqrt(diag(fit$var))) 
      # # returning these summary stats are possible
      # # this is a note for future improvements 
      # list( loglik = fit$loglik, Wald = fit$wald.test[1], score = fit$score, rscore = fit$rscore)
    }
  sm = cbind(sm, z =sm[,1]/sm[,2])
  
  structure( cbind( sm, p = unname(-2*expm1(pnorm(sm[,"z"], log.p = T)))), f = func )
}

# convert right censored to tied events
f_tied_censored <- function(S){
  S1r = S
  S1r[,1] = max(S1r[,1])+1-S1r[,1]
  S1r[,2] = 1
  S1r
}

