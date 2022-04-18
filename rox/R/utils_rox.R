# test new implementation
f_estiweight <- function(xh, y, switching = T){
  if(!any(is.na(xh))) return(NA)

  S = f_NA2Surv(xh, Inf, "right")

  w = 1- sum(is.na(xh))/length(xh)
  if(!switching)
    return(
      list(w =  w, fits = NULL, dS = NULL, S= S)
    )

  dy_fit = concordance(S[,2]~y)
  d1_fit = concordance(S[S[,2]==1,1]~y[S[,2]==1], reverse = T)

  gong = seesaw(d1_fit, dy_fit)

  # number of pairs and corresponding ds
  n1 = sum(d1_fit$count[1:3])
  ny = sum(dy_fit$count[1:3])
  d1 = d1_fit$concordance
  dy = dy_fit$concordance
  if(is.na(d1)) d1 = 0.5

  # calculate final possible d
  dS = (n1*d1 + ny*w*dy)/(n1+ny*w)

  # for direction to be mattered
  opt2 = (1-dS) > d1
  opt2 = T
  flag = if( (dS < 0.5) & opt2 ) (1-dy) >= (1-d1) else  dy >= d1

  w =if(flag | dS==0.5 ) w else 0
    
  # for speciall cases 
  #if( (gong==0 & w !=0 ) | (gong!=0 & w ==0) ){
  #  #warning("Estimation might be biased!", call. = F)
  #}
  #
  #if(gong>0 & switching){
  #  #w = 0
  #  #warning("Missigness might be due to right censoring, so estimation might be biased!", call. = F)
  #}
  #

  list(w =  w,
       fits = list(d1_fit = d1_fit, dy_fit = dy_fit),
       dS = dS, S= S)
}

# function to create desired Surv object
f_NA2Surv<- function(xh, direction = c(Inf, -Inf), type = c("right", "count")){
  direction = match.arg(type)
  type = match.arg(type)

  inds = is.na(xh)
  if(direction > 0){
    x = max(xh,na.rm = T) - xh + 1
    x[inds] = max(x,na.rm = T) + 1 # add at the end
    S <- if(type=="right") Surv(x, !inds) else Surv(x*0+1, x+1, !inds)
  }else{
    x = xh - min(xh,na.rm = T) + 1 + sum(inds) # lag for adding at the beginning
    t1 = x*0 + 0.5
    t2 = x+1
    t1[inds] = seq(sum(inds))-0.5
    t2[inds] = seq(sum(inds))
    S = Surv(t1, t2, t1>0)
  }
  return(S)
}

# returns d and correponding p-val
summary_dfit <- function(d, se = NULL){
  # se = "var" or "cvar" or harmonized
  se<- sqrt( if(is.null(se)) (d$cvar/2 + d$var/2) else d[[se]] )
  # # --------------
  z =  (d$concordance - 0.5)/se
  # be careful p value is biased and not exact
  p = -expm1( pnorm( abs(z), log.p = T) )
  p = min(1, 2*p) # one sided p-value?
  c(d = d$concordance, se = se, z=z, p=p)
}

# warn if there is an anti-LOD
seesaw <- function(d1, d0){
  d1 = d1$concordance - 2*sign(d0$concordance - 0.5)*sqrt(d1$cvar+d1$var)
  d0 = d0$concordance

  if(max(d0, 1-d0) > max(d1, 1-d1)){
    if( d1==0.5 ) return(-2*(d0>0.5) +1)
    if( (d0>0.5) == (d1>0.5) ) return(-1)
    return(1)
  }
  return(0)
}



