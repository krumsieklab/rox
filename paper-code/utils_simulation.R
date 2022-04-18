# induces LOD-based missingness, given percentage
fnaify <- function(x, perc = 0.3, rand =F){
  inds <- 
    if(!rand) order(x)[  seq(length(x)*perc) ] 
  else sample( length(x)*perc )
  x[inds] = NA
  x
}

# induce missingness with probabilistic LOD modeled by the parameter pLOD
# pLOD: defines the LODness as follow
# a: (area under the sigmoid curve for the region x = [0,4], y = [0.5,1]) 
# b: (area above the sigmoid curve for the region x = [0,4], y = [0.5,1])
# pLOD = a/(a+b)
fnaify_plod <- function(x, pLOD, d=0){
  # find plogis scale for pLOD
  si = find_plogis_scale(pLOD)$par
  
  x =  as.numeric(scale(x)) # needed for plogis  
  p_na =  plogis(-x+d/2, scale = si)
  
  # realization of missingness based on px
  isNA =  sapply( p_na, function(prob) rbinom(1,1,prob))==1
  # missing values coded as NA 
  xh = x
  xh[isNA]=NA
  
  xh
}

# finds the sigmoid scale parameter for given plod paramrter 
# which models the LODness
find_plogis_scale <- function( plod,q = 4){
  if(plod==1) return(list (par = 1e-10) )
  if(plod==0) return(list( par = Inf) )
  
  plod = (plod+1)/2 # to convert into 0.5-1 scale
  
  ff <- function(a){
    ph = integrate(function(x)plogis(-x,scale = a), lower = -q, upper = 0)$value/q
    abs(ph-plod)
  }
  optim(0, ff, method = "Brent", lower = 1e-10, upper = 100)
}

# finds the noise_factor given concordance, d
#   where y = x + noise_factor*z
#         d = concordance(x~y)
#         {x, z} ~ N(0,1)
#
# if fast = TRUE, analytical approximation used 
#         = FALSE, it is found emprically 
find_noise_factor <- function(d, fast = TRUE, N = 1e4){
  if(fast){
    # returns the beta, b, where y = x + b*z
    # given that r = cor(x, y), and {x, z} ~ N(0,1)
    r2beta <- function(r) sqrt(1/(r^2) -1) 
    
    # returns the pearson r given the concordane
    # based on Greiner relations
    # assume there is no ties
    d2r <- function(d) sin(pi/2 *(2*d-1))
    
    # combine d2r and r2beta to implement d2beta 
    d2beta <- function(d) r2beta(d2r(d))
    
    return(d2beta(d))
  }else{
    if(is.null(x) & is.null(z)){
      x = rnorm(N)
      z = rnorm(N)
    }
    # # pick the parameter such that desired concordance is obtained: 
    # #     d = concordance(y~x), where y = x + nfs*z
    # # noise factors
    # nfs = seq(sqrt(0.5),sqrt(6.5),length.out = 601)^2
    # 
    # or using optim function
    
    ff <- function(b){
      y = z*b + x
      abs(d-concordance(x~y)$concordance)
    }
    optimize( ff, interval = c(0.5, 7), tol = 0.01)$minimum
  }
}
