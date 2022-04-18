#' rox: estimates concordance based on (R)ank (O)rder with (X)missing values
#'
#' @param xh numeric. A vector of measurements that might include missing values, such as a metabolite measurement where missing values are due to left-truncation of the distribution, e.g. because of a Limit of Detection (LOD)
#' @param y numeric. A variable of interest with no missing values, such as age, gender, or disease status
#' @param switching logical. TRUE by default. Indicates whether, if no evidence in support of an LOD effect is found, switching to complete case analysis is allowed
#' @param w numeric (optional). Value of the weight to be assigned to missing values. By default it is defined as: w = 1 - (missingness ratio), as described in Buyukozkan et al.
#'
#' @return estimated concordance index and its corresponding standard error and p-value
#' @export
#'
#' @examples
#' # sample size
#' n = 300
#' set.seed(42)
#' y = rnorm(n)
#' x =  y + rnorm(n)
#'
#' # induce missing values
#' xh = x
#' xh[xh<0] = NA
#'
#' # rox with missing values
#' rox::rox(xh, y)
#' # no missing values
#' rox::rox(x, y)
#'
rox <- function(xh, y, switching = T, w = NULL, ...){ # cvar: based on null
  # get the rox test
  # to do: clean rox-cox fit from environment for memory efficiency
  exargs = list(...)
  fits = NULL

  # no missing value return the fit
  if(!any(is.na(xh))){
    # write summary printer
    dfit = concordance(xh~y)
  }else{
    if(is.null(w)){ # estimate plod
      nod = f_estiweight(xh, y, switching)
      fits = nod[-length(nod)]
      w = nod$w
      S = nod$S
      # fit dS for its variance
      dfit = if(w >0) concordance(S~y, weights = c(w,1)[S[,2]+1], reverse = T) else nod$fits$d1_fit
    }else{
      stopifnot( length(w)==1 & w <= 1 & w > 0 )
      S = f_NA2Surv(xh, Inf, "right")
      dfit = concordance(S~y, weights = c(w,1)[S[,2]+1], reverse = T)
    }
  }

  re = list()
  re$stats = summary_dfit(dfit, se = exargs$se)
  re$w = w
  re$dfit = dfit
  re$dfit$call <-NULL

  # if details to be returned
  if(!is.null(exargs$details) & !is.null(fits$fits))
    re$details = c(fits$fits, list(dS = list(concordance = unname(fits$dS))))

  # write print function for this
  structure( re, class = "rox" )
}

print.rox <- function(obj){

  sm = unlist(obj$stats)
  cat("\n ", paste(paste(c("  ","  ","    ", "   "), names(sm)),collapse = " "))
  cat("\n")
  cat(" ", c(round(sm[-4],4), format.pval(sm[4], digits = 2)))
  cat("\n")
  cat("---\n")
  if(!is.null(obj$w)) cat("w =", obj$w,"\n") else cat("no missing value.\n")
  if(!is.null(obj$details)){
    cat("\n")
    print( round(sapply(obj$details, `[[`, "concordance"),3) )
  }
  cat("\n")

  invisible(NULL)

}
