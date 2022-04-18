
#' multivaariable rox modeling
#'
#' @param xh numeric. A vector of measurements that might include missing values, such as a metabolite measurement where missing values are due to left-truncation of the distribution, e.g. because of a Limit of Detection (LOD)
#' @param ys data.frame with numeric variables. A data.frame where columns are the variables of interest with no missing values, such as age, gender, or disease status
#' @param switching logical. TRUE by default. Indicates whether, if no evidence in support of an LOD effect is found, switching to complete case analysis is allowed
#' @param w numeric (optional). Value of the weight to be assigned to missing values. By default it is defined as: w = 1 - (missingness ratio), as described in Buyukozkan et al.
#' @param func function. This defines how to model the multivariable relation. The default is 'concreg', which uses an exponential link function to model concordance probabilities.
#' Another option is 'coxph', which uses Cox's partial likelihood and provides an approximate solution to 'concreg' but is faster for big datasets.
#'
#' @return results and details of the fit
#' @export
#'
#' @examples
#'
#'  # sample size for simulation
#'  n = 500
#'
#'  set.seed(42)
#'  y = scale(rnorm(n))
#'  z = scale(rnorm(n))
#'  x = y + 2*z + rnorm(n)
#'
#'  # induce missing values
#'  xh = x
#'  xh[ xh < -0.75] = NA
#'
#'  # with missing values
#'  rox::rox_mv(xh, data.frame(z,y), func = "concreg")
#'  # no missing values
#'  rox::rox_mv(x, data.frame(z,y), func = "concreg")
#'
rox_mv <- function(xh, ys, switching = T, w = NULL, func = c("concreg", "coxph"),...){ # cvar: based on null
  # get the rox fit
  # to do: clean rox-cox fit from environment for memory efficiency

  exargs = list(...)
  func = match.arg(func)
  fits = NULL

  if(!is.data.frame(ys)) stop("ys is not data.frame, use rox instead!")

  # # prepare the data
  # does not matter how surv object is created
  S = f_NA2Surv(xh, Inf, "right")
  # do matrix
  ysh <- f_deparse(ys)
  # keep colnames
  cnames = colnames(ysh)

  # check if var =0 for any of the column
  sdi = apply(ysh, 2, function(x) abs(max(x) - min(x)) < .Machine$double.eps ^ 0.5)
  if(sum(sdi) == ncol(ysh)) stop("variance of outcome is 0!")
  else if(any(sdi)) warning("variance of some of the outcomes are 0!")
  ysh = ysh[, !sdi, drop =F]
  # # end of data preparation

  #---estimate w
  fit0 = NULL
  if(!any(is.na(xh))){ # no missing value return the fit
    w = 1
  }else{
    if(is.null(w)){ # estimate plod
      fit0 = get_fit(S, ysh, sum(!is.na(xh))/length(xh), func = func) #-sum(is.na(xh))/length(xh)
      nod = f_estiweight(xh,  ysh %*% (fit0[,"B"]), switching)
      fits = nod[-length(nod)]
      w = nod$w
    }else{
      stopifnot( length(w)==1 & w <= 1 & w > 0 )
    }
  }
  fit1 = get_fit(S, ysh, w, func = func)

  # summary of multivariate fit, fix the names
  sm = matrix(NA, length(cnames), 4, dimnames = list(cnames, c( "B", "se","z", "p")))
  sm[colnames(ysh), ] = fit1

  # overall fit
  ws = c(w,1)[S[,2]+1]
  B = fit1[,"B"]
  B[is.na(B)] = 0
  dfit = concordance(S~yh, reverse = T,
                     data = data.frame(S=S, yh = ysh %*% B)[ws>0,],
                     weights = ws[ws>0])
  dfit$call <-NULL

  re = list()
  re$stats_overall = summary_dfit(dfit, se = exargs$se)
  re$stats_model = sm
  re$func = func
  re$w = if(any(is.na(xh))) w else NULL
  re$dfit_overall = dfit

  # if details to be returned
  if(!is.null(exargs$details) & !is.null(fits$fits))
    re$details = c(fits$fits, list(dS = list(concordance = unname(fits$dS)), m0func = fit0))

  # write print function for this
  structure( re, class = "roxfit" )

}

print.roxfit <- function(obj){

  # summary of model fit
  ndigit = options()$digits
  options(digits = 4)
  cat("\n(coefficients for multivariable model) \n")
  cat("model:", obj$func,"\n")
  print(as.data.frame(obj$stats_model))
  options(digits = ndigit)
  cat("---------------------------------------\n")
  # summary of model fit


  sm = unlist(obj$stats_overall)
  cat("\n ", paste(paste(c("  ","  ","    ", "   "), names(sm)),collapse = " "))
  cat("\n")
  cat(" ", c(round(sm[-4],4), format.pval(sm[4], digits = 2)))
  cat("\n")
  cat("---\n")
  if(!is.null(obj$w)) cat("w =", obj$w,"\n") else cat("no missing value.\n")


  if(!is.null(obj$details)){
    cat("\n") # last one is just coefs of fit0 not concordance
    print( round(sapply(obj$details[-length(obj$details)], `[[`, "concordance"),3) )
  }
  cat("\n")

  invisible(NULL)
}


# understand how we can run fortran code of other package without including
# into the package, and implement for concreg



