#'
#'  This code is adapted from concreg package
#'  Please refer to https://cran.r-project.org/web/packages/concreg/
#'
# concreg edited for our purpose
# credit -> package::concreg
concreg2 <- function
(
  formula=attr(data, "formula"),         # formula
  data=sys.parent(),                     #
  w_p_LOD,
  id=NULL,                               # identifier: numeric or character or factor
  normalize=TRUE,
  scale.weights=1,
  offset=NULL,
  alpha=0.05,                            # confidence limit
  maxit=50,                              # max. iterations
  maxhs=5,                               # half steps
  epsilon=1e-6,                          #
  maxstep=2.5,                           #
  x=TRUE,                                # for output
  y=TRUE,                                # for output
  print=TRUE,                            # print fitting information on screen
  c.risk=NULL,                           # Vector encoding competing risk: 0=censored, 1=event, 2=competing risk
  strata.var=NULL,                       # strata variable name
  trunc.weights=1,                       # quantile for weight truncation: all weights greater than that quantile will be truncated to that value
  npar=FALSE,                            # nonparametric estimation?
  # robust=FALSE,
  # jack=FALSE,
  ...                                    # N -> breslow; km -> prentice
)

{

  alpha.fp=c(0.20, 0.05)
  fp.iter=10                            # maximum number of iterations of large <fp> loop
  n <- nrow(data)

  ## generate or reorder id's such that values are within 1:n
  if (is.null(id)) { id <- 1:n } else id <- as.numeric(as.factor(id))
  maxid <- max(id)

  decomposeSurv <- function( formula, data, sort=FALSE, offset=NULL){
    ### decomposes complex survival formula
    ### trans: I(), powM2, powM1, powM0.5, sqrt, pow2, pow3, log,
    ### 2008-04

    orig.formula <- formula

    ## expand formula if needed:
    repeat {
      terms <- terms(formula, "fp", data=data)
      fac <- attr(terms, "factors")

      needed <- rownames(fac)[!rownames(fac) %in% colnames(fac)][-1]

      if(length(needed) == 0) break

      formula <- as.formula(paste(as.character(formula)[2], "~",
                                  as.character(formula)[3], "+",
                                  paste(needed, sep="+")))
    }

    ## determine position of <fp> terms (of all orders) in
    ## not expanded formula
    spec <- untangle.specials(terms, "fp", order=1:3)
    fpnames <- gsub("fp\\(([a-zA-Z0-9]*)\\)", "\\1", spec$vars) # c("a", "b")
    fpiden <- paste("(", fpnames, ")", sep="") # c("(a)", "(b)")

    formulaOrig <- formula
    ## replace fp() by the 8x2 possible terms (in formula[3])
    ## attention: fp-variables must consist of a-z,A-Z,0-9, otherwise
    ## must be made modifications here (gsub) ...
    ## FURTHER: 1st term is the linear one
    ## FURTHER: 1st half (=8) are the usual terms, 2nd half (9-16) are the repeated powers
    ## FURTHER: problems possible if one variable name is part of the other, e.g. ABC & ABCDE
    ##CODE <- paste("(I(\\1) + powM2(\\1) + powM1(\\1) + powM0.5(\\1) + log(\\1) + sqrt(\\1) + pow2(\\1) + pow3(\\1) + ",
    ##              "  RI(\\1) + RpowM2(\\1) + RpowM1(\\1) + RpowM0.5(\\1) + Rlog(\\1) + Rsqrt(\\1) + Rpow2(\\1) + Rpow3(\\1) )")
    CODE <- paste("(I(PT(\\1)) + powM2(PT(\\1)) + powM1(PT(\\1)) + powM0.5(PT(\\1)) + log(PT(\\1)) + ",
                  "sqrt(PT(\\1)) + pow2(PT(\\1)) + pow3(PT(\\1)) + ",
                  "RI(PT(\\1)) + RpowM2(PT(\\1)) + RpowM1(PT(\\1)) + RpowM0.5(PT(\\1)) + ",
                  "Rlog(PT(\\1)) + Rsqrt(PT(\\1)) + Rpow2(PT(\\1)) + Rpow3(PT(\\1)) )")
    sub3 <- gsub("fp\\(([a-zA-Z0-9]*)\\)", CODE, formulaOrig[3])
    formula <- as.formula(paste(as.character(formula)[2], "~", sub3))

    ## define simple transformations
    powM2 <- function(z) z^(-2)
    powM1 <- function(z) z^(-1)
    powM0.5 <- function(z) z^(-0.5)
    pow2 <- function(z) z^2
    pow3 <- function(z) z^3

    ## define repeated powers
    RI <- function(z) z * log(z)
    RpowM2 <- function(z) z^(-2) * log(z)
    RpowM1 <- function(z) z^(-1) * log(z)
    RpowM0.5 <- function(z) z^(-0.5) * log(z)
    Rlog <- function(z) log(z) * log(z)
    Rsqrt <- function(z) sqrt(z) * log(z)
    Rpow2 <- function(z) z^2 * log(z)
    Rpow3 <- function(z) z^3 * log(z)

    ## pretransformation function
    PT <- function(z) {
      obj <- fp.scale(z)
      (z + obj$shift) / obj$scale }

    ## construct 3-col response:
    resp <- model.extract(model.frame(formula, data = data), "response")
    if(is.null(dim(resp))) 	 resp<-cbind(resp-min(min(resp-0.01),0), rep(1,length(resp)))
    if(ncol(resp) == 2)      resp <- cbind(start=rep(0, nrow(resp)), resp)

    ## sortieren nach STOPzeit und -Cens
    if(sort) {
      sort <- order(resp[, 2],  -resp[, 3])
      data <- data[sort, , drop=FALSE]
      resp <- resp[sort, ]
    }

    mm <- model.matrix(formula, data = data) ## Model-Matrix
    mm1 <- mm[, -1, drop=FALSE]	# w/o intercept

    ## *** RETREIVE PRETRANS COEFS
    sub3 <- gsub("fp\\(([a-zA-Z0-9]*)\\)", "PT(\\1)", formulaOrig[3])
    formulaPT <- as.formula(paste(as.character(formulaOrig)[2], "~", sub3))
    ## varied pretransformation function to get coefficients:
    PT <- function(z) {
      obj <- fp.scale(z)
      c(shift=obj$shift, scale=obj$scale, rep(-999, length(z) - 2)) }
    mmPT <- model.matrix(formulaPT, data = data)
    PTcoefs <- mmPT[1:2, !is.na(mmPT[1,]) & (mmPT[3,]==-999), drop=FALSE]
    rownames(PTcoefs) <- c("shift", "scale")
    ## *** END

    ## offset ..
    if(length(offset) != 0)
      offset.values <- offset
    else
      offset.values <- NA

    terms <- terms(formula, "fp", data=data)
    fac <- attr(terms, "factors")
    labels <- attr(terms, "term.labels")

    ## splits by special chars
    f <- function(str)
      for(chars in c("(", ")", ":", " ", ",", "*", "^"))
        str <- unlist(strsplit(str, split=chars, fixed=TRUE))

    rowSplit <- sapply(rownames(fac), f, simplify=FALSE)	# splitted effects
    stopName <- tail(rowSplit[[1]], 2)[1]	# name of stoptime
    rowInter <- unlist(lapply(rowSplit[-1], function(z) any(z == stopName)))
    ##  rowOffset <- unlist(lapply(rowSplit[-1], function(z) any(z == "offset")))
    ##  rowOffset


    fac <- fac[-1, , drop=FALSE]	# omit Surv

    colSplit <- lapply(colnames(fac), f)
    colInter <- unlist(lapply(colSplit, function(z) any(z == stopName)))

    nTimes <- colSums(fac[rowInter, , drop=FALSE])
    nFac   <- colSums(fac[!rowInter, , drop=FALSE])

    inters <- (nFac>0) & (nTimes>0)
    NTDE <- sum(inters)


    timedata <- matrix(0, nrow(data), 0)
    timeind <- c()

    ## loop for (time x effect)
    for(i in which(inters)) {
      ## search pure time:
      ind <- (colSums(fac[rowInter, i] != fac[rowInter, , drop=FALSE]) == 0) & (nFac==0)
      timedata <- cbind(timedata, mm1[, ind, drop=FALSE])

      ## search pure effect:
      ind <- (colSums(fac[!rowInter, i] != fac[!rowInter, , drop=FALSE]) == 0) & (nTimes == 0)
      timeind <- c(timeind, which(ind[!colInter]))
    }
    mm1 <- mm1[, !colInter, drop=FALSE]

    covnames <- c(colnames(mm1),
                  paste(colnames(timedata), colnames(mm1)[timeind], sep=":")
    )
    alter <- c(colnames(mm1),
               paste(colnames(mm1)[timeind], colnames(timedata), sep=":")
    )

    ## indicator to identify the original formula:
    ind <- covnames %in% colnames(attr(terms(orig.formula, "fp", data=data), "factors")) |
      alter %in% colnames(attr(terms(orig.formula, "fp", data=data), "factors"))


    ## FP indicator matrix (1:16)
    NFP <- length(fpnames)
    fpind <- matrix(0, NFP, length(covnames), dimnames=
                      list(fpnames, covnames))
    for(i in seq(length=NFP)) {
      inds <- grep(fpiden[i], covnames)
      ## when interactions occur, 1:16,1:16 instead of 1:32 is needed
      fpind[i, inds] <- ((seq(along=inds) - 1) %% 16) + 1
    }

    ## return object
    list(
      fac=fac,                    # factor matrix ..
      resp=resp,                  # N x 3 - response matrix
      mm1=mm1,                    # model matrix without time effects

      NTDE=NTDE,                  # number time dep. effects
      timedata=timedata,          # matrix with time functions as columns
      timeind=timeind, 		       # indicator of time-dependend effects

      NFP=NFP,                    # number of frac.polys
      fpnames=fpnames,            # names of the variables used as fractional poly's
      fpind=fpind,                # matrix with frac.polyn. in each row, numbered by 1 to 5
      PTcoefs=PTcoefs,            # coefficients of pretransformations

      covnames=covnames,          # names of covariates
      ind=ind,			               # indicator vector:
      # which terms are really part of the formula
      offset.values=offset.values # offset values
    )
  }


  ## here only ONCE the full model matrix is spanned with all possible fp-effects
  obj.full <- decomposeSurv(formula, data, sort=FALSE, offset)

  # change Daniela Competing Risk
  if (is.null(c.risk)) {
    crisk <- obj.full$resp[,3]
    obj.full$resp <- cbind(obj.full$resp, crisk)
  }

  ## stratagruppen bilden
  if (is.null(strata.var)[1]) {
    obj.full$stratum <- rep(1, n)
  }

  max.strata <- max(obj.full$stratum)

  # change Daniela Competing Risk Weights (getrennt nach Strata)
  # geht nicht so einfach, Problem mit ties (gleiche Zeiten nur 1 mal in der Liste)
  G<-rep(1,nrow(obj.full$resp))     #### Georg 110517
  # new time variable for competing risk
  time.crisk <- obj.full$resp[,2]
  time.crisk[obj.full$resp[,"crisk"]==2] <- max(obj.full$resp[,2])+1
  obj.full$resp <- cbind(obj.full$resp, time.crisk)


  # calculate weights
  W <- w_p_LOD#rep(1, nrow(obj.full$resp) ) #concreg.wei(resp=obj.full$resp, max.strata=max.strata, stratum=obj.full$stratum, trunc.weights=trunc.weights)
  #browser()

  obj <- obj.full

  kk <- ncol(obj.full$mm1) # should be k2 - NTDE
  obj$mm1 <- obj.full$mm1[, obj$ind[1:kk], drop=FALSE]
  obj$covnames <- obj.full$covnames[obj$ind]

  obj$timeind  <- obj.full$timeind[obj$ind[-(1:kk)]]                      #! weg?
  obj$timedata <- obj.full$timedata[, obj$ind[-(1:kk)], drop=FALSE]       #! weg?
  ## re-index $timeind
  obj$timeind <- match(obj$timeind, (1:kk)[obj$ind[1:kk]])                #! weg?
  NTDE <- obj$NTDE <- length(obj$timeind)                                 #! weg?

  ind.offset <- sum(length(offset) != 0)

  k <- ncol(obj$mm1)    # number covariates w/o time-dep effects          # Anzahl Variablen
  k2 <- k + NTDE                                                          #! weg?
  if (npar) nonpar<-1
  else nonpar<-0

  PARMS <- c(n, k, 0, maxit, maxhs, maxstep, epsilon,  0, 0, 0, 0, 0, 0, max.strata, maxid, ind.offset, nonpar)

  # alles sortieren (obj, W, G, id ist sortiert)
  ord <- order(obj.full$stratum, obj.full$resp[,"time.crisk"],(1-obj.full$resp[,3]))
  obj$mm1     <- obj$mm1[ord,]
  obj$resp    <- obj$resp[ord,]
  obj$stratum <- obj$stratum[ord]
  G           <- G[ord]
  W           <- W[ord]
  id          <- id[ord]

  ## **************** fit model *********************
  concreg.fit2 <- function(obj, id, W, G, PARMS, npar){

    k <- PARMS[2]
    k2 <- k + obj$NTDE
    maxid <- PARMS[15]

    ## standardize model matrix, but only in semiparametric mode
    if(!npar) {
      sd1 <- apply(as.matrix(obj$mm1),2,sd)
      sd2 <- apply(as.matrix(obj$timedata),2,sd)
      Z.sd <- c(sd1, sd2 * sd1[obj$timeind])
      ZxZ <- as.matrix(Z.sd) %*% t(as.matrix(Z.sd)) # used often to restandardize ...
      obj$mm1 <- scale(obj$mm1, FALSE, sd1)
    } else {
      sd1<-1
      sd2 <- 1
      Z.sd<-1
      ZxZ <-1
    }

    ##if(ind.offset)
    obj$mm1o <- if(PARMS[16] != 0) cbind(obj$offset.values, obj$mm1) else obj$mm1

    CARDS <- cbind(obj$mm1o, obj$resp[,c(2, 4)], W, G, id, obj$stratum)

    if(!npar) obj$timedata <- scale(obj$timedata, FALSE, sd2)                         # weg?
    mmm <- cbind(obj$mm1, obj$timedata) # model matrix inc. time data       #! ohne timedata?

    DFBETA <- matrix(0, maxid, k2)                                          #! k2 durch k ersetzen?
    IOARRAY <- rbind(rep(1, k2), matrix(0, 2+2*k2, k2))                     #! k2 durch k ersetzen?
    if(obj$NTDE >0)                                                         #! weg?
      IOARRAY[4, (k+1):k2] <- obj$timeind                                   #! weg?

    ## --------------- Aufruf Fortran-Routine ????????? ----------
    storage.mode(CARDS) <- storage.mode(PARMS) <- storage.mode(IOARRAY) <- storage.mode(DFBETA) <- "double"
    #       dyn.load("concreg_dll.dll")
    #        dyn.load("D:\\WORK\\concreg_dll.dll")
    value <- .Fortran("concregfit", #"CONCREG",                                            #! anpassen
                      cards=CARDS,
                      outpar = PARMS,
                      outtab = IOARRAY)
    #                          PACKAGE=concreg.fp)
    # browser()
    if(value$outpar[8])
      warning("Error in Fortran routine concreg; parms8 <> 0")

    coefs <- value$coefs / Z.sd                                             #! value$coefs gibt es nicht
    cov.mb <- matrix(value$outtab[4:(k2+3), ], ncol=k2) / ZxZ            # "model-based" varianz
    cov.rob <- matrix(value$outtab[(3+k2+1):(3+2*k2), ], ncol=k2) / ZxZ  # "robuste" varianz (standard)

    res <- list(                                                            #! anpassen
      cards=value$cards,
      outpar=value$outpar,
      outtab=matrix(value$outtab, nrow=3+2*k2),

      coef.orig=value$outtab[3,  ],
      coefs=value$outtab[3,  ] / Z.sd, # coefficients
      cov.rob=cov.rob,                 # covariances
      cov.mb=cov.mb,
      Z.sd=Z.sd,
      ZxZ=ZxZ,                                                    #! weg?
      mmm=mmm                          # model matrix
    )
    res
  }

  value0 <- concreg.fit2(obj=obj, id=id, W=W, G=G, PARMS=PARMS, npar=npar)            #! Aufruf der Funktion mit dem Fortran Aufruf

  vars <- as.matrix(diag(value0$cov.rob))

  ## DECIDE THE MODEL: USE PVALS ##############
  probs <- 1 - pchisq((value0$coefs^2/vars), 1)

  ## ########## NOW ONLY FINAL MODEL IS CONSIDERED ############
  names(value0$coefs) <- obj$covnames
  if(value0$outpar[10]>=maxit)
    cat("No convergence attained in ", value0$outpar[10], " iterations.\n", sep="")

  Means <- colMeans(value0$mmm)

  ## return object
  fit <- list(coefficients = value0$coefs,     # coefficients of the fit
              cards    = value0$cards,         #
              parms    = value0$outpar,
              ioarray  = value0$outtab,
              loglik   = value0$outpar[12:11],
              #                    dfbeta.resid = dfbeta.resid,
              alpha    = alpha,                # significance level
              var      = value0$cov.rob,                 # covariance matrix
              df       = k2,                   # degrees of freedom (k + NTDE)
              iter     = value0$outpar[10],
              method.ties = "no",         #
              n = n,                           # number observations

              y = obj$resp,                    # responses
              formula = formula,               # original formula
              exit.code=value0$outpar[8],

              call    = match.call(),
              cov.mb  = value0$cov.mb,

              Wald =  (t(value0$coefs) %*% solve(value0$cov.rob)) %*% value0$coefs,
              means   = Means,               # means of <X> (model matrix)
              linear.predictors= as.vector(scale(value0$mmm, Means, scale=FALSE) %*% value0$coefs),
              method  = "Weighted Estimation",
              method.ci= "Wald",             #
              ci.lower= exp(value0$coefs + qnorm(alpha/2) * vars^0.5), #
              ci.upper= exp(value0$coefs + qnorm(1 - alpha/2) * vars^0.5), #
              prob    = probs,               # p-values
              G       = G,
              #                    W       = W,
              W       = cbind(obj$stratum,obj$resp[,2],W)[obj$resp[,3]==1,],
              offset.values= obj$offset.values, #
              x       = if(x) obj$mm1 else NA,   # return original model matrix if requested
              npar = npar
  )

  names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- obj$covnames
  attr(fit, "class") <- c("concreg")
  fit
}










