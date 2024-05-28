#' Spatially-dependent normalisation for spatial transcriptomics datas
#'
#' Performs normalisation of spatial transcriptomics data using spatially-dependent spot-specific size factor.
#'
#' @param spe a SpatialExperiment or Seurat object, with the count data stored in 'counts' or 'data' assays respectively.
#' @param gene.model a character, specifying the model to use for gene/protein abundances (default 'nb'). This should be 'nb' for count based datasets.
#' @param sample.p a numeric, specifying the (maximum) proportion of cells/spots to sample for model fitting (default is 0.25).
#' @param scale.factor a numeric, specifying the sample-specific scaling factor to scale the adjusted count.
#' @param df.tps a numeric, specifying the degrees of freedom for the thin-plate spline (default is 6).
#' @param lambda.a a numeric, specifying the smoothing parameter for regularizing regression coefficients (default is 0.0001). Actual lambda.a used is lambda.a * ncol(spe).
#' @param step.factor a numeric, specifying the multiplicative factor to decrease IRLS step by when log-likelihood diverges (default is 0.5).
#' @param inner.maxit a numeric, specifying the maximum number of IRLS iteration for estimating mean parameters for a given dispersion parameter (default is 50).
#' @param outer.maxit a numeric, specifying the maximum number of IRLS iteration (default is 25).
#' @param ... other parameters to pass to SpaNorm.

#' @return a SpatialExperiment or Seurat object with the adjusted data stored in 'logcounts' or ''
#' @export
setGeneric("SpaNorm", function(
    spe,
    gene.model = c("nb"),
    sample.p = 0.25,
    scale.factor = 1,
    df.tps = 6,
    lambda.a = 0.0001,
    ...) {
  standardGeneric("SpaNorm")
})

#' @rdname filterGenes
setMethod(
  "SpaNorm",
  signature("SpatialExperiment"),
  function(spe, gene.model, ...) {
    checkSPE(spe)
    gene.model = match.arg(gene.model)

    # extract counts and coords
    emat = SummarizedExperiment::assay(spe, "counts")
    coords = SpatialExperiment::spatialCoords(spe)

    # fit SpaNorm model
    fit.spanorm = SpaNorm_intl(emat, coords, gene.model, ...)

    # normalise using logPAC

    return(spe)
  }
)

sampleRandom <- function(coords, nsub) {
  idx = rep(FALSE, nrow(coords))
  idx[sample.int(nrow(coords), size = nsub)] = TRUE
  return(idx)
}

normaliseLogPAC <- function(Y, fit.spanorm) {
}

SpaNorm_intl <- function(Y, coords, gene.model, sample.p = 0.25,  scale.factor = 1, df.tps = 6, lambda.a = 0.0001, ...) {
  # parameter checks
  if (sample.p <= 0 | sample.p > 1) {
    stop("'sample.p' should be in the interval (0,1]")
  }
  if (lambda.a <= 0) {
    stop("'lambda.a' should be greater than 0")
  }

  # setting-up variables for the NB regression models
  lambda.a = lambda.a * ncol(Y)
  cl = scran::quickCluster(Y)
  logLS = log(pmax(1e-08, scran::calculateSumFactors(Y, clusters = cl)))
  bs.xy = bs.tps(coords[, 1], coords[, 2], df.tps = 6) # get basis for the thin-plate spline
  W = model.matrix(~ logLS * bs.xy)[, -1]

  # sample data for computational efficiency
  nsub = round(sample.p * ncol(Y))
  message(sprintf("%d cells/spots being sampled", nsub))
  if (nsub > 3000) {
    warning(sprintf("consider reducing 'sample.p' to %.2f to increase computational efficiency", floor(3000 / ncol(Y) * 100) / 100))
  } else if (nsub == 0) {
    stop(sprintf("'sample.p' is too small, consider using %.2f", min(1, floor(3000 / ncol(Y) * 100) / 100)))
  }
  # random sampling
  idx = sampleRandom(coords, nsub)

  # fit model
  if (gene.model == "nb") {
    fit.spanorm = SpaNormNB_intl(Y, W, idx, lambda.a, ...)
  }
}

SpaNormNB_intl <- function(Y, W, idx, lambda.a, step.factor = 0.5, inner.maxit = 50, outer.maxit = 25) {
  # parameter checks
  if (step.factor <= 0 | step.factor >= 1) {
    stop("'step.factor' should be in the interval (0,1)")
  }
  if (inner.maxit <= 0 | outer.maxit <= 0) {
    stop("'inner.maxit' and 'outer.maxit' should be greater than 0")
  }
  if (inner.maxit - as.integer(inner.maxit) != 0 | outer.maxit - as.integer(outer.maxit) != 0) {
    stop("'inner.maxit' and 'outer.maxit' should be integers")
  }

  # subset data points used for model fitting
  Ysub = as.matrix(Y[, idx, drop = FALSE])
  Wsub = W[idx, , drop = FALSE]
  nW = ncol(Wsub)
  nsub = sum(idx)

  # setting initial regression parameter estimates for all genes
  gmean = rowMeans(log(Ysub + 1)) # rowMeans(Z)
  alpha = matrix(0, nrow(Ysub), ncol(W))
  # alpha for logLS
  alpha[, 1] = 1

  # subset to a maximum of 50 cells/spots for dispersion parameter estimation (for speed-up)
  nsub.psi = min(50, nsub)
  psi.idx = sample.int(nsub, size = nsub.psi)
  Ysub.psi = Ysub[, psi.idx, drop = FALSE]
  psi = rep(0, nrow(Ysub))
  
  # convergence trackers
  conv = FALSE
  logl.psi = c()
  iter = 1

  # outer loop of IRLS for estimating disp parameter
  while (!conv) {
    # calculate/update dispersion parameter (psi) estimates for all genes
    Wa = Matrix::tcrossprod(alpha, Wsub) # offsets
    offs.psi = Wa[, psi.idx, drop = FALSE]
    psi.tmp = tryCatch(
      {
        edgeR::estimateDisp(Ysub.psi, as.matrix(rep(1, nsub.psi)), offset = offs.psi, tagwise = TRUE, robust = TRUE)$tagwise.dispersion
      },
      error = function(e) {
        rep(NA, length(psi))
      }
    )
    valid.psi = !is.na(psi.tmp) & !is.infinite(psi.tmp)
    psi[valid.psi] = psi.tmp[valid.psi]

    # calculate initial loglik for this iteration
    loglik = colSums(dnbinom(Ysub, mu = exp(gmean + Wa), size = 1 / psi, log = TRUE))

    # fit NB given dispersion estimates (psi) and extract required components
    fit.nb = fitNBGivenPsi(Ysub, Wsub, gmean, alpha, psi, lambda.a, step.factor, inner.maxit, loglik)
    Wsub = fit.nb$Wsub
    gmean = fit.nb$gmean
    alpha = fit.nb$alpha
    loglik = fit.nb$loglik

    # check convergence
    logl.psi = c(sum(loglik), logl.psi) # push to stack
    conv.logl = FALSE
    iter = iter + 1 # similar post-increment to fitNBGivenPsi
    if (iter > 2) { # should be 2 because iter is post-incremented
      conv.logl = ifelse((logl.psi[1] - logl.psi[2]) / abs(logl.psi[2]) < 1e-04, TRUE, FALSE)
    }
    conv = conv.logl | iter > outer.maxit
  }

  # final values
  fit.spanorm = list(
    gmean = gmean,
    alpha = alpha,
    psi = psi,
    loglik = loglik,
    loglik.iter = logl.psi
  )

  return(fit.spanorm)

  # # calculate mu
  # mu = gmean + Matrix::tcrossprod(alpha, W) # log(mu)
  # lmu.max = matrixStats::rowMedians(mu) + 4 * matrixStats::rowMads(mu)
  # mu = exp(pmin(mu, lmu.max)) # winsorize

  # # calculate mu without library size effect
  # mu.2 = gmean + Matrix::tcrossprod(alpha[, -1], W[, -1]) # log(mu)
  # lmu.2.max = matrixStats::rowMedians(mu.2) + 4 * matrixStats::rowMads(mu.2)
  # mu.2 = exp(pmin(mu.2, lmu.2.max)) # winsorize

  # # winsorize dispersion parameters (large disp slow the process)
  # psi.max = exp(median(log(psi)) + 3 * mad(log(psi)))
  # psi = pmin(psi, psi.max)
}

fitNBGivenPsi <- function(Ysub, Wsub, gmean, alpha, psi, lambda.a, step.factor, maxit, loglik = NULL) {
  stopifnot(ncol(Ysub) == nrow(Wsub))
  stopifnot(nrow(Ysub) == length(gmean))
  stopifnot(ncol(Wsub) == ncol(alpha))
  stopifnot(nrow(Ysub) == nrow(alpha))
  stopifnot(nrow(Ysub) == length(psi))

  if (is.null(loglik))  {
    # if not pre-computed (e.g., by outer loop), compute
    # helpful pattern for external use of the fitting function
    lmu.hat = gmean + Matrix::tcrossprod(alpha, Wsub)
    sig.inv = 1 / (exp(-lmu.hat) + outer(psi, rep(1, nsub), "/"))
    loglik = colSums(dnbinom(Ysub, mu = exp(lmu.hat), size = 1 / psi, log = TRUE))
  }

  # set the IRLS step size to 1
  nsub = ncol(Ysub)
  nW = ncol(Wsub)
  step = rep(1, nsub)

  # convergence trackers
  conv = FALSE
  logl.beta = c() # stack of loglik values - push values
  iter = 1

  # initialise halving counters
  halving = 0

  while (!conv) {
    # saving new best estimate (in case we need to go back in our search)
    best.gmean = gmean
    best.a = alpha

    # calc working vector Z for a subset of cells
    lmu.hat = gmean + Matrix::tcrossprod(alpha, Wsub)
    sig.inv = 1 / (exp(-lmu.hat) + outer(psi, rep(1, nsub), "/"))
    Z = lmu.hat + sweep(((Ysub + 0.01) / (exp(lmu.hat) + 0.01) - 1), 2, step, "*")

    # store current alpha est
    alpha.old = alpha

    # update alpha for all genes
    wt.cell = matrixStats::colMeans2(sig.inv)
    # prevent outlier large weight
    wt.cell = pmin(wt.cell, quantile(wt.cell, probs = 0.98))
    b = (Z - gmean) %*% diag(wt.cell) %*% Wsub
    alpha = b %*% Matrix::chol2inv(chol(Matrix::crossprod(Wsub * wt.cell, Wsub)))
    # set first column of alpha to be the same for all genes (see Model specification)
    alpha[, 1] = mean(alpha[, 1])
    b = (Z - gmean - Matrix::tcrossprod(alpha[, 1, drop = FALSE], Wsub[, 1, drop = FALSE])) %*% diag(wt.cell) %*% Wsub[, -1]
    alpha[, -1] = b %*% Matrix::chol2inv(chol(Matrix::crossprod(Wsub[, -1] * wt.cell, Wsub[, -1]) + lambda.a * diag(nW - 1)))

    # if new alpha est has missing/inf values, revert to previous estimate
    if (any(is.na(alpha) | is.infinite(alpha))) alpha <- alpha.old

    # reduce outliers
    a.med = matrixStats::colMedians(alpha)
    a.mad = matrixStats::colMads(alpha)
    a.max = a.med + 4 * a.mad
    a.min = a.med - 4 * a.mad
    alpha = pmin(t(alpha), a.max) # +ve outliers, orig: t(pmin(t(alpha), a.max))
    alpha = t(pmax(alpha, a.min)) # -ve outliers, orig: t(pmax(t(alpha), a.min))

    # step2: update gmean
    Z.res = Z - Matrix::tcrossprod(alpha, Wsub)
    gmean = matrixStats::rowSums2(Z.res * sig.inv) / matrixStats::rowSums2(sig.inv)

    # calculate current logl
    lmu.hat = gmean + Matrix::tcrossprod(alpha, Wsub)
    loglik.tmp = colSums(dnbinom(Ysub, mu = exp(lmu.hat), size = 1 / psi, log = TRUE))

    # check degenerate case
    degener = sum(loglik) > sum(loglik.tmp)
    if (is.na(degener) | degener) {
      # alter stepsize of genes
      check.gene = loglik > loglik.tmp
      step[check.gene] = step[check.gene] * step.factor
      # revert to previous best
      gmean = best.gmean
      alpha = best.a
      # increase halving counter
      halving = halving + 1
      if (halving >= 3) {
        loglik.tmp = loglik
        degener = !degener
      }
    }

    # check degener (again in case halving exceeds maximum, i.e., 3)
    if (!degener) {
      loglik = loglik.tmp
      logl.beta = c(sum(loglik), logl.beta) # push to stack
      iter = iter + 1
      halving = 0 # reset
    }

    # check convergence
    conv.logl = FALSE
    if (iter > 2) { # should be 2 because iter is post-incremented
      conv.logl = ifelse((logl.beta[1] - logl.beta[2]) / abs(logl.beta[2]) < 1e-04, TRUE, FALSE)
    }
    conv = conv.logl | iter > maxit
  }

  # final values
  fit.nb = list(
    gmean = gmean,
    alpha = alpha,
    loglik = loglik,
    loglik.iter = logl.beta
  )

  return(fit.nb)
}

bs.tps <- function(x, y, df.tps = 6) {
  # checks
  if (df.tps <= 0) {
    stop("'df.tps' should be greater than 0")
  }
  if (df.tps - as.integer(df.tps) != 0) {
    stop("'df.tps' should be an integer")
  }

  bs.x = splines::ns(x, df = df.tps)
  bs.y = splines::ns(y, df = df.tps)
  bs.xy = matrix(0, nrow = length(x), ncol = df.tps ^ 2)
  for (i in 1:df.tps) {
    for (j in 1:df.tps) {
      bs.xy[, (i - 1) * ncol(bs.x) + j] <- bs.x[, i] * bs.y[, j]
    }
  }
  bs.xy = scale(bs.xy, scale = FALSE)
  return(bs.xy)
}

SpaNorm_old <- function(spe,pCells.touse=0.25,lambda.a=0.0001,scale.factor=1,step.fac=0.5,inner.maxit=50,outer.maxit=25) {
  lambda.a <- lambda.a*ncol(spe)
  cl <- scran::quickCluster(spe)
  spe <- scran::computeSumFactors(spe, clusters = cl)
  spe$sizeFactor <- pmax(1e-08,spe$sizeFactor)
  spe$logLS <- log(spe$sizeFactor)
  x <- spatialCoords(spe)[,1]
  y <- spatialCoords(spe)[,2]
  bs.x <- splines::ns(x,df=6)
  bs.y <- splines::ns(y,df=6)
  bs.xy<- matrix(0,nrow=nrow(bs.x),ncol=ncol(bs.x)*ncol(bs.y))
  for(i in 1:ncol(bs.x)) {
    for(j in 1:ncol(bs.y)) {
      bs.xy[,(i-1)*ncol(bs.x)+j] <- bs.x[,i]*bs.y[,j]
    }
  }
  bs.xy <- scale(bs.xy,scale=F)

  Y   <- assays(spe)$counts
  idx <- sample(ncol(Y),size=min(3000,round(pCells.touse*ncol(Y))))
  Ysub<- as.matrix(Y[,idx]) 
  nsub<- ncol(Ysub)
  W   <- model.matrix(~spe$logLS*bs.xy)[,-1]
  nW  <- ncol(W)
  Wsub<- W[idx,]

  # initial estimate
  lmu.hat<- Z <- log(Ysub+1) 
  gmean  <- rowMeans(Z)
  alpha  <- matrix(0,nrow(Ysub),ncol(W))
  # alpha for logLS
  alpha[,1] <- 1

  # estimate initial psi
  psi.idx <- sample(ncol(Ysub),size=min(50,round(0.1*ncol(Y))))
  bef <- Sys.time() 
  offs   <- Wa <- Matrix::tcrossprod(alpha,Wsub) 
  offs.psi <- offs[,psi.idx]
  Ysub.psi <- Ysub[,psi.idx]
  nsub.psi  <- ncol(Ysub.psi)
  psi <- edgeR::estimateDisp(Ysub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE)$tagwise.dispersion

  # initiate
  conv <- FALSE
  iter.outer <- 0
  logl.outer <- NULL
  step <- rep(1,nsub)
  message('Start...')
  while(!conv) {
    iter.outer <- iter.outer + 1
    logl.beta <- NULL
    
    ### calculate initial loglik
    lmu.hat    <- gmean + Wa
    sig.inv<- 1/(exp(-lmu.hat) +outer(psi,rep(1,nsub),'/'))

    temp <- dnbinom(Ysub,mu=exp(lmu.hat),size=1/psi,log=TRUE)
    loglik <- colSums(temp)

    #############################
    step <- rep(1,nsub)
    conv.beta <- FALSE
    # saving temporary best estimate
    best.W <- Wsub ; best.psi <- psi ; best.a <- alpha; best.gmean <- gmean
    
    iter <- 1
    halving <- 0
    while(!conv.beta) {
      # print(paste0('Inner iter:', iter))
      # calc working vector Z for a subset of cells 
      lmu.hat    <- gmean + Matrix::tcrossprod(alpha,Wsub) 

      sig.inv<- 1/(exp(-lmu.hat) +outer(psi,rep(1,nsub),'/'))
      Z      <- lmu.hat + sweep( ((Ysub+0.01)/(exp(lmu.hat)+0.01) - 1),2,step,'*')

      # store current alpha est
      alpha.old <- alpha

      # get alpha 
      bef <- Sys.time()
      wt.cell  <- matrixStats::colMeans2(sig.inv)
      # prevent outlier large weight
      wt.cell  <- pmin(wt.cell, quantile(wt.cell,probs=0.98))
      b         <- (Z-gmean) %*% diag(wt.cell) %*% Wsub
      alpha<- b %*% Matrix::chol2inv(chol(Matrix::crossprod(Wsub*wt.cell,Wsub)))
      
      #print(dim(alpha))
      #print(dim(b))
      #print(dim(Wsub))
      #print(length(wt.cell))

      alpha[,1] <- mean(alpha[,1])
      b         <- (Z-gmean-Matrix::tcrossprod(alpha[,1,drop=FALSE],Wsub[,1,drop=FALSE])) %*% diag(wt.cell) %*% Wsub[,-1]
      alpha[,-1]<- b %*% Matrix::chol2inv(chol(Matrix::crossprod(Wsub[,-1]*wt.cell,Wsub[,-1])+lambda.a*diag(nW-1)))

      aft <- Sys.time()

      # if new alpha est has missing/inf values, revert to previous estimate
      if(any(is.na(alpha) | is.infinite(alpha) )) alpha <- alpha.old

      # reduce outliers
      a.med <- matrixStats::colMedians(alpha)
      a.mad <- matrixStats::colMads(alpha)
      for(i in 1:nW) {
        alpha[which(alpha[,i]> (a.med[i]+4*a.mad[i])),i] <- a.med[i]+4*a.mad[i]
        alpha[which(alpha[,i]< (a.med[i]-4*a.mad[i])),i] <- a.med[i]-4*a.mad[i]
      } 

      #step2: update gmean
      Z.res  <- Z -  Matrix::tcrossprod(alpha,Wsub)  
      bef <- Sys.time()
      gmean<- matrixStats::rowSums2(Z.res*sig.inv)/matrixStats::rowSums2(sig.inv)
      aft <- Sys.time()

      # calculate current logl
      lmu.hat    <- gmean + Matrix::tcrossprod(alpha,Wsub) 
      temp <- dnbinom(Ysub,mu=exp(lmu.hat),size=1/psi,log=TRUE)
      loglik.tmp <- colSums(temp) 

      # check degenerate case
      if(iter>=1) {
        degener<- sum(loglik)>sum(loglik.tmp)
        degener[is.na(degener)] <- TRUE
        if(degener) {
          check.gene <- loglik>loglik.tmp
          step[check.gene] <- step[check.gene]*step.fac
          Wsub <- best.W; alpha <- best.a ; psi <- best.psi ; gmean <- best.gmean 
          halving <- halving + 1
          if(halving>=3) {
            loglik.tmp <- loglik
            degener <- !degener
          }
        }
        if(!degener) {
          loglik     <- loglik.tmp
          logl.beta  <- c(logl.beta,sum(loglik))
          #print(paste0('Outer Iter ',iter.outer, ', Inner iter ', iter, ' logl-likelihood:', logl.beta[length(logl.beta)]))
          # saving temporary best estimate
          best.W <- Wsub; best.psi <- psi ; best.a <- alpha; best.gmean <- gmean 
          iter <- iter + 1
          halving <- 0
        }
      }

      conv.logl <- FALSE
      if(iter>2) {
        conv.logl <- ifelse( (logl.beta[length(logl.beta)] - logl.beta[length(logl.beta)-1])/abs(logl.beta[length(logl.beta)-1]) < 1e-04,TRUE,FALSE)
      }
      conv.beta <- iter>=inner.maxit  | conv.logl
    } # end of IRLS inner loop

    #update psi
    logl.outer.tmp <- sum(loglik,na.rm=T)
    updt.psi <- TRUE
    if(iter.outer>1)
      updt.psi <- abs(logl.outer[length(logl.outer)]-logl.outer.tmp)/abs(logl.outer[length(logl.outer)]) > 1e-08  

    Wa <- offs <- Matrix::tcrossprod(best.a,best.W) 
    offs.psi <- offs[,psi.idx]
    if(updt.psi) 
      psi.new <- tryCatch({ edgeR::estimateDisp(Ysub.psi,as.matrix(rep(1,nsub.psi)),offset=offs.psi,tagwise=TRUE,robust=TRUE)$tagwise.dispersion},error = function(e) {NA})
    # update psi
    psi[!is.na(psi.new) & !is.infinite(psi.new)] <- psi.new[!is.na(psi.new) & !is.infinite(psi.new)]


    # record outer logl
    Wsub <- best.W
    alpha <- best.a
    gmean <- best.gmean
    # new changess
    best.psi <- psi

    # recalculate logl
    lmu.hat    <- gmean + Matrix::tcrossprod(alpha,Wsub) 
    temp       <- dnbinom(Ysub,mu=exp(lmu.hat),size=1/psi,log=TRUE)
    logl.outer <- c(logl.outer, sum(temp)) #-lambda.a*sum(scale(alpha,scale=FALSE)^2)-lambda.b*sum(Mb^2))
    message(paste0('Outer Iter ',iter.outer, ' logl-likelihood:', logl.outer[length(logl.outer)]))

    if(iter.outer>1) {
      conv <- ifelse( (logl.outer[iter.outer] - logl.outer[iter.outer-1])/abs(logl.outer[iter.outer-1]) < 1e-04 | iter.outer>=outer.maxit,TRUE,FALSE)
    }
  }
  # calculate residual
  mu      <- exp( gmean + Matrix::tcrossprod(alpha,W))
  # winsorize mean (to avoid slow process for extremely large mean)
  max.val  <- exp(matrixStats::rowMedians(log(mu)) + 4*matrixStats::rowMads(log(mu)))
  winsorize<- mu>max.val
  mu <- mu*(1-winsorize) + winsorize*max.val

  pred2.mu<- exp(gmean + Matrix::tcrossprod(alpha[,2:37,drop=FALSE],W[,2:37,drop=FALSE]))
  max.val  <- exp(matrixStats::rowMedians(log(pred2.mu)) + 4*matrixStats::rowMads(log(pred2.mu)))
  winsorize<- pred2.mu>max.val
  pred2.mu <- pred2.mu*(1-winsorize) + winsorize*max.val

  # bound dispersion parameters (large disp slow the process)
  max.val  <- exp( median(log(psi)) + 3*mad(log(psi)) )
  winsorize<- psi>max.val
  psi <- psi*(1-winsorize) + winsorize*max.val

  #logPAC
  lb <- pnbinom(as.matrix(Y)-1,mu=mu, size= 1/psi)
  ub <- dnbinom(as.matrix(Y),mu=mu, size=1/psi) + lb
  p  <- (lb + ub)/2
  p  <- pmax(pmin(p, 0.999), 0.001)
  # return logPAC
  PAC <- qnbinom(p,mu=scale.factor*pred2.mu,size=1/psi)
  assays(spe,withDimnames=FALSE)$logPAC <- log(PAC+1)
  # return pearson residuals
  assays(spe,withDimnames=FALSE)$pearson <- (Y- exp(Matrix::tcrossprod(alpha[,-2:-37,drop=FALSE],W[,-2:-37,drop=FALSE])))/sqrt(mu+mu^2 * psi )
  # return log mean 'BIO' expression
  assays(spe,withDimnames=FALSE)$medBIO  <- log(qnbinom(0.5,mu=scale.factor*pred2.mu,size=1/psi)+1)
  assays(spe,withDimnames=FALSE)$meanBIO <- log(pred2.mu)
  return(spe=spe)
}

