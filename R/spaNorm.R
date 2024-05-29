#' Spatially-dependent normalisation for spatial transcriptomics datas
#'
#' Performs normalisation of spatial transcriptomics data using spatially-dependent spot- and gene- specific size factors.
#'
#' @param spe a SpatialExperiment or Seurat object, with the count data stored in 'counts' or 'data' assays respectively.
#' @param sample.p a numeric, specifying the (maximum) proportion of cells/spots to sample for model fitting (default is 0.25).
#' @param gene.model a character, specifying the model to use for gene/protein abundances (default 'nb'). This should be 'nb' for count based datasets.
#' @param adj.method a character, specifying the method to use to adjust the data (default 'auto', see details)
#' @param scale.factor a numeric, specifying the sample-specific scaling factor to scale the adjusted count.
#' @param df.tps a numeric, specifying the degrees of freedom for the thin-plate spline (default is 6).
#' @param lambda.a a numeric, specifying the smoothing parameter for regularizing regression coefficients (default is 0.0001). Actual lambda.a used is lambda.a * ncol(spe).
#' @param step.factor a numeric, specifying the multiplicative factor to decrease IRLS step by when log-likelihood diverges (default is 0.5).
#' @param maxit.nb a numeric, specifying the maximum number of IRLS iteration for estimating NB mean parameters for a given dispersion parameter (default is 50).
#' @param maxit.psi a numeric, specifying the maximum number of IRLS iterations to estimate the dispersion parameter (default is 25).
#' @param tol a numeric, specifying the tolerance for convergence (default is 1e-4).
#' @param verbose a logical, specifying wether to show update messages (default TRUE).
#' @param ... other parameters to pass to SpaNorm.
#' 
#' @details SpaNorm works by first fitting a spatial regression model for library size to the data. Normalised data can then be computed using various adjustment approaches. When a negative binomial gene-model is used, the data can be adjusted using the following approaches: 'logpac', 'pearson', 'medbio', and 'meanbio'.
#' 
#' @return a SpatialExperiment or Seurat object with the adjusted data stored in 'logcounts' or 'data', respectively.
#' 
#' @examples 
#' data(HumanDLPFC)
#' HumanDLPFC = HumanDLPFC[filterGenes(HumanDLPFC, 0.5), ]
#' SpaNorm(HumanDLPFC, df.tps = 2, tol = 1e-2)
#' @export
#' 
setGeneric("SpaNorm", function(
    spe,
    sample.p = 0.25,
    gene.model = c("nb"),
    adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
    scale.factor = 1,
    verbose = TRUE,
    ...) {
  standardGeneric("SpaNorm")
})

#' @rdname filterGenes
setMethod(
  "SpaNorm",
  signature("SpatialExperiment"),
  function(spe, sample.p, gene.model, adj.method, scale.factor, verbose = TRUE, ...) {
    checkSPE(spe)
    adj.method = match.arg(adj.method)
    gene.model = match.arg(gene.model)
    
    # message function depending on verbose param
    msgfun = ifelse(verbose, message, \(...){})

    # extract counts and coords
    emat = SummarizedExperiment::assay(spe, "counts")
    coords = SpatialExperiment::spatialCoords(spe)

    # fit SpaNorm model
    msgfun("(1/2) Fitting SpaNorm model")
    fit.spanorm = fitSpaNorm(emat, coords, sample.p, gene.model, msgfun = msgfun, ...)

    # computed normalised data
    msgfun("(2/2) Normalising data")
    adj.fun = getAdjustmentFun(gene.model, adj.method)
    normmat = adj.fun(emat, scale.factor, fit.spanorm)
    
    # add model to assay
    SingleCellExperiment::logcounts(spe) = normmat
    SummarizedExperiment::metadata(spe) = fit.spanorm

    return(spe)
  }
)

getAdjustmentFun <- function(gene.model, adj.method) {
  if (gene.model == "nb") {
    adj.fun = switch(
      adj.method,
      auto = normaliseLogPAC,
      logpac = normaliseLogPAC,
      pearson = normalisePearson,
      medbio = normaliseMedianBio,
      meanbio = normaliseMeanBio,
      stop("invalid argument for 'adj.method'")
    )
  } else {
    stop(sprintf("'%s' gene model not supported", gene.model))
  }
  return(adj.fun)
}

sampleRandom <- function(coords, nsub) {
  idx = rep(FALSE, nrow(coords))
  idx[sample.int(nrow(coords), size = nsub)] = TRUE
  return(idx)
}

fitSpaNorm <- function(Y, coords, sample.p, gene.model, df.tps = 6, lambda.a = 0.0001, msgfun = message, ...) {
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
  bs.xy = bs.tps(coords[, 1], coords[, 2], df.tps = df.tps) # get basis for the thin-plate spline
  W = model.matrix(~ logLS * bs.xy)[, -1]

  # sample data for computational efficiency
  nsub = round(sample.p * ncol(Y))
  msgfun(sprintf("%d cells/spots sampled to fit model", nsub))
  if (nsub > 3000) {
    warning(sprintf("consider reducing 'sample.p' to %.2f to increase computational efficiency", floor(3000 / ncol(Y) * 100) / 100))
  } else if (nsub == 0) {
    stop(sprintf("'sample.p' is too small, consider using %.2f", min(1, floor(3000 / ncol(Y) * 100) / 100)))
  }
  # random sampling
  idx = sampleRandom(coords, nsub)

  # fit model
  if (gene.model == "nb") {
    fit.spanorm = fitSpaNormNB(Y, W, idx, ..., lambda.a = lambda.a, msgfun = msgfun)
  }

  # add W
  fit.spanorm$W = W
  # mark factors representing biology of interest
  fit.spanorm$isbio = rep(FALSE, ncol(W))
  fit.spanorm$isbio[seq(2, df.tps ^ 2 + 1)] = TRUE

  return(fit.spanorm)
}

fitSpaNormNB <- function(Y, W, idx, maxit.psi = 25, tol = 1e-4, ..., msgfun = message) {
  # parameter checks
  if (maxit.psi <= 0) {
    stop("'maxit.psi' should be greater than 0")
  }
  if (maxit.psi - as.integer(maxit.psi) != 0) {
    stop("'maxit.psi' should be integers")
  }

  # message function for inner loop
  msgfunNB = \(x, ...) {\(...) msgfun(x, ...)}

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
  psi.idx = rep(FALSE, nsub)
  psi.idx[sample.int(nsub, size = nsub.psi)] = TRUE
  Ysub.psi = Ysub[, psi.idx, drop = FALSE]
  psi = rep(0, nrow(Ysub))
  
  # convergence trackers
  conv = FALSE
  logl.psi = c()
  iter = 1

  # outer loop of IRLS for estimating disp parameter
  while (!conv) {
    msgfun(sprintf("iter:%3d, estimating gene-wise dispersion", iter))
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
    loglik = colSums(dnbinom(Ysub, mu = exp(gmean + Wa), size = 1/psi, log = TRUE))
    msgfun(sprintf("iter:%3d, log-likelihood: %f", iter, sum(loglik)))

    # fit NB given dispersion estimates (psi) and extract required components
    msgfun(sprintf("iter:%3d, fitting NB model", iter))
    fit.nb = fitNBGivenPsi(Ysub, Wsub, gmean, alpha, psi, ..., loglik = loglik, msgfun = msgfunNB(sprintf("iter:%3d, ", iter)))
    gmean = fit.nb$gmean
    alpha = fit.nb$alpha
    loglik = fit.nb$loglik

    # check convergence
    logl.psi = c(sum(loglik), logl.psi) # push to stack
    conv.logl = FALSE
    iter = iter + 1 # similar post-increment to fitNBGivenPsi
    if (iter > 2) { # should be 2 because iter is post-incremented
      conv.logl = ifelse((logl.psi[1] - logl.psi[2]) / abs(logl.psi[2]) < tol, TRUE, FALSE)
    }
    conv = conv.logl | iter > maxit.psi
  }

  # convergence message
  if (iter > maxit.psi) {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (did not converged)", iter, logl.psi[1]))
  } else {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (converged)", iter, logl.psi[1]))
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
}

fitNBGivenPsi <- function(Ysub, Wsub, gmean, alpha, psi, lambda.a, step.factor = 0.5, maxit.nb = 50, tol = 1e-4, loglik = NULL, msgfun = message) {
  stopifnot(ncol(Ysub) == nrow(Wsub))
  stopifnot(nrow(Ysub) == length(gmean))
  stopifnot(nrow(Ysub) == nrow(alpha))
  stopifnot(nrow(Ysub) == length(psi))
  stopifnot(ncol(Wsub) == ncol(alpha))

  # parameter checks
  if (lambda.a <= 0) {
    stop("'lambda.a' should be greater than 0")
  }
  if (step.factor <= 0 | step.factor >= 1) {
    stop("'step.factor' should be in the interval (0,1)")
  }
  if (maxit.nb <= 0) {
    stop("'maxit.nb' should be greater than 0")
  }
  if (maxit.nb - as.integer(maxit.nb) != 0) {
    stop("'maxit.nb' should be integers")
  }

  # set the IRLS step size to 1
  nsub = ncol(Ysub)
  nW = ncol(Wsub)
  step = rep(1, nsub)
  
  if (is.null(loglik))  {
    # if not pre-computed (e.g., by outer loop), compute
    # helpful pattern for external use of the fitting function
    lmu.hat = gmean + Matrix::tcrossprod(alpha, Wsub)
    sig.inv = 1 / (exp(-lmu.hat) + outer(psi, rep(1, nsub), "/"))
    loglik = colSums(dnbinom(Ysub, mu = exp(lmu.hat), size = 1 / psi, log = TRUE))
  }

  # convergence trackers
  conv = FALSE
  logl.beta = c(sum(loglik)) # stack of loglik values - push values
  iter = 1

  # initialise halving counters
  halving = 0

  while (!conv) {
    # message
    msgfun(sprintf("iter:%3d, log-likelihood: %f", iter, logl.beta[1]))

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
      conv.logl = ifelse((logl.beta[1] - logl.beta[2]) / abs(logl.beta[2]) < tol, TRUE, FALSE)
    }
    conv = conv.logl | iter > maxit.nb
  }

  # convergence message
  if (iter > maxit.nb) {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (did not converged)", iter, logl.beta[1]))
  } else {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (converged)", iter, logl.beta[1]))
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

normaliseLogPAC <- function(Y, scale.factor, fit.spanorm) {
  # calculate mu with all effects
  mu = calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  # calculate mu with biology only
  isbio = fit.spanorm$isbio
  mu.2 = calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio], fit.spanorm$W[, isbio])
  # winsorize dispersion parameters (large disp slow the process)
  psi = fit.spanorm$psi
  psi.max = exp(median(log(psi)) + 3 * mad(log(psi)))
  psi = pmin(psi, psi.max)

  #logPAC
  lb = pnbinom(as.matrix(Y) - 1, mu = mu, size = 1/psi)
  ub = dnbinom(as.matrix(Y), mu = mu, size = 1/psi) + lb
  p = (lb + ub)/2
  p = pmax(pmin(p, 0.999), 0.001)

  # return logPAC
  normmat = log(qnbinom(p, mu = scale.factor * mu.2, size = 1/psi) + 1)
  colnames(normmat) = colnames(Y)
  rownames(normmat) = rownames(Y)

  return(normmat)
}

normaliseMeanBio <- function(Y, scale.factor, fit.spanorm) {
  isbio = fit.spanorm$isbio
  # calculate mu without library size effect
  normmat = log(calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio], fit.spanorm$W[, isbio])) # log(mu.2)
  colnames(normmat) = colnames(Y)
  rownames(normmat) = rownames(Y)

  return(normmat)
}

normaliseMedianBio <- function(Y, scale.factor, fit.spanorm) {
  # calculate mu without library size effect
  isbio = fit.spanorm$isbio
  mu.2 = calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio], fit.spanorm$W[, isbio])
  # winsorize dispersion parameters (large disp slow the process)
  psi = fit.spanorm$psi
  psi.max = exp(median(log(psi)) + 3 * mad(log(psi)))
  psi = pmin(psi, psi.max)

  # median Bio
  normmat = log(qnbinom(0.5, mu = scale.factor * mu.2, size = 1/psi) + 1)
  colnames(normmat) = colnames(Y)
  rownames(normmat) = rownames(Y)

  return(normmat)
}

normalisePearson <- function(Y, scale.factor, fit.spanorm) {
  isbio = fit.spanorm$isbio
  # calculate mu
  mu = calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  # winsorize dispersion parameters (large disp slow the process)
  psi = fit.spanorm$psi
  psi.max = exp(median(log(psi)) + 3 * mad(log(psi)))
  psi = pmin(psi, psi.max)

  # Pearson
  normmat = exp(Matrix::tcrossprod(fit.spanorm$alpha[, !isbio, drop=FALSE], fit.spanorm$W[, !isbio, drop = FALSE]))
  normmat = (Y - normmat) / sqrt(mu + mu ^ 2 * psi)
  colnames(normmat) = colnames(Y)
  rownames(normmat) = rownames(Y)

  return(normmat)
}

calculateMu <- function(gmean, alpha, W) {
  # calculate mu
  mu = gmean + Matrix::tcrossprod(alpha, W) # log(mu)
  # winsorise
  lmu.max = matrixStats::rowMedians(mu) + 4 * matrixStats::rowMads(mu)
  mu = exp(pmin(mu, lmu.max))

  return(mu)
}
