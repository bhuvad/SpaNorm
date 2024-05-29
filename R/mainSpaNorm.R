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
#' @name SpaNorm
#' 
#' @examples 
#' data(HumanDLPFC)
#' \dontrun{
#' SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#' }
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

#' @rdname SpaNorm
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

    # fit/retrieve SpaNorm model
    precomputed = validSpaNormSPE(spe)
    if (precomputed) {
      msgfun("(1/2) Retrieve precomputed SpaNorm model")
      fit.spanorm = S4Vectors::metadata(spe)$SpaNorm
    } else{
      msgfun("(1/2) Fitting SpaNorm model")
      fit.spanorm = fitSpaNorm(emat, coords, sample.p, gene.model, msgfun = msgfun, ...)
      # add model to assay
      S4Vectors::metadata(spe)$SpaNorm = fit.spanorm
    }

    # computed normalised data
    msgfun("(2/2) Normalising data")
    adj.fun = getAdjustmentFun(gene.model, adj.method)
    normmat = adj.fun(emat, scale.factor, fit.spanorm)
    SummarizedExperiment::assay(spe, "logcounts") = normmat

    return(spe)
  }
)

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

  # add params
  fit.spanorm$W = W
  fit.spanorm$df.tps = df.tps
  fit.spanorm$sample.p = sample.p
  fit.spanorm$gene.model = gene.model
  fit.spanorm$lambda.a = lambda.a
  # mark factors representing biology of interest
  fit.spanorm$isbio = rep(FALSE, ncol(W))
  fit.spanorm$isbio[seq(2, df.tps ^ 2 + 1)] = TRUE

  return(fit.spanorm)
}

validSpaNormSPE <- function(spe) {
  valid = FALSE

  # is object present in metadata
  if ("SpaNorm" %in% names(S4Vectors::metadata(spe))) {
    # is object valid
    fit.spanorm = S4Vectors::metadata(spe)$SpaNorm
    valid = checkNBParams(
      nrow(spe),
      ncol(spe),
      fit.spanorm$W,
      fit.spanorm$gmean,
      fit.spanorm$alpha,
      fit.spanorm$psi
    )
    # check isbio vector
    valid = valid & "isbio" %in% names(fit.spanorm) & length(fit.spanorm$isbio) == ncol(fit.spanorm$W)
    # check df.tps is valis
    valid = valid & "df.tps" %in% names(fit.spanorm) & ncol(fit.spanorm$W) >= fit.spanorm$df.tps^2 * 2 + 1 & sum(fit.spanorm$isbio) >= fit.spanorm$df.tps

    if (!valid) {
      warning("an invalid SpaNorm fit exists and will be replaced")
    }
  }

  return(valid)
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
