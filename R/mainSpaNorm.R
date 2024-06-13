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
#' @param batch a vector or numeric matrix, specifying the batch design to regress out (default NULL, representing no batch effects). See details for more information on how to define this variable.
#' @param step.factor a numeric, specifying the multiplicative factor to decrease IRLS step by when log-likelihood diverges (default is 0.5).
#' @param maxit.nb a numeric, specifying the maximum number of IRLS iteration for estimating NB mean parameters for a given dispersion parameter (default is 50).
#' @param maxit.psi a numeric, specifying the maximum number of IRLS iterations to estimate the dispersion parameter (default is 25).
#' @param tol a numeric, specifying the tolerance for convergence (default is 1e-4).
#' @param verbose a logical, specifying wether to show update messages (default TRUE).
#' @param ... other parameters fitting parameters.
#' 
#' @details SpaNorm works by first fitting a spatial regression model for library size to the data. Normalised data can then be computed using various adjustment approaches. When a negative binomial gene-model is used, the data can be adjusted using the following approaches: 'logpac', 'pearson', 'medbio', and 'meanbio'.
#' 
#' Batch effects can be specified using the `batch` parameter. If this parameter is a vector, a design matrix will be created within the function using `model.matrix`. If a custom design is provided in the form of a numeric matrix, this should ideally be created using `model.matrix`. The batch matrix should be created with an intercept term. The SpaNorm function will automatically detect the intercept term and remove the relevant column. Alternatively, users can subset the model matrix to remove this column manually. Please note that the model formula should include the intercept term and that the intercept column should be subset out after.
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
    df.tps = 6,
    lambda.a = 0.0001,
    batch = NULL,
    tol = 1e-4,
    step.factor = 0.5,
    maxit.nb = 50,
    maxit.psi = 25,
    verbose = TRUE,
    ...) {
  standardGeneric("SpaNorm")
})

#' @rdname SpaNorm
setMethod(
  "SpaNorm",
  signature("SpatialExperiment"),
  function(spe, sample.p, gene.model, adj.method, scale.factor, df.tps, lambda.a, batch, tol, step.factor, maxit.nb, maxit.psi, verbose, ...) {
    checkSPE(spe)
    adj.method = match.arg(adj.method)
    gene.model = match.arg(gene.model)
    
    # message function depending on verbose param
    msgfun = ifelse(verbose, message, \(...){})

    # extract counts, coords, and size factors
    emat = SummarizedExperiment::assay(spe, "counts")
    coords = SpatialExperiment::spatialCoords(spe)
    logLS = log(pmax(1e-08, SingleCellExperiment::sizeFactors(spe)))

    # fit/retrieve SpaNorm model
    precomputed = validSpaNormSPE(spe, gene.model)
    fit = S4Vectors::metadata(spe)$SpaNorm
    if (precomputed &&
        fit$sample.p == sample.p &&
        fit$gene.model == gene.model &&
        fit$df.tps == df.tps &&
        all.equal(fit$batch, batch)
      ) {
      msgfun("(1/2) Retrieve precomputed SpaNorm model")
      fit.spanorm = S4Vectors::metadata(spe)$SpaNorm
    } else{
      msgfun("(1/2) Fitting SpaNorm model")
      fit.spanorm = fitSpaNorm(Y = emat, coords = coords, sample.p = sample.p, gene.model = gene.model, msgfun = msgfun, df.tps = df.tps, lambda.a = lambda.a, batch = batch, logLS = logLS, tol = tol, step.factor = step.factor, maxit.nb = maxit.nb, maxit.psi = maxit.psi)
      # add model to assay
      S4Vectors::metadata(spe)$SpaNorm = fit.spanorm
    }

    # computed normalised data
    msgfun("(2/2) Normalising data")
    adj.fun = getAdjustmentFun(gene.model, adj.method)
    normmat = adj.fun(emat, scale.factor, fit.spanorm)
    if ("logcounts" %in% SummarizedExperiment::assayNames(spe)) {
      warning("'logcounts' exists and will be overwritten")
    }
    SummarizedExperiment::assay(spe, "logcounts") = normmat

    return(spe)
  }
)

sampleRandom <- function(coords, nsub) {
  idx = rep(FALSE, nrow(coords))
  idx[sample.int(nrow(coords), size = nsub)] = TRUE
  return(idx)
}

fitSpaNorm <- function(Y, coords, sample.p, gene.model, df.tps = 6, lambda.a = 0.0001, batch, logLS, msgfun = message, ...) {
  # parameter checks
  if (sample.p <= 0 | sample.p > 1) {
    stop("'sample.p' should be in the interval (0,1]")
  }
  if (lambda.a <= 0) {
    stop("'lambda.a' should be greater than 0")
  }
  lambda.a = lambda.a * ncol(Y)

  # prepare splines
  bs.xy = bs.tps(coords[, 1], coords[, 2], df.tps = df.tps) # get basis for the thin-plate spline

  # calculate effective library size if not precomputed
  if (is.null(logLS)) {
    cl = scran::quickCluster(Y)
    logLS = log(pmax(1e-08, scran::calculateSumFactors(Y, clusters = cl)))
  }

  # setting-up variables for the NB regression models
  W = model.matrix(~ logLS * bs.xy)[, -1]
  # add batch design matrix
  W = cbind(W, checkBatch(batch, ncol(Y)))

  # sample data for computational efficiency
  maxn = 3000
  nsub = round(sample.p * ncol(Y))
  msgfun(sprintf("%d cells/spots sampled to fit model", nsub))
  if (nsub > maxn) {
    warning(sprintf("consider reducing 'sample.p' to %.2f to increase computational efficiency", max(floor(maxn / ncol(Y) * 100) / 100, 0.01)))
  } else if (nsub == 0) {
    stop(sprintf("'sample.p' is too small, consider using %.2f", min(1, floor(maxn / ncol(Y) * 100) / 100)))
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
  fit.spanorm$isbio[seq(2, df.tps^2 + 1)] = TRUE
  fit.spanorm$batch = batch

  return(fit.spanorm)
}

validSpaNormSPE <- function(spe, gene.model) {
  valid = FALSE

  # check validity based on model
  if (gene.model == "nb") {
    valid = validSpaNormNBSPE(spe)
  }

  return(valid)
}

validSpaNormNBSPE <- function(spe) {
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

checkBatch <- function(batch, nobs) {
  # if null, do nothing
  if (is.null(batch)) {
    batch = c()
  }
  
  # if matrix, check and return
  if (is.matrix(batch)) {
    # check dimensions
    if (nrow(batch) != nobs) {
      stop("number of rows in the 'batch' matrix do not match number of cells/spots")
    }

    # check type
    if (!is.numeric(batch)) {
      stop("'batch' should be a numeric matrix (consider using 'model.matrix' to define the design)")
    }

    # check for intercept
    isintercept = grepl("intercept", colnames(batch), ignore.case = TRUE)
    isintercept = isintercept | matrixStats::colAlls(batch == 1)
    if (any(isintercept)) {
      warning("'intercept' term detected and will be removed")
      batch = batch[, !isintercept, drop = FALSE]
    }
  }
  
  # if vector
  if (is.vector(batch)) {
    # check dimensions
    if (length(batch) != nobs) {
      stop("length of 'batch' vector does not match number of cells/spots")
    }
    batch = model.matrix(~batch)[, -1, drop = FALSE]
  }

  return(batch)
}
