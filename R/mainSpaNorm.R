#' Spatially-dependent normalisation for spatial transcriptomics data
#'
#' Performs normalisation of spatial transcriptomics data using spatially-dependent spot- and gene- specific size factors.
#'
#' @param spe a SpatialExperiment or Seurat object, with the count data stored in 'counts' or 'data' assays respectively.
#' @param sample.p a numeric, specifying the (maximum) proportion of cells/spots to sample for model fitting (default is 0.25).
#' @param gene.model a character, specifying the model to use for gene/protein abundances (default 'nb'). This should be 'nb' for count based datasets.
#' @param adj.method a character, specifying the method to use to adjust the data (default 'auto', see details)
#' @param scale.factor a numeric, specifying the sample-specific scaling factor to scale the adjusted count.
#' @param df.tps a numeric, of length 1 or 2, specifying the maximum degrees of freedom for the thin-plate spline for the biology and library size effects, respectively (default is 6, see details).
#' @param lambda.a a numeric, of length 1 or 2, specifying the smoothing parameter for regularizing regression coefficients (default is 0.0001, see details). Actual lambda.a used is lambda.a * ncol(spe).
#' @param batch a vector or numeric matrix, specifying the batch design to regress out (default NULL, representing no batch effects). See details for more information on how to define this variable.
#' @param step.factor a numeric, specifying the multiplicative factor to decrease IRLS step by when log-likelihood diverges (default is 0.5).
#' @param maxit.nb a numeric, specifying the maximum number of IRLS iteration for estimating NB mean parameters for a given dispersion parameter (default is 50).
#' @param maxit.psi a numeric, specifying the maximum number of IRLS iterations to estimate the dispersion parameter (default is 25).
#' @param maxn.psi a numeric, specifying the maximum number of cells/spots to sample for dispersion estimation (default is 500).
#' @param tol a numeric, specifying the tolerance for convergence (default is 1e-4).
#' @param overwrite a logical, specifying whether to force recomputation and overwrite an existing fit (default FALSE). Note that if df.tps, batch, lambda.a, or gene.model are changed, the model is recomputed and overwritten.
#' @param verbose a logical, specifying whether to show update messages (default TRUE).
#' @param ... other parameters fitting parameters.
#' 
#' @details SpaNorm works by first fitting a spatial regression model for library size to the data. Normalised data can then be computed using various adjustment approaches. When a negative binomial gene-model is used, the data can be adjusted using the following approaches: 'logpac', 'pearson', 'medbio', and 'meanbio'.
#' 
#' The `df.tps` parameter specifies the degrees of freedom for the thin-plate spline. If only 1 value is provided, it specifies the degrees of freedom of the biology with the degrees of freedom of the library size being half of that. If 2 values are provided, the first value specifies the degrees of freedom of the biology and the second value specifies the degrees of freedom of the library size. For rectangular tissues, df.tps specifies the degrees of freedom along the length, with the degrees of freedom along the width calculated ceiling(width / length * df.tps).
#' 
#' Similarly, the `lambda.a` parameter specifies the smoothing parameter for regularizing regression coefficients. If only 1 value is provided, it specifies the lambda.a for both the biology and library size functions. If 2 values are provided, the first value specifies the lambda.a for the biology and the second value specifies the lambda.a for the library size. Batch effects are not regularised.
#' 
#' Batch effects can be specified using the `batch` parameter. If this parameter is a vector, a design matrix will be created within the function using `model.matrix`. If a custom design is provided in the form of a numeric matrix, this should ideally be created using `model.matrix`. The batch matrix should be created with an intercept term. The SpaNorm function will automatically detect the intercept term and remove the relevant column. Alternatively, users can subset the model matrix to remove this column manually. Please note that the model formula should include the intercept term and that the intercept column should be subset out after.
#' 
#' @return a SpatialExperiment or Seurat object with the adjusted data stored in 'logcounts' or 'data', respectively.
#' @name SpaNorm
#' 
#' @examples 
#' data(HumanDLPFC)
#' \donttest{
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
    maxn.psi = 500,
    overwrite = FALSE,
    verbose = TRUE,
    ...) {
  standardGeneric("SpaNorm")
})

#' @rdname SpaNorm
setMethod(
  "SpaNorm",
  signature("SpatialExperiment"),
  function(spe, sample.p, gene.model, adj.method, scale.factor, df.tps, lambda.a, batch, tol, step.factor, maxit.nb, maxit.psi, overwrite, verbose, ...) {
    checkSPE(spe)
    adj.method = match.arg(adj.method)
    gene.model = match.arg(gene.model)
    df.tps = as.integer(df.tps)
    
    # message function depending on verbose param
    msgfun = ifelse(verbose, message, \(...){})

    # extract counts, coords, and size factors
    emat = SummarizedExperiment::assay(spe, "counts")
    coords = SpatialExperiment::spatialCoords(spe)
    LS = SingleCellExperiment::sizeFactors(spe)

    # fit/retrieve SpaNorm model
    fit.spanorm = getSpaNormFit(spe, validate = FALSE)
    if (!overwrite &&
        !is.null(fit.spanorm) &&
        fit.spanorm$ngenes == nrow(spe) &&
        fit.spanorm$ncells == ncol(spe) &&
        fit.spanorm$gene.model == gene.model &&
        matchDftps(fit.spanorm$df.tps, df.tps) &&
        matchLambda(fit.spanorm$lambda.a, lambda.a) &&
        all.equal(fit.spanorm$batch, batch)
      ) {
      msgfun("(1/2) Retrieve precomputed SpaNorm model")
    } else{
      msgfun("(1/2) Fitting SpaNorm model")
      fit.spanorm = fitSpaNorm(Y = emat, coords = coords, sample.p = sample.p, gene.model = gene.model, msgfun = msgfun, df.tps = df.tps, lambda.a = lambda.a, batch = batch, LS = LS, tol = tol, step.factor = step.factor, maxit.nb = maxit.nb, maxit.psi = maxit.psi, maxn.psi = maxn.psi)
      # add model to assay
      S4Vectors::metadata(spe)$SpaNorm = fit.spanorm
    }

    # computed normalised data
    msgfun("(2/2) Normalising data")
    if (!any(fit.spanorm$wtype == "biology")) {
      stop("'SpaNorm' fit should have at least one column representing 'biology'")
    }
    adj.fun = getAdjustmentFun(gene.model, adj.method)
    normmat = adj.fun(emat, scale.factor, fit.spanorm)
    if ("logcounts" %in% SummarizedExperiment::assayNames(spe)) {
      warning("'logcounts' exists and will be overwritten")
    }
    SummarizedExperiment::assay(spe, "logcounts") = methods::as(normmat, "sparseMatrix")

    return(spe)
  }
)

sampleRandom <- function(coords, nsub) {
  idx = rep(FALSE, nrow(coords))
  idx[sample.int(nrow(coords), size = nsub)] = TRUE
  return(idx)
}

fitSpaNorm <- function(Y, coords, sample.p, gene.model, df.tps = 6, lambda.a = 0.0001, batch, LS, msgfun = message, ...) {
  # parameter checks
  if (!gene.model %in% getGeneModels()) {
    stop(sprintf("'gene.model' should be one of: %s", paste(getGeneModels(), collapse = ", ")))
  }
  if (sample.p <= 0 | sample.p > 1) {
    stop("'sample.p' should be in the interval (0,1]")
  }
  if (any(lambda.a < 0)) {
    stop("'lambda.a' should be positive")
  }
  if (!length(lambda.a) %in% c(1, 2)) {
    stop("'lambda.a' should be a single value or a vector of length 2")
  }
  if (length(lambda.a) == 1) {
    lambda.a = rep(lambda.a, 2)
  }

  # scale coordinates
  coords = apply(coords, 2, \(x) {
    (x - min(x)) / (max(x) - min(x)) - 0.5
  })

  # get basis for the thin-plate spline
  if (length(df.tps) == 1) {
    df.tps.bio = df.tps
    df.tps.ls = ceiling(df.tps.bio / 2)
  } else if (length(df.tps) == 2) {
    df.tps.bio = df.tps[1]
    df.tps.ls = df.tps[2]
  } else {
    stop("'df.tps' should be a single integer or a vector of length 2")
  }
  bs.xy.bio = bs.tps(coords[, 1], coords[, 2], df.tps = df.tps.bio) # biology
  bs.xy.ls = bs.tps(coords[, 1], coords[, 2], df.tps = df.tps.ls) # library size
  df.tps = c(as.integer(attr(bs.xy.bio, "df.tps")), as.integer(attr(bs.xy.ls, "df.tps")))

  # calculate effective library size if not precomputed
  if (is.null(LS)) {
    cl = scran::quickCluster(Y)
    logLS = log(pmax(1e-08, scran::calculateSumFactors(Y, clusters = cl)))
  } else {
    logLS = log(pmax(1e-08, LS))
  }

  # setting-up variables for the NB regression models
  W = model.matrix(~ logLS + bs.xy.bio + logLS:bs.xy.ls)[, -1]
  # add batch design matrix
  W = cbind(W, checkBatch(batch, ncol(Y)))

  # mark factors representing biology of interest
  wtype = rep("batch", ncol(W))
  wtype[seq(2, prod(df.tps[1:2]) + 1)] = "biology"
  wtype[c(1, seq(prod(df.tps[1:2]) + 2, prod(df.tps[1:2]) + prod(df.tps[3:4]) + 1))] = "ls"

  # create lambda.a vector
  lambda.a.vec = rep(0, length(wtype))
  lambda.a.vec[wtype == "biology"] = lambda.a[1]
  lambda.a.vec[wtype == "ls"] = lambda.a[2]
  lambda.a.vec = lambda.a.vec[-1]

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
    fit.spanorm.nb = fitSpaNormNB(Y, W, idx, ..., lambda.a = lambda.a.vec * ncol(Y), msgfun = msgfun)
  }

  # create object
  fit.spanorm = SpaNormFit(
    ngenes = nrow(Y),
    ncells = ncol(Y),
    gene.model = gene.model,
    df.tps = df.tps,
    sample.p = sample.p,
    lambda.a = lambda.a,
    batch = batch,
    W = W,
    alpha = fit.spanorm.nb$alpha,
    gmean = fit.spanorm.nb$gmean,
    psi = fit.spanorm.nb$psi,
    wtype = wtype,
    loglik = fit.spanorm.nb$loglik,
    sampling = fit.spanorm.nb$sampling
  )

  return(fit.spanorm)
}

matchDftps <- function(df.fit, df.par) {
  stopifnot(length(df.par) %in% c(1, 2))
  stopifnot(length(df.fit) == 4)
  
  if (length(df.par) == 2) {
    if (df.par[1] == max(df.fit[1:2]) && df.par[2] == max(df.fit[3:4])) {
      return(TRUE)
    }
  } else if (length(df.par) == 1) {
    if (df.par == max(df.fit[1:2]) && ceiling(df.par / 2) == max(df.fit[3:4])) {
      return(TRUE)
    }
  }

  return(FALSE)
}

matchLambda <- function(lambda.fit, lambda.par) {
  stopifnot(length(lambda.par) %in% c(1, 2))
  stopifnot(length(lambda.fit) == 2)
  
  if (all(lambda.fit == lambda.par)) {
    return(TRUE)
  }

  return(FALSE)
}

bs.tps <- function(x, y, df.tps = 6) {
  stopifnot(length(df.tps) == 1)
  stopifnot(length(x) == length(y))

  # checks
  if (df.tps <= 0) {
    stop("'df.tps' should be greater than 0")
  }
  if (df.tps - as.integer(df.tps) != 0) {
    stop("'df.tps' should be an integer")
  }

  # determine df along each axis
  xrng = diff(range(x))
  yrng = diff(range(y))
  gap = max(xrng, yrng) / df.tps
  df.tps.x = ceiling(xrng / gap)
  df.tps.y = ceiling(yrng / gap)

  # construct spline
  bs.x = splines::ns(x, df = df.tps.x)
  bs.y = splines::ns(y, df = df.tps.y)
  bs.xy = matrix(0, nrow = length(x), ncol = df.tps.x * df.tps.y)
  for (i in seq_len(df.tps.x)) {
    for (j in seq_len(df.tps.y)) {
      bs.xy[, (i - 1) * ncol(bs.x) + j] <- bs.x[, i] * bs.y[, j]
    }
  }
  bs.xy = scale(bs.xy, scale = FALSE)

  # add spline dimensions to attributes
  attr(bs.xy, 'df.tps') = c(df.tps.x, df.tps.y)

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

getGeneModels <- function() {
  c("nb")
}

getSpaNormFit <- function(spe, null = FALSE, validate = TRUE) {
  name = ifelse(null, "SpaNormNull", "SpaNorm")

  # Retrieve model
  if (is(spe, "SpatialExperiment")) {
    fit.spanorm <- S4Vectors::metadata(spe)[[name]]
  } else {
    stop("'spe' must be a SpatialExperiment")
  }

  # Basic validation
  if (validate) {
    # Check if model exists
    if (is.null(fit.spanorm)) {
      stop(sprintf("%s model not found in 'spe'. Please run '%s' on the `spe` object first.", name, ifelse(null, "SpaNormSVG", "SpaNorm")))
    }

    # Check dimensions match
    if ((fit.spanorm$ngenes != nrow(spe) ||
    fit.spanorm$ncells != ncol(spe))) {
      stop(sprintf("%s model dimensions do not match the data. Please rerun '%s'.", name, ifelse(null, "SpaNormSVG", "SpaNorm")))
    }
  }

  return(fit.spanorm)
}
