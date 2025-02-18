#' @importFrom stats dnbinom qnbinom pnbinom mad median quantile model.matrix pf p.adjust
NULL

#' Model-based spatially variable gene (SVG) calling
#'
#' Spatially variable gene (SVG) calling using the SpaNorm model.
#'
#' @param spe a SpatialExperiment or Seurat object, with the count data stored in 'counts' or 'data' assays respectively, and a SpaNorm model fit.
#' @param verbose a logical, specifying whether to show update messages (default TRUE).
#'
#' @details SpaNorm SVG calling works by using the SpaNorm model fit for data normalisation to perform a likelihood ratio test (LRT). The model used for normalisation is considered to be the full model. A second nested model is fit without the splines representing biology. These nested models are then compared using a LRT to identify genes where the splines representing biology contain strong signal.
#'
#' @return a SpatialExperiment or Seurat object with F-statistics, false discovery rates (FDRs). For SpatialExperiment objects, these are stored in the rowData.
#' @name SpaNormSVG
#'
#' @examples
#' 
#' library(SpatialExperiment)
#' library(ggplot2)
#' 
#' data(HumanDLPFC)
#' 
#' HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#' HumanDLPFC = SpaNormSVG(HumanDLPFC)
#' head(rowData(HumanDLPFC))
#' 
#' @export
#'
setGeneric("SpaNormSVG", function(
    spe,
    verbose = TRUE) {
  standardGeneric("SpaNormSVG")
})

#' @rdname SpaNormSVG
setMethod(
  "SpaNormSVG",
  signature("SpatialExperiment"),
  function(spe, verbose) {
    checkSPE(spe)

    # message function depending on verbose param
    msgfun = ifelse(verbose, message, \(...){})
    # Add progress tracking
    total_steps = 3
    current_step = 0
    
    report_progress <- function(msg) {
      if (verbose) {
        current_step <<- current_step + 1
        msgfun(sprintf("(%d/%d) %s", current_step, total_steps, msg))
      }
    }

    # Extract counts
    emat = SummarizedExperiment::assay(spe, "counts")
    
    # Check previous results
    svg.cols = c("svg.F", "svg.p", "svg.fdr")
    if (any(svg.cols %in% colnames(SummarizedExperiment::rowData(spe)))) {
      warning("SVG results exist in 'spe' and will be overwritten")
      cols = setdiff(colnames(SummarizedExperiment::rowData(spe)), svg.cols)
      SummarizedExperiment::rowData(spe) = SummarizedExperiment::rowData(spe)[, cols]
    }

    # Retrieve and validate SpaNorm model
    report_progress("Retrieving SpaNorm model")
    fit.spanorm = getSpaNormFit(spe)

    # Fit nested model
    report_progress("Fitting nested SpaNorm model") 
    fit.technical = fitSpaNormTechnical(emat, fit.spanorm, msgfun)

    # F-test
    report_progress("Finding SVGs")
    df.svg = svgTest(emat, fit.spanorm, fit.technical)
    SummarizedExperiment::rowData(spe) = cbind(SummarizedExperiment::rowData(spe), df.svg[rownames(spe), ])
    msgfun(sprintf("%d SVGs found (FDR < 0.05)", sum(df.svg$svg.fdr < 0.05)))

    spe
  }
)

fitSpaNormTechnical <- function(Y, fit.spanorm, msgfun) {
  msgfun(sprintf("%d cells/spots sampled to fit model", sum(fit.spanorm$sampling != "all")))

  # select technical covariates
  wtype = fit.spanorm$wtype
  W = fit.spanorm$W[, wtype != "biology", drop = FALSE]
  wtype = wtype[wtype != "biology"]

  # fit model
  fit.technical = fitSpaNormNB(
    Y,
    W,
    fit.spanorm$sampling != "all",
    maxn.psi = sum(fit.spanorm$sampling == "dispersion"),
    lambda.a = fit.spanorm$lambda.a,
    msgfun = msgfun
  )

  # create object
  fit.technical = SpaNormFit(
    ngenes = nrow(Y),
    ncells = ncol(Y),
    gene.model = fit.spanorm$gene.model,
    df.tps = fit.spanorm$df.tps,
    sample.p = fit.spanorm$sample.p,
    lambda.a =fit.spanorm$lambda.a,
    batch = fit.spanorm$batch,
    W = W,
    alpha = fit.technical$alpha,
    gmean = fit.technical$gmean,
    psi = fit.technical$psi,
    wtype = wtype,
    loglik = fit.technical$loglik,
    sampling = fit.technical$sampling
  )

  return(fit.technical)
}

svgTest <- function(Y, fit.spanorm, fit.technical) {
  # Add validation for input types
  if (!is.matrix(Y) && !methods::is(Y, "Matrix")) {
    stop("Y must be a matrix or Matrix object")
  }
  if (!methods::is(fit.spanorm, "SpaNormFit")) {
    stop("fit.spanorm must be a SpaNormFit object") 
  }
  if (!methods::is(fit.technical, "SpaNormFit")) {
    stop("fit.technical must be a SpaNormFit object")
  }

  # Existing dimension checks
  if (length(unique(c(nrow(Y), fit.spanorm$ngenes, fit.technical$ngenes))) != 1) {
    stop("number of genes differ between SpaNorm fits and/or data")
  }
  if (length(unique(c(ncol(Y), fit.spanorm$ncells, fit.technical$ncells))) != 1) {
    stop("number of cells differ between SpaNorm fits and/or data")
  }
  if (!all(colnames(fit.technical$W) %in% colnames(fit.spanorm$W)) | ncol(fit.technical$W) == ncol(fit.spanorm$W)) {
    stop("technical model is not nested in the full model")
  }

  # convert to regular matrix
  Y = as.matrix(Y)

  # nested model residual deviance
  # calculate mu
  mu = calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  # winsorize dispersion parameters (large disp slow the process)
  psi = fit.spanorm$psi
  psi.max = exp(median(log(psi)) + 3 * mad(log(psi)))
  psi = pmin(psi, psi.max)
  loglik.spanorm = rowSums(dnbinom(Y, mu = mu, size = 1 / psi, log = TRUE))

  # nested model residual deviance
  # calculate mu
  mu = calculateMu(fit.technical$gmean, fit.technical$alpha, fit.technical$W)
  # winsorize dispersion parameters (large disp slow the process)
  psi = fit.technical$psi
  psi.max = exp(median(log(psi)) + 3 * mad(log(psi)))
  psi = pmin(psi, psi.max)
  loglik.technical = rowSums(dnbinom(Y, mu = mu, size = 1 / psi, log = TRUE))

  # F-test
  df1 = ncol(fit.spanorm$W) - ncol(fit.technical$W)
  df2 = ncol(Y) - ncol(fit.spanorm$W)
  F.lrt = 2 * (loglik.spanorm - loglik.technical) / df1
  # Threshold to 0 due to convergence issues
  F.lrt = pmax(F.lrt, 0)
  p.val = pf(F.lrt, df1, df2, lower.tail = FALSE)
  fdr = p.adjust(p.val, method = "fdr")

  # create results table
  df.svg = data.frame(
    svg.F = F.lrt,
    svg.p = p.val,
    svg.fdr = fdr
  )

  return(df.svg)
}

#' Export top SVG results to a data frame
#' 
#' @param spe a SpatialExperiment object with SVG results from SpaNormSVG.
#' @param n a numeric, specifying the number of top SVGs to call.
#' @param fdr a numeric, specifying the false discovery rate (FDR) threshold for calling SVGs.
#' 
#' @return A data frame containing the top SVGs from F-test results including F-statistics, p-values and FDR.
#' @examples
#' 
#' library(SpatialExperiment)
#' library(ggplot2)
#'
#' data(HumanDLPFC)
#'
#' HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#' HumanDLPFC = SpaNormSVG(HumanDLPFC)
#' topSVGs = topSVGs(HumanDLPFC, n = 10)
#' @export
topSVGs <- function(spe, n = 10, fdr = 1) {
  stopifnot(n > 0)
  stopifnot(fdr >= 0 && fdr <= 1)
  checkSPE(spe)
  if (!all(c("svg.F", "svg.p", "svg.fdr") %in% colnames(SummarizedExperiment::rowData(spe)))) {
    stop("SVG results not found. Please run 'SpaNormSVG' first.")
  }
  
  cols = union(c("svg.F", "svg.p", "svg.fdr"), colnames(SummarizedExperiment::rowData(spe)))
  results = as.data.frame(SummarizedExperiment::rowData(spe)[, cols])
  results = results[order(results$svg.fdr), ]
  results = results[results$svg.fdr <= fdr, , drop = FALSE]
  n = min(n, nrow(results))
  results = results[1:n, , drop = FALSE]

  return(results)
}
