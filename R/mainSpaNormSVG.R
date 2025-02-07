#' @importFrom stats dnbinom qnbinom pnbinom mad median quantile model.matrix pf p.adjust
NULL

#' Model-based spatially variable gene (SVG) calling
#'
#' Performs normalisation of spatial transcriptomics data using spatially-dependent spot- and gene- specific size factors.
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

#' @rdname SpaNorm
setMethod(
  "SpaNormSVG",
  signature("SpatialExperiment"),
  function(spe, verbose) {
    checkSPE(spe)

    # message function depending on verbose param
    msgfun = ifelse(verbose, message, \(...){})

    # extract counts, coords, and size factors
    emat = SummarizedExperiment::assay(spe, "counts")

    # check whether previous results are present
    svg.cols = c("svg.F", "svg.p", "svg.fdr")
    if (any(svg.cols %in% colnames(SummarizedExperiment::rowData(spe)))) {
      warning("SVG results exist in 'spe' and will be overwritten")
      cols = setdiff(colnames(SummarizedExperiment::rowData(spe)), svg.cols)
      SummarizedExperiment::rowData(spe) = SummarizedExperiment::rowData(spe)[, cols]
    }

    # retrieve SpaNorm model
    fit.spanorm = S4Vectors::metadata(spe)$SpaNorm
    if (!is.null(fit.spanorm) &&
      fit.spanorm$ngenes == nrow(spe) &&
      fit.spanorm$ncells == ncol(spe)
    ) {
      msgfun("(1/3) Retrieve precomputed SpaNorm model")
    } else {
      stop("SVG calling requires a SpaNorm model. Please run 'SpaNorm()' on the `spe` object first.")
    }

    # fit nested model
    msgfun("(2/3) Fitting nested SpaNorm model")
    fit.technical = fitSpaNormTechnical(emat, fit.spanorm, msgfun)

    # F-test
    msgfun("(3/3) Finding SVGs")
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

