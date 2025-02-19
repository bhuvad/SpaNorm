#' @importFrom BiocSingular bsparam
#' @importFrom BiocParallel SerialParam
NULL

#' GLM-based (SpaNorm) PCA
#'
#' GLM-based PCA using the SpaNorm model. The null model is considered to consist of the library size effects, batch effects, and the gene mean. GLM-PCA is approximated by regressing the null model from the data, and performing PCA on the residuals (Pearson or deviance).
#'
#' @param spe a SpatialExperiment or Seurat object, with the count data stored in 'counts' or 'data' assays respectively, and a SpaNorm model fit.
#' @param nsvgs the number of SVGs to use for PCA.
#' @param ncomponents the number of components to compute.
#' @param svg.fdr the FDR threshold for SVG calling.
#' @param BSPARAM a BiocSingularParam object specifying which algorithm should be used to perform the PCA.
#' @param BPPARAM a BiocParallelParam object specifying whether the PCA should be parallelized.
#' @param name the name of the reducedDim to store the PCA results.
#' @param residuals the type of residuals to use for PCA. Either "deviance" (default) or "pearson".
#' 
#' @details SpaNorm PCA works by using the SpaNorm model fit for data normalisation to approximate a GLM-based PCA as described in Townes et al. (Genome Biology, 2019). The model used for normalisation represents the library size effects and the gene mean. Regressing these covariates, we remain with the deviance or Pearson residuals, upon which PCA can be performed to approximate the GLM-PCA.
#'
#' @return a SpatialExperiment or Seurat object with PCA results. For SpatialExperiment objects, these are stored in the reducedDims.
#' @name SpaNormPCA
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
#' HumanDLPFC = SpaNormPCA(HumanDLPFC)
#' reducedDims(HumanDLPFC)
#' 
#' @export
#'
setGeneric("SpaNormPCA", function(
    spe,
    nsvgs = 3000,
    ncomponents = 50,
    svg.fdr = 1,
    BSPARAM = bsparam(),
    BPPARAM = SerialParam(),
    residuals = c("deviance", "pearson"),
    name = "PCA"
    ) {
  standardGeneric("SpaNormPCA")
})

#' @rdname SpaNormPCA
setMethod(
  "SpaNormPCA",
  signature("SpatialExperiment"),
  function(spe, nsvgs, ncomponents, svg.fdr, BSPARAM, BPPARAM, residuals, name) {
    checkSPE(spe)
    stopifnot(nsvgs > 0, ncomponents > 0, ncomponents <= ncol(spe))
    residuals = match.arg(residuals)

    # Get and validate SpaNorm model
    fit.spanorm = getSpaNormFit(spe)
    fit.technical = getSpaNormFit(spe, null = TRUE)
    # get SVGs
    svgs = rownames(topSVGs(spe, nsvgs, svg.fdr))

    # compute residuals
    emat = SummarizedExperiment::assay(spe, "counts")
    emat = getResiduals(emat, fit.technical, type = residuals)[svgs, , drop = FALSE]

    # add PCA results to reducedDims
    SingleCellExperiment::reducedDim(spe, name) = calculatePCA(emat, ncomponents, BSPARAM, BPPARAM)

    spe
  }
)

getResiduals <- function(emat, fit.technical, type = c("deviance", "pearson")) {
  type = match.arg(type)
  if (type == "pearson") {
    emat = normalisePearson(emat, 1, fit.technical)
  } else {
    emat = devianceResiduals(emat, fit.technical)
  }
  return(emat)
}

devianceResiduals <- function(Y, fit.technical) {
  # calculate mu
  mu <- calculateMu(fit.technical$gmean, fit.technical$alpha, fit.technical$W)
  mu <- as(mu, "sparseMatrix")
  # winsorize dispersion parameters (large disp slow the process)
  psi <- fit.technical$psi
  psi.max <- exp(median(log(psi)) + 3 * mad(log(psi)))
  psi <- pmin(psi, psi.max)

  # calculate deviance residuals
  dev <- Y * log(Y / mu)
  dev <- dev - (Y + 1 / psi) * log((1 + Y * psi) / (1 + mu * psi))
  dev[Y == 0] <- (-1 / psi * log(1 / (1 + mu * psi)))[Y == 0]
  dev <- sign(Y - mu) * sqrt(2 * dev)

  return(dev)
}

calculatePCA <- function(emat, ncomponents, BSPARAM, BPPARAM) {
  # compute PCA
  res = BiocSingular::runPCA(
    Matrix::t(emat),
    rank = ncomponents,
    center = TRUE,
    BSPARAM = BSPARAM,
    BPPARAM = BPPARAM
  )
  res$var = res$sdev^2
  
  # create attributes
  pca_attr = list(
    percentVar = res$var / sum(res$var) * 100,
    varExplained = res$var,
    rotation = res$rotation
  )

  # create results
  res = res$x
  attributes(res) = c(attributes(res), pca_attr)

  return(res)
}