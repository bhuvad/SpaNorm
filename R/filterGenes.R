#' @importFrom methods is
#' @importFrom stats dnbinom mad median model.matrix pnbinom qnbinom quantile
NULL
 
#' @title  Filter genes based on expression
#' @description This function removes genes that are very lowly expressed.
#'
#' @param spe a SpatialExperiment, Seurat, or SpatialFeatureExperiment object containing count data.
#' @param prop a numeric, indicating the proportion of loci/cells where the gene should be expressed (default is 0.1, i.e., genes should be expressed in at least 10% of the loci/cells).
#'
#' @name filterGenes
#'
#' @return a logical vector encoding which genes should be kept for further analysis.
#' @examples
#' data(HumanDLPFC)
#' keep <- filterGenes(HumanDLPFC)
#' table(keep)
#'
#' @export
setGeneric("filterGenes", function(spe, prop = 0.1) standardGeneric("filterGenes"))

#' @rdname filterGenes
setMethod(
  "filterGenes",
  signature("SpatialExperiment", "ANY"),
  function(spe, prop) {
    stopifnot(prop >= 0 & prop <= 1)
    checkSPE(spe)
    emat = SummarizedExperiment::assay(spe, "counts")
    keep = filterGenes_intl(emat, prop)
    return(keep)
  }
)

#' @rdname filterGenes
setMethod(
  "filterGenes",
  signature("Seurat", "ANY"),
  function(spe, prop) {
    stopifnot(prop >= 0 & prop <= 1)
    checkSeurat(spe)
    emat = SeuratObject::LayerData(spe, "counts")
    keep = filterGenes_intl(emat, prop)
    return(keep)
  }
)

filterGenes_intl <- function(emat, prop) {
  keep = Matrix::rowMeans(emat > 0) >= prop
  return(keep)
}
