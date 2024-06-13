#' @title  Filter genes based on expression
#' @description This function computes the size factors using a fast but inaccurate approach. Size factors are computed using the direct estimate of library sizes (sum of all counts). Though fast, this approach does not cater for compositional biases in the data and therefore is less accurate than scran-based estimates.
#'
#' @param spe a SpatialExperiment, Seurat, or SpatialFeatureExperiment object containing count data.
#'
#' @name fastSizeFactors
#'
#' @return a SpatialExperiment, Seurat, or SpatialFeatureExperiment, containing size factors in the 'sizeFactor' column of the column annotation.
#' @examples
#' data(HumanDLPFC)
#' HumanDLPFC <- fastSizeFactors(HumanDLPFC)
#' head(HumanDLPFC$sizeFactor)
#'
#' @export
setGeneric("fastSizeFactors", function(spe) standardGeneric("fastSizeFactors"))

setMethod(
  "fastSizeFactors",
  signature("SpatialExperiment"),
  function(spe) {
    checkSPE(spe)
    emat <- SummarizedExperiment::assay(spe, "counts")
    SingleCellExperiment::sizeFactors(spe) <- fastSizeFactors_intl(emat)
    return(spe)
  }
)

#' @rdname filterGenes
fastSizeFactors_intl <- function(Y) {
  LS = Matrix::colSums(Y)
  LS = LS / mean(LS)
  return(LS)
}