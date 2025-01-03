#' Diagnostic plot of predicted expression for a covariate
#' 
#' This function can be used to spatially visualise the library size, biology or batch specific effect modelled for each gene.
#' 
#' @param spe a SpatialExperiment object.
#' @param covariate a character, specifying the type of covariate to be plot: "biology" (default), "ls" to plot the library size effect, and "batch" to plot the batch-specific effect.
#' @param ... additional parameters to be passed to the [plotSpatial] function.
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' library(ggplot2)
#' 
#' data(HumanDLPFC)
#' \donttest{
#' HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#' # plot spatial region annotations
#' p1 <- plotCovariate(HumanDLPFC, covariate = "biology", colour = ENSG00000075624) +
#'   scale_colour_viridis_c(option = "F")
#' p1
#' 
#' p2 <- plotCovariate(HumanDLPFC, covariate = "ls", colour = ENSG00000075624) +
#'   scale_colour_viridis_c(option = "F")
#' p2
#' }
#'
plotCovariate <- function(spe, covariate = c("biology", "ls", "batch"), ...) {
  checkSPE(spe)
  covariate = match.arg(covariate)

  #----extract aes----
  aesmap = rlang::enquos(...)
  # compute plot
  aesmap = aesmap[!names(aesmap) %in% c("what", "assay", "dimred")]
  pltcols = as.character(sapply(aesmap, rlang::quo_get_expr))
  genes = intersect(pltcols, rownames(spe))

  if (length(genes) == 0) {
    stop("no genes selected for the plot")
  }

  #----compute adjustment----
  if (!"SpaNorm" %in% names(S4Vectors::metadata(spe))) {
    stop("SpaNorm fit not found in the SpatialExperiment object")
  }
  fit = S4Vectors::metadata(spe)$SpaNorm

  # check for batch covariates if plotting
  if (covariate == "batch" & !"batch" %in% fit@wtype) {
    stop("no batch covariates found in the SpaNorm fit")
  }

  # compute adjustment
  normaliseFun = switch(covariate,
    biology = normaliseMeanBio,
    ls = normaliseMeanLS,
    batch = normaliseMeanBatch
  )
  Y = SingleCellExperiment::counts(spe)
  SingleCellExperiment::logcounts(spe) = normaliseFun(Y, 1, fit)

  # plot
  p1 = plotSpatial(spe, what = "expression", assay = "logcounts", !!!aesmap) +
    ggplot2::labs(title = covariate)  

  return(p1)
}
