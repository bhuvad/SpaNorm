plotCovariates <- function(spe, covariate = c("biology", "ls", "batch"), ...) {
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
