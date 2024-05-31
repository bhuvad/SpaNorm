#' Human dorsolateral prefrontal cortex (DLPFC) visium sample
#'
#' Human DLPFC sample profiled using the 10x Visium platform. Lowly expressed genes were filtered out using the `filterGenes` function to retain genes expressed in at least 10% of samples. This version was obtained from the SpatialLIBD R/Bioconductor package.
#'
#' @format A SpatialExperiment containing region annotations in the `colData` column 'AnnotatedCluster'.
#' @docType data
#' @references Maynard KR, Collado-Torres L, Weber LM, Uytingco C, Barry BK, Williams SR, Catallini JL, Tran MN, Besich Z, Tippani M, Chew J. Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex. Nature neuroscience. 2021 Mar;24(3):425-36.
#' @source {SpatialLIBD} R/Bioconductor package.
"HumanDLPFC"