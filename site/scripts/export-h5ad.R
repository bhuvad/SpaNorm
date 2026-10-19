# Export the HumanDLPFC example data (gene symbols as names) to h5ad so the Python
# tabs operate on exactly the same data as the R tabs. Run in its own R process
# because zellkonverter initialises its own (basilisk) Python.
args <- commandArgs(trailingOnly = TRUE)
out_file <- args[[1]]
if (dir.exists("_rlib")) .libPaths(c(normalizePath("_rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(SpaNorm)
  library(SpatialExperiment)
})
data(HumanDLPFC, package = "SpaNorm")
se <- HumanDLPFC
rownames(se) <- rowData(se)$gene_name
out <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = as.matrix(counts(se))),
  colData = S4Vectors::DataFrame(
    AnnotatedCluster = as.character(colData(se)$AnnotatedCluster),
    cell_count = colData(se)$cell_count,
    row.names = colnames(se)
  )
)
SingleCellExperiment::reducedDim(out, "spatial") <- as.matrix(spatialCoords(se))
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
zellkonverter::writeH5AD(out, out_file, X_name = "counts", verbose = FALSE)
