# Build the shared fitted R object ONCE per site build (not once per reference page). Every
# executed example on the R reference pages loads site/_cache/HumanDLPFC-fit.rds instead of
# refitting. Also writes the gene symbols of the subset (_cache/shared-genes.txt) so that the
# Python shared fit (build-shared-fit-py.py) uses exactly the same genes.
#
# The object: HumanDLPFC filtered with filterGenes(prop = 0.2), cut to the 300 most highly
# expressed genes (plus MOBP and ACTB, which the examples plot), with cheap settings
# (sample.p = 0.05, df.tps = 2, tol = 1e-2), SpaNormSVG results and a SpaNormPCA.
# Cached only inside one build (the file is git-ignored): delete _cache/ to force a refit.
if (dir.exists("_rlib")) .libPaths(c(normalizePath("_rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(SpaNorm)
  library(SpatialExperiment)
})
t0 <- Sys.time()
out <- "_cache/HumanDLPFC-fit.rds"
dir.create("_cache", showWarnings = FALSE)
data(HumanDLPFC)
spe <- HumanDLPFC[filterGenes(HumanDLPFC, prop = 0.2), ]
tot <- Matrix::rowSums(counts(spe))
must <- rowData(spe)$gene_name %in% c("MOBP", "ACTB")
keep <- order(-tot)[seq_len(300)]
keep <- sort(union(keep, which(must)))
spe <- spe[keep, ]
set.seed(36)
spe <- SpaNorm(spe, sample.p = 0.05, df.tps = 2, tol = 1e-2)
spe <- SpaNormSVG(spe)
spe <- SpaNormPCA(spe)
saveRDS(spe, out)
writeLines(rowData(spe)$gene_name, "_cache/shared-genes.txt")
message(sprintf("shared R fit: %d genes x %d spots, %.0f s -> %s", nrow(spe), ncol(spe),
                as.numeric(difftime(Sys.time(), t0, units = "secs")), out))
