# SpaNorm 1.5.4

* Added a `BPPARAM` argument to `SpaNorm()` to parallelise the normalisation step over gene-blocks via `BiocParallel`, speeding up the (expensive) logpac transform on large datasets. Defaults to `SerialParam()` (no parallelisation); results are identical regardless of blocking.
* Sparse count inputs are kept sparse throughout model fitting and normalisation (only the sampled sub-matrix and one gene-block at a time are ever densified), and the logpac transform now densifies the counts once rather than twice, trimming a redundant full-matrix allocation. Results are unchanged.

# SpaNorm 1.5.3

* Switched the optional GPU backend from TensorFlow to the torch package, removing the Python/reticulate dependency. GPU acceleration now runs on NVIDIA CUDA and Apple Silicon (Metal/MPS) devices, falling back to CPU when neither is available. Users of `backend = "gpu"` should install `torch` in place of `tensorflow`.

# SpaNorm 1.2.0

* Added model-based spatially variable gene (SVG) calling.
* Added spatial visualisation funciton `plotSpatial` to visualise colData, gene expression, and reduced dimensions.
* Added spatial visualisation function `plotCovariate` to visualise the biolgy, batch, and library size functions estimated by SpaNorm.
* Dynamic calculation of df.tps for rectangular tissue sections.
* Allow separate specification of df.tps for biology and library size.
* Added GLM-PCA approximation through the `SpaNormPCA` function. The null model is considered to consist of the library size effects, batch effects, and the gene mean.

# SpaNorm 1.0.0

* Initial Bioconductor submission.
