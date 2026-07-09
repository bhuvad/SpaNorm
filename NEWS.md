# SpaNorm 1.5.5

* `SpaNorm()` now detects a `DelayedArray`-backed counts assay (e.g. disk-backed via `HDF5Array`) and normalises it block-wise so the full matrix is never realised in memory at once, capping peak memory for out-of-core datasets. Results are identical to the in-memory path. Block size follows `DelayedArray::getAutoBlockSize()`.

# SpaNorm 1.5.4

* Added a `BPPARAM` argument to `SpaNorm()` to parallelise the normalisation step over gene-blocks via `BiocParallel`, speeding up the (expensive) logpac transform on large datasets. Defaults to `SerialParam()` (no parallelisation); results are identical regardless of blocking.
* Sparse count inputs are kept sparse throughout model fitting and normalisation (only the sampled sub-matrix and one gene-block at a time are ever densified), and the logpac transform now densifies the counts once rather than twice, trimming a redundant full-matrix allocation. Results are unchanged.
* Fixed a severe GPU (torch/MPS) slowdown in the IRLS model fitting: device tensors created each iteration were not being reclaimed (R's garbage collector does not see GPU memory pressure), so the MPS allocator degraded and iterations grew progressively slower. The IRLS loops now release device memory each iteration, keeping per-iteration time flat, and the first regression coefficient's mean is reduced on-device instead of copying the whole coefficient matrix to the host each iteration.

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
