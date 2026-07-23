# SpaNorm 1.8.0

## New Features

* `SpaNorm()` gains a `BPPARAM` argument to parallelise normalisation across workers via `BiocParallel`, accelerating large datasets. It defaults to `BiocParallel::SerialParam()` (no parallelisation), and results are identical regardless of the backend used.
* `SpaNorm()` now normalises `DelayedArray`-backed count assays (e.g. disk-backed via `HDF5Array`) block-wise, so out-of-core datasets are processed without ever loading the full matrix into memory. Results match the in-memory path.
* Exported `fitNB()`, which fits a per-gene negative binomial GLM over an arbitrary design matrix using SpaNorm's IRLS engine (with optional ridge regularisation and adjustable outlier winsorisation). This exposes the model-fitting machinery for reuse independently of SpaNorm's spatial model.

## Improvements

* The optional GPU backend now uses the `torch` package instead of TensorFlow, adding native support for NVIDIA CUDA and Apple Silicon (Metal/MPS) devices and removing the Python/reticulate dependency. Users of `backend = "gpu"` should install `torch` in place of `tensorflow`.
* The GPU backend now automatically detects available accelerator memory and fits large datasets in gene-blocks so peak GPU memory stays bounded, avoiding out-of-memory failures on GPUs with limited VRAM. This requires no additional arguments; the detected budget can be overridden via the new `gpu.mem.budget` parameter, and results match the CPU backend within a small numerical tolerance.
* The dispersion winsorisation used during normalisation now clamps at 4 MAD (previously 3), matching the coefficient and mean winsorisation, and is configurable via the winsorisation controls on the fitting/normalisation helpers.

# SpaNorm 1.2.0

* Added model-based spatially variable gene (SVG) calling.
* Added spatial visualisation funciton `plotSpatial` to visualise colData, gene expression, and reduced dimensions.
* Added spatial visualisation function `plotCovariate` to visualise the biolgy, batch, and library size functions estimated by SpaNorm.
* Dynamic calculation of df.tps for rectangular tissue sections.
* Allow separate specification of df.tps for biology and library size.
* Added GLM-PCA approximation through the `SpaNormPCA` function. The null model is considered to consist of the library size effects, batch effects, and the gene mean.

# SpaNorm 1.0.0

* Initial Bioconductor submission.
