# Changelog

## SpaNorm 1.8.0

### New Features

- [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
  gains a `BPPARAM` argument to parallelise normalisation across workers
  via `BiocParallel`, accelerating large datasets. It defaults to
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  (no parallelisation), and results are identical regardless of the
  backend used.
- [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
  now normalises `DelayedArray`-backed count assays (e.g. disk-backed
  via `HDF5Array`) block-wise, so out-of-core datasets are processed
  without ever loading the full matrix into memory. Results match the
  in-memory path.

### Improvements

- The optional GPU backend now uses the `torch` package instead of
  TensorFlow, adding native support for NVIDIA CUDA and Apple Silicon
  (Metal/MPS) devices and removing the Python/reticulate dependency.
  Users of `backend = "gpu"` should install `torch` in place of
  `tensorflow`.

## SpaNorm 1.2.0

- Added model-based spatially variable gene (SVG) calling.
- Added spatial visualisation funciton `plotSpatial` to visualise
  colData, gene expression, and reduced dimensions.
- Added spatial visualisation function `plotCovariate` to visualise the
  biolgy, batch, and library size functions estimated by SpaNorm.
- Dynamic calculation of df.tps for rectangular tissue sections.
- Allow separate specification of df.tps for biology and library size.
- Added GLM-PCA approximation through the `SpaNormPCA` function. The
  null model is considered to consist of the library size effects, batch
  effects, and the gene mean.

## SpaNorm 1.0.0

- Initial Bioconductor submission.
