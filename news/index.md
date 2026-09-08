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

- Exported
  [`fitNB()`](https://bhuvad.github.io/spaNorm/reference/fitNB.md),
  which fits a per-gene negative binomial GLM over an arbitrary design
  matrix using SpaNorm’s IRLS engine (with optional ridge regularisation
  and adjustable outlier winsorisation). This exposes the model-fitting
  machinery for reuse independently of SpaNorm’s spatial model.

- [`fitNB()`](https://bhuvad.github.io/spaNorm/reference/fitNB.md) and
  [`calculateMu()`](https://bhuvad.github.io/spaNorm/reference/calculateMu.md)
  gain an `offset` argument: a genes x cells matrix added to the linear
  predictor with its coefficient fixed at 1, so
  `log(mu) = gmean + tcrossprod(alpha, W) + offset`. Use it for an
  effect that is already known rather than adding a column to the
  design, whose coefficient would float and could absorb signal
  correlated with that effect. The offset is subset alongside the
  counts, applied to the dispersion estimation as well as the mean, and
  sliced in step with the counts when the fit is carved into
  gene-blocks. It defaults to `NULL`, which is a strict no-op, so all
  existing behaviour is unchanged. The GPU path is implemented
  backend-agnostically but has not yet been executed on an accelerator.

- [`fitNB()`](https://bhuvad.github.io/spaNorm/reference/fitNB.md) also
  gains a `psi` argument for supplying per-gene NB dispersions instead
  of estimating them. Supplied dispersions are used as-is (no
  re-estimation, no winsorisation) and the outer dispersion loop is
  bypassed: the coefficients come from a single IRLS fit at the given
  `psi`, so
  [`edgeR::estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)
  is never called and the returned `sampling` factor carries no
  `"dispersion"` level. The dispersions need not come from an identical
  design – values estimated on a design nested in (or equal to) the
  fitting design are appropriate, e.g. pooled across a coarser model,
  which errs conservative. `psi` composes with `offset`, and
  `psi = NULL` (the default) is a strict no-op.

- Exported the quasi-likelihood dispersion machinery of edgeR v4 as
  generic NB-GLM functions with a CPU and a torch backend:
  [`nbUnitDeviance()`](https://bhuvad.github.io/spaNorm/reference/nbUnitDeviance.md)
  (the NB unit deviance, elementwise over genes x cells),
  [`nbDevianceMoments()`](https://bhuvad.github.io/spaNorm/reference/nbDevianceMoments.md)
  (the unit deviance’s first two moments under the fitted NB by direct
  summation over the pmf, the definition edgeR’s Chebyshev tables
  approximate) and
  [`qlDispersion()`](https://bhuvad.github.io/spaNorm/reference/qlDispersion.md)
  (per-gene adjusted deviance, effective residual df and their ratio).
  The moments depend only on (mu, phi), so
  [`qlDispersion()`](https://bhuvad.github.io/spaNorm/reference/qlDispersion.md)
  evaluates them once on a shared (log mu, log phi) table and
  interpolates onto every gene and cell; the per-gene work is then
  elementwise and runs on the accelerator when the counts and means are
  torch tensors. Oracle-tested against
  [`edgeR::glmQLFit()`](https://rdrr.io/pkg/edgeR/man/glmQLFit.html);
  used by spiDE for the standard-error scale of its niche tests.

### Improvements

- The optional GPU backend now uses the `torch` package instead of
  TensorFlow, adding native support for NVIDIA CUDA and Apple Silicon
  (Metal/MPS) devices and removing the Python/reticulate dependency.
  Users of `backend = "gpu"` should install `torch` in place of
  `tensorflow`.
- The GPU backend now automatically detects available accelerator memory
  and fits large datasets in gene-blocks so peak GPU memory stays
  bounded, avoiding out-of-memory failures on GPUs with limited VRAM.
  This requires no additional arguments; the detected budget can be
  overridden via the new `gpu.mem.budget` parameter, and results match
  the CPU backend within a small numerical tolerance.
- The dispersion winsorisation used during normalisation now clamps at 4
  MAD (previously 3), matching the coefficient and mean winsorisation,
  and is configurable via the winsorisation controls on the
  fitting/normalisation helpers.
- Exported the generic GPU device/memory-budget/tensor-conversion layer
  ([`checkGPU()`](https://bhuvad.github.io/spaNorm/reference/checkGPU.md),
  [`getBackendDevice()`](https://bhuvad.github.io/spaNorm/reference/getBackendDevice.md),
  [`getBackendDtype()`](https://bhuvad.github.io/spaNorm/reference/getBackendDtype.md),
  [`gpuDtypeBytes()`](https://bhuvad.github.io/spaNorm/reference/gpuDtypeBytes.md),
  [`getGPUMemoryBudget()`](https://bhuvad.github.io/spaNorm/reference/getGPUMemoryBudget.md),
  [`resetGPUCache()`](https://bhuvad.github.io/spaNorm/reference/resetGPUCache.md),
  [`is_torch_tensor()`](https://bhuvad.github.io/spaNorm/reference/is_torch_tensor.md),
  [`toGPUMatrix()`](https://bhuvad.github.io/spaNorm/reference/toGPUMatrix.md),
  [`toGPUVector()`](https://bhuvad.github.io/spaNorm/reference/toGPUVector.md),
  [`toRMatrix()`](https://bhuvad.github.io/spaNorm/reference/toRMatrix.md),
  [`diag_mat()`](https://bhuvad.github.io/spaNorm/reference/diag_mat.md),
  [`tcrossprod_gpu()`](https://bhuvad.github.io/spaNorm/reference/tcrossprod_gpu.md),
  [`matmul_gpu()`](https://bhuvad.github.io/spaNorm/reference/matmul_gpu.md),
  [`add_vec_mat_gpu()`](https://bhuvad.github.io/spaNorm/reference/add_vec_mat_gpu.md),
  [`mult_vec_mat_gpu()`](https://bhuvad.github.io/spaNorm/reference/mult_vec_mat_gpu.md),
  [`rowSums_gpu()`](https://bhuvad.github.io/spaNorm/reference/rowSums_gpu.md),
  [`dnbinom_gpu()`](https://bhuvad.github.io/spaNorm/reference/dnbinom_gpu.md),
  [`hasBadValues()`](https://bhuvad.github.io/spaNorm/reference/hasBadValues.md))
  and added a new
  [`invert_mat_batched()`](https://bhuvad.github.io/spaNorm/reference/invert_mat_batched.md),
  so downstream packages (e.g. spiDE) can build their own GPU-blocked
  per-gene computation on the same device/dtype/memory-budget machinery,
  without depending on unexported internals.
  [`calculateMu()`](https://bhuvad.github.io/spaNorm/reference/calculateMu.md)
  now accepts a `backend` argument (`"cpu"` by default, unchanged
  behaviour) and dispatches to the accelerator when requested.
- Fixed
  [`getGPUMemoryBudget()`](https://bhuvad.github.io/spaNorm/reference/getGPUMemoryBudget.md)’s
  CUDA auto-detection reporting the whole physical GPU’s free memory
  rather than the process’s assigned MIG instance’s, causing a many-fold
  budget overestimate (and downstream out-of-memory errors) on
  MIG-partitioned GPUs; it now resolves the correct device via
  `CUDA_VISIBLE_DEVICES`. Added
  [`setGPUMemoryBudget()`](https://bhuvad.github.io/spaNorm/reference/setGPUMemoryBudget.md)
  to explicitly set (and cache) the budget for the session, for cases
  where auto-detection remains unreliable.

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
