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

- Exported
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md),
  [`nbProfilePsi()`](https://bhuvad.github.io/spaNorm/reference/nbProfilePsi.md)
  and the Newton solvers behind them
  ([`nbNewtonSolver()`](https://bhuvad.github.io/spaNorm/reference/nbNewtonSolver.md),
  [`nbNewtonSolverBatch()`](https://bhuvad.github.io/spaNorm/reference/nbNewtonSolver.md),
  [`nbGramBatch()`](https://bhuvad.github.io/spaNorm/reference/nbNewtonSolver.md),
  [`nbAbsorbGramBatch()`](https://bhuvad.github.io/spaNorm/reference/nbNewtonSolver.md),
  [`nbMuFloor()`](https://bhuvad.github.io/spaNorm/reference/nbMuFloor.md)),
  moved here from spiDE’s per-gene convergence stage.
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)
  takes an already-fitted
  [`fitNB()`](https://bhuvad.github.io/spaNorm/reference/fitNB.md)-style
  model (design, coefficients, dispersion) and converges every gene
  separately to its own penalised negative binomial optimum by damped
  Newton, gene-blocked and dispatched via `BiocParallel`;
  [`nbProfilePsi()`](https://bhuvad.github.io/spaNorm/reference/nbProfilePsi.md)
  re-estimates just the dispersion by profile maximum likelihood at
  fixed coefficients. Both take an `offset` (a known effect added to the
  linear predictor with its coefficient fixed at 1, rather than a
  floating design column that could absorb correlated signal);
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)
  also takes `absorb`/`absorb.batch` (a wide indicator block – such as
  per-sample random-effect columns – folded into the per-gene solve by a
  Schur complement, at the cost of one dense-column gram per iteration
  regardless of block width). spiDE calls these for its own per-gene
  convergence stage.

- Added
  [`polishSpaNorm()`](https://bhuvad.github.io/spaNorm/reference/polishSpaNorm.md)
  and
  [`isPolished()`](https://bhuvad.github.io/spaNorm/reference/isPolished.md).
  [`polishSpaNorm()`](https://bhuvad.github.io/spaNorm/reference/polishSpaNorm.md)
  maps a fitted
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
  (or
  [`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md)
  null) model onto
  [`polishNB()`](https://bhuvad.github.io/spaNorm/reference/polishNB.md)‘s
  generic problem – the shared library-size coefficient `a1` pulled out
  of the design as a per-cell offset, the per-gene mean `gmean` an
  unpenalised intercept column, and the biology/library-size/batch
  columns ridge-penalised as
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
  penalises them (`lambda.a` times the number of cells/spots, by column
  `wtype`) – converges every gene, and rewrites the normalised assay
  from the polished fit. `psi.method` sets how each gene’s dispersion is
  handled at the converged mean (`"fixed"`, the default, keeps
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)’s
  value; `"profile"` re-estimates it by profile ML and re-polishes);
  `ls` sets whether the shared `a1` is held fixed (`"fixed"`, the
  default) or re-estimated jointly with the per-gene coefficients at the
  optimum of the total penalised likelihood over the polished genes
  (`"joint"`); `cells` chooses all cells/spots or only the ones
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
  sampled to fit; a gene with no counts in the polished cells, or with
  none in the polished cells of one batch level (a level of the
  unpenalised 0/1 batch columns), is held out and keeps its input fit
  (`polished = FALSE`, the reason in the diagnostics’ `held_out` column,
  and a warning counting the batch-level cases) rather than driven to an
  unbounded coefficient; a warning also counts any gene the Newton
  engine could not polish; and `overwrite` allows re-polishing an
  already-polished fit from the unpolished one kept alongside it.
  [`isPolished()`](https://bhuvad.github.io/spaNorm/reference/isPolished.md)
  reports whether a stored fit has been polished. Measured on four real
  cores (YTMA CosMx WTA, 948-16,350 genes x 395-2,645 cells,
  `psi.method = "fixed"` only): `ls = "joint"` moved the shared `a1` by
  7.9-79.0 of its own profiled standard error and took the SVG calls
  (FDR \< 0.05) from 137 to 123 on the 395-cell core, and from 5 to 5,
  180 to 178 and 545 to 551 on the other three. On two of the cores
  (10,738 x 395 and 948 x 646), over a random 50 genes each, it shifted
  the genes’ fitted log-means by a median 0.09 and 0.03 of their own
  standard error (max 0.68 and 0.17), concentrated in the
  lowest-library-size cells. The joint polish took 4-7x the wall time of
  a fixed one in that comparison, which ran with an earlier stop on `a1`
  (`tol.ls = 1e-6`) under which every joint fit took all 10 steps; the
  shipped stop (`1e-3`) ends sooner, and the cost has not been
  re-measured. `ls = "fixed"` is the default because
  [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)’s
  own `a1` is an unweighted mean over genes where the joint estimate is
  information-weighted and so dominated by the brightest ones, changing
  the estimand rather than only its precision. That comparison covers
  four cores only, with no ground truth for which `a1` is correct and no
  independent validation of the resulting SVG-call changes.

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

- [`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md)
  now fits, retrieves or re-polishes the null (technical-only) model so
  that it is paired with the full fit’s polish state AND settings
  (`psi.method`, `ls`, `cells`, and under `psi.method = "profile"` the
  dispersion’s search interval `psi.range`), not only whether each is
  [`isPolished()`](https://bhuvad.github.io/spaNorm/reference/isPolished.md):
  two fits can both be “polished” at different settings, which changes
  the objective just as much as one being unpolished. `svgTest()`
  carries the same check and refuses (rather than silently scoring) a
  full/null pair that is not paired this way. Either polish the full
  model before calling
  [`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md)
  (which then polishes the null to match) or call
  [`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md)
  directly and let it pair the null itself.

### Bug Fixes

- Fixed the null (technical-only) model fitted by
  [`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md)
  being penalised less than the full model. `fitSpaNorm()` scales the
  ridge penalty `lambda.a` by the number of cells/spots before fitting,
  but `fitSpaNormTechnical()` rebuilt the penalty from the stored,
  unscaled `lambda.a` and did not rescale it, so the nested null was
  under-penalised by a factor of `ncol(Y)` on the library-size terms.
  The two nested fits now use the same penalty, making the
  likelihood-ratio test consistent. Because the less-penalised null
  generally fitted slightly better, `svg.F` was slightly deflated for
  most genes; SVG statistics will typically increase slightly after
  refitting the null model, although individual genes can move in either
  direction (a cached `SpaNormNull` fit in `metadata(spe)` is reused
  as-is, so remove it to refit). `fitSpaNormTechnical()` now also stops
  if the data and the full fit have different numbers of cells.

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
