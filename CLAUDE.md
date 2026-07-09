# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this is

SpaNorm is a Bioconductor R package for spatially-aware library-size
normalisation of spatial transcriptomics data. It fits a gene-wise
negative binomial (NB) regression whose covariates are spatially-smooth
thin-plate-spline (TPS) functions, decomposes variation into
library-size (technical) vs. library-size-independent (biological)
components, and produces adjusted data (default: log percentile-
adjusted counts, “logpac”).

## Common commands

This is a pure-R package with no compiled code (`src/` is empty). Run
from the package root.

``` r

# Build & install locally
devtools::document()          # regenerate NAMESPACE + man/*.Rd from roxygen (roxygen2 v8.0.0)
devtools::load_all()          # load without installing (fastest dev loop)
devtools::install()

# Tests (testthat 3e)
devtools::test()                                    # all tests
testthat::test_file("tests/testthat/test-mainSpaNorm.R")   # single file
devtools::test(filter = "fitSpaNormNB")             # files matching test-<filter>.R

# Checks (must pass before a Bioconductor push)
devtools::check()             # R CMD check
BiocCheck::BiocCheck()        # Bioconductor-specific checks

# Vignette / docs
devtools::build_vignettes()
pkgdown::build_site()
```

Prefer the **`deploy-bioc` skill** for the full pre-push ritual
(document → vignette → review → R CMD check → BiocCheck → version bump →
commit/push). CI (`.github/workflows/R-CMD-check.yaml`) runs R CMD check
on macOS/Windows/Ubuntu (release + devel).

## Architecture

### Public API (see `NAMESPACE`)

All entry points are S4 generics with methods for **both
`SpatialExperiment` and `Seurat`** objects (and some for
`SpatialFeatureExperiment`). The container-specific methods extract
counts/coords/size factors, delegate to a shared internal
implementation, then write results back into the container:

- [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md) —
  fit the model and write adjusted values (`logcounts` assay / Seurat
  `data` layer).
- [`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md)
  — spatially-variable-gene calling via a likelihood-ratio test (full
  vs. biology-free nested model).
- [`SpaNormPCA()`](https://bhuvad.github.io/spaNorm/reference/SpaNormPCA.md)
  — GLM-PCA approximation on Pearson/deviance residuals of the null
  (LS + batch + gene mean) model.
- [`filterGenes()`](https://bhuvad.github.io/spaNorm/reference/filterGenes.md),
  [`fastSizeFactors()`](https://bhuvad.github.io/spaNorm/reference/fastSizeFactors.md),
  [`plotSpatial()`](https://bhuvad.github.io/spaNorm/reference/plotSpatial.md),
  [`plotCovariate()`](https://bhuvad.github.io/spaNorm/reference/plotCovariate.md),
  [`topSVGs()`](https://bhuvad.github.io/spaNorm/reference/topSVGs.md).

### Core fit pipeline (`R/mainSpaNorm.R`)

[`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
method → **`.spaNormCore()`** (shared by SPE + Seurat methods) → decides
cache-reuse vs. refit → **`fitSpaNorm()`** → **`fitSpaNormNB()`**
(`R/fitSpaNormNB.R`, the IRLS engine).

- `fitSpaNorm()` builds the design matrix `W` from: `logLS`, biology TPS
  basis, `logLS × LS-TPS` interaction, and optional `batch` columns.
  Each column is tagged in **`wtype`** as `"biology"`, `"ls"`, or
  `"batch"` — this labelling is what lets the adjustment/SVG/PCA steps
  selectively zero-out technical variation.
- TPS bases are built by `bs.tps()` (tensor product of
  [`splines::ns`](https://rdrr.io/r/splines/ns.html) along scaled x/y).
  `df.tps` may differ for biology vs. LS; for rectangular tissue the
  per-axis df is derived from the aspect ratio.
- Size factors: if absent, computed internally via
  [`scran::quickCluster`](https://rdrr.io/pkg/scran/man/quickCluster.html) +
  `calculateSumFactors`.
- The result is a **`SpaNormFit`** S4 object (`R/AllClasses.R`) with a
  `$` accessor and a `setValidity` invariant. It stores `W`,
  coefficients `alpha`, `gmean`, dispersion `psi`, `wtype`, `loglik`,
  `sampling`.

### Fit caching

The `SpaNormFit` is stored in `metadata(spe)$SpaNorm` (SPE) or
`spe@misc$SpaNorm` (Seurat) and **reused on subsequent calls** unless
`overwrite=TRUE` or a fit-defining parameter changed (`df.tps`,
`lambda.a`, `batch`, `gene.model`, dimensions) — see the guard in
`.spaNormCore()` and `matchDftps`/`matchLambda`.
`SpaNormSVG`/`SpaNormPCA` retrieve this cached fit (via `getSpaNormFit`)
rather than refitting.

### Adjustment / normalisation (`R/fitSpaNormNB.R`)

`getAdjustmentFun()` maps `adj.method` → one of `normaliseLogPAC`
(default/“auto”), `normalisePearson`, `normaliseMedianBio`,
`normaliseMeanBio`. These are **row-wise (per-gene)** transforms.

`normaliseBlocked()` applies the chosen transform over contiguous
gene-blocks and `rbind`s the results — identical to the whole-matrix
result. Blocking engages when either: - input counts are a
**`DelayedArray`** (disk-backed, e.g. HDF5Array) → finite `block.size`
keeps the array out-of-core (block size follows
[`DelayedArray::getAutoBlockSize()`](https://rdrr.io/pkg/DelayedArray/man/AutoBlock-global-settings.html)),
or - a multi-worker **`BPPARAM`** is supplied → blocks run in parallel
via
[`BiocParallel::bplapply`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html).
In-memory + serial takes the direct whole-matrix path
(`block.size = Inf`). The dispersion-winsorisation threshold is a
*global* statistic, so `subsetFitGenes()` carries it as a `psi.max`
attribute per block.

### GPU backend (`R/gpuFunctions.R`)

Optional acceleration via the **`torch`** package (not TensorFlow —
migrated in 1.8.0; no reticulate/Python). `backend = "auto"|"cpu"|"gpu"`
on `SpaNorm`/`SpaNormSVG`. Device is resolved once per session and
cached in the `.spanorm_env` environment (`resetGPUCache()` clears it,
used in tests): - Device precedence: CUDA → MPS (Apple Silicon) → CPU
(`getBackendDevice`). - **dtype differs by device**: CUDA/CPU use
float64 for exact CPU parity; **MPS cannot do float64 and uses float32**
(`getBackendDtype`). Expect small numerical differences on MPS. - The
IRLS math is written against `*_gpu` helpers (`tcrossprod_gpu`,
`dnbinom_gpu`, `colSums_gpu`, …) that operate on either base-R matrices
or `torch_tensor`s. **Watch for mixing base-R and tensor objects** — use
`toGPUMatrix`/`toRMatrix`/`is_torch_tensor` at the boundaries.
`dnbinom_gpu`/`pnbinom`/`qnbinom` are real torch implementations of the
NB density.

## Conventions

- **S4 argument flow**: extra fitting parameters (e.g. `maxn.psi`,
  `step.factor`, `maxit.*`) are declared on the `SpaNorm` generic but
  ride through `...` and are unpacked only where consumed
  (`fitSpaNormNB`). When adding a tuning parameter, thread it through
  the generic signature and the `...` forwarding in `.spaNormCore()` →
  `fitSpaNorm()` → `fitSpaNormNB()`, not through every intermediate
  signature.
- Only **single-sample** SpatialExperiment objects are supported
  (`checkSPE` errors on multiple `sample_id`s).
- `counts` must be non-negative integers (`checkSPE`/`checkSeurat`
  enforce this).
- ggplot NSE symbols used in `aes()` are declared in `R/zzz.R` via
  [`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)
  to silence R CMD check “undefined global” NOTEs — add new plotting
  symbols there.
- Example dataset: `data(HumanDLPFC)` — a small filtered Visium
  SpatialExperiment used across man examples, the vignette, and parity
  tests.
