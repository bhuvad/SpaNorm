# Package index

## All functions

- [`HumanDLPFC`](https://bhuvad.github.io/spaNorm/reference/HumanDLPFC.md)
  : Human dorsolateral prefrontal cortex (DLPFC) visium sample

- [`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md) :
  Spatially-dependent normalisation for spatial transcriptomics data

- [`` `$`( ``*`<SpaNormFit>`*`)`](https://bhuvad.github.io/spaNorm/reference/SpaNormFit.md)
  : An S4 class to store a SpaNorm model fit

- [`SpaNormPCA()`](https://bhuvad.github.io/spaNorm/reference/SpaNormPCA.md)
  : GLM-based (SpaNorm) PCA

- [`SpaNormSVG()`](https://bhuvad.github.io/spaNorm/reference/SpaNormSVG.md)
  : Model-based spatially variable gene (SVG) calling

- [`add_vec_mat_gpu()`](https://bhuvad.github.io/spaNorm/reference/add_vec_mat_gpu.md)
  : Broadcast-add a vector to a matrix, GPU-aware

- [`calculateMu()`](https://bhuvad.github.io/spaNorm/reference/calculateMu.md)
  : Compute fitted means from a negative binomial GLM fit

- [`checkGPU()`](https://bhuvad.github.io/spaNorm/reference/checkGPU.md)
  : Is an accelerated (GPU) torch device in use?

- [`diag_mat()`](https://bhuvad.github.io/spaNorm/reference/diag_mat.md)
  : Build a diagonal matrix from a vector, GPU-aware

- [`dnbinom_gpu()`](https://bhuvad.github.io/spaNorm/reference/dnbinom_gpu.md)
  : Negative binomial density (mean/size parameterisation), GPU-aware

- [`fastSizeFactors()`](https://bhuvad.github.io/spaNorm/reference/fastSizeFactors.md)
  : Filter genes based on expression

- [`filterGenes()`](https://bhuvad.github.io/spaNorm/reference/filterGenes.md)
  : Filter genes based on expression

- [`fitNB()`](https://bhuvad.github.io/spaNorm/reference/fitNB.md) : Fit
  a per-gene negative binomial GLM

- [`getBackendDevice()`](https://bhuvad.github.io/spaNorm/reference/getBackendDevice.md)
  : Resolve the active torch backend device

- [`getBackendDtype()`](https://bhuvad.github.io/spaNorm/reference/getBackendDtype.md)
  : Resolve the active torch backend dtype

- [`getGPUMemoryBudget()`](https://bhuvad.github.io/spaNorm/reference/getGPUMemoryBudget.md)
  : Determine the accelerator memory budget for blocked fitting

- [`gpuDtypeBytes()`](https://bhuvad.github.io/spaNorm/reference/gpuDtypeBytes.md)
  : Bytes per element for the active backend dtype

- [`hasBadValues()`](https://bhuvad.github.io/spaNorm/reference/hasBadValues.md)
  : Does an object hold any NA/NaN/Inf?

- [`invert_mat()`](https://bhuvad.github.io/spaNorm/reference/invert_mat.md)
  : Invert a symmetric positive-definite matrix

- [`invert_mat_batched()`](https://bhuvad.github.io/spaNorm/reference/invert_mat_batched.md)
  : Invert a batch of symmetric positive-definite matrices

- [`is_torch_tensor()`](https://bhuvad.github.io/spaNorm/reference/is_torch_tensor.md)
  : Is an object a torch tensor?

- [`matmul_gpu()`](https://bhuvad.github.io/spaNorm/reference/matmul_gpu.md)
  : Tensor-aware matrix multiplication

- [`mult_vec_mat_gpu()`](https://bhuvad.github.io/spaNorm/reference/mult_vec_mat_gpu.md)
  : Broadcast-multiply a vector and a matrix, GPU-aware

- [`plotCovariate()`](https://bhuvad.github.io/spaNorm/reference/plotCovariate.md)
  : Diagnostic plot of predicted expression for a covariate

- [`plotSpatial()`](https://bhuvad.github.io/spaNorm/reference/plotSpatial.md)
  : Plot spatial transcriptomic annotations per spot

- [`resetGPUCache()`](https://bhuvad.github.io/spaNorm/reference/resetGPUCache.md)
  : Clear the cached GPU device/memory-budget state

- [`rowSums_gpu()`](https://bhuvad.github.io/spaNorm/reference/rowSums_gpu.md)
  : Row sums, GPU-aware

- [`tcrossprod_gpu()`](https://bhuvad.github.io/spaNorm/reference/tcrossprod_gpu.md)
  :

  Tensor-aware `tcrossprod`

- [`toGPUMatrix()`](https://bhuvad.github.io/spaNorm/reference/toGPUMatrix.md)
  : Convert a matrix to a torch tensor on the active backend

- [`toGPUVector()`](https://bhuvad.github.io/spaNorm/reference/toGPUVector.md)
  : Convert a vector to a torch tensor on the active backend

- [`toRMatrix()`](https://bhuvad.github.io/spaNorm/reference/toRMatrix.md)
  : Convert a torch tensor (or R object) to a base R matrix

- [`topSVGs()`](https://bhuvad.github.io/spaNorm/reference/topSVGs.md) :
  Export top SVG results to a data frame
