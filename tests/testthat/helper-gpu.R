# The torch R package can be installed while its libtorch/lantern backend is
# not: CI installs Suggests but never runs torch::install_torch(), and every
# tensor call then errors with "Lantern is not loaded". requireNamespace() (and
# so skip_if_not_installed("torch")) cannot see this, so probe by actually
# building a tensor. Resolved once per run.
torchIsUsable <- local({
  usable <- NULL
  function() {
    if (is.null(usable)) {
      usable <<- requireNamespace("torch", quietly = TRUE) &&
        isTRUE(tryCatch({
          torch::torch_tensor(1, dtype = torch::torch_float64())
          TRUE
        }, error = function(e) FALSE))
    }
    usable
  }
})

skip_if_no_torch <- function() {
  testthat::skip_if_not(torchIsUsable(), "torch backend (libtorch/lantern) not installed")
}

skip_if_no_gpu <- function() {
  testthat::skip_if_not(checkGPU(), "torch GPU/MPS not available")
}

# accelerator tolerance: near machine precision on float64 devices (cuda/cpu),
# looser on MPS which can only run float32
gpu_tol <- function() {
  if (getBackendDevice() == "mps") 1e-6 else 1e-10
}

# a memory budget (bytes) small enough to force geneBlockCount() to split
# an n_genes x n_cells fit into multiple blocks, for exercising the blocked
# GPU fitting path in tests
tinyGpuBudget <- function(n_genes, n_cells, divisor = 6) {
  (n_genes * n_cells * gpuDtypeBytes() * GPU_BLOCK_TENSOR_MULTIPLIER) / divisor
}
