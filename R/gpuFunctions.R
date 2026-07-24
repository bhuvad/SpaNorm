#' @importFrom stats dnbinom qnbinom pnbinom mad median quantile model.matrix
NULL

# session-level cache for the torch device; probing the torch runtime on every
# call is expensive and this is invoked many times per IRLS iteration
.spanorm_env <- new.env(parent = emptyenv())

#' Clear the cached GPU device/memory-budget state
#'
#' Clears the session-level cache of the resolved torch device and memory
#' budget (set by \code{\link{getBackendDevice}}/\code{\link{getGPUMemoryBudget}}).
#' Useful in tests, or if the available device/memory changes within a session.
#'
#' @return \code{NULL}, invisibly.
#'
#' @examples
#' resetGPUCache()
#'
#' @export
resetGPUCache <- function() {
  for (key in c("gpu", "device", "gpu.mem.budget")) {
    if (exists(key, envir = .spanorm_env, inherits = FALSE)) {
      rm(list = key, envir = .spanorm_env)
    }
  }
  invisible(NULL)
}

#' Resolve the active torch backend device
#'
#' Resolves (once per session, then cached -- see \code{\link{resetGPUCache}})
#' the torch device to use for accelerated operations: \code{"cuda"} if a CUDA
#' GPU is available, else \code{"mps"} on Apple Silicon, else \code{"cpu"}.
#' Exposed for downstream packages (e.g. spiDE) that push their own tensors
#' onto the same device SpaNorm is using.
#'
#' @return a character, one of \code{"cuda"}, \code{"mps"} or \code{"cpu"}.
#'
#' @examples
#' getBackendDevice()
#'
#' @export
getBackendDevice <- function() {
  if (exists("device", envir = .spanorm_env, inherits = FALSE)) {
    return(get("device", envir = .spanorm_env, inherits = FALSE))
  }

  dev <- "cpu"
  if (requireNamespace("torch", quietly = TRUE)) {
    dev <- tryCatch({
      if (isTRUE(torch::cuda_is_available())) {
        "cuda"
      } else if (isTRUE(torch::backends_mps_is_available())) {
        "mps"
      } else {
        "cpu"
      }
    }, error = function(e) "cpu")
  }

  assign("device", dev, envir = .spanorm_env)
  return(dev)
}

#' Resolve the active torch backend dtype
#'
#' The torch floating-point dtype to use for accelerated tensors: float64
#' everywhere except MPS, which cannot represent float64 and so is capped at
#' float32 (its maximum precision).
#'
#' @return a torch dtype (\code{torch::torch_float64()} or
#'   \code{torch::torch_float32()}).
#'
#' @examples
#' if (checkGPU()) getBackendDtype()
#'
#' @export
getBackendDtype <- function() {
  if (getBackendDevice() == "mps") torch::torch_float32() else torch::torch_float64()
}

#' Bytes per element for the active backend dtype
#'
#' The element width (in bytes) of \code{\link{getBackendDtype}()}'s dtype --
#' 4 for float32 (MPS), 8 for float64 (cuda/cpu). Used for sizing memory
#' budgets against genes x cells tensors; exposed so downstream packages (e.g.
#' spiDE) doing their own GPU memory-budget accounting don't need to duplicate
#' this device/dtype mapping.
#'
#' @return a numeric, 4 or 8.
#'
#' @examples
#' gpuDtypeBytes()
#'
#' @export
gpuDtypeBytes <- function() {
  if (getBackendDevice() == "mps") 4 else 8
}

#' Is an accelerated (GPU) torch device in use?
#'
#' @return a logical, \code{TRUE} iff \code{\link{getBackendDevice}()} resolves
#'   to a non-\code{"cpu"} device.
#'
#' @examples
#' checkGPU()
#'
#' @export
checkGPU <- function() {
  getBackendDevice() != "cpu"
}

# ---- GPU memory budget detection and gene-block sizing -------------------
# Blocking exists only to bound peak accelerator memory for a genes x cells
# fit that would otherwise not fit on-device; it never engages for
# backend = "cpu" (see geneBlockCount()) so the CPU path is untouched.

# Conservative estimate of how many big (genes_block x cells) tensors are
# co-resident on-device within a single block-processing pass: the block's
# counts, fitted log-mean, working residual, dispersion weight, plus several
# temporaries inside dnbinom_gpu()'s log-pmf computation. Deliberately
# generous rather than tight -- this needs empirical validation on real
# hardware; if blocks still OOM or are needlessly small, adjust this value.
GPU_BLOCK_TENSOR_MULTIPLIER <- 10

# fallback budget (bytes) used when neither device-specific detection method
# below succeeds (e.g. nvidia-smi not on PATH, sysctl unavailable)
GPU_MEM_BUDGET_FALLBACK_BYTES <- 2 * 1024^3

# raw (unparsed) free-memory report from nvidia-smi for the active CUDA
# device, or NA_character_ if nvidia-smi is unavailable/errors. Split from
# parsing (.parseNvidiaSmiFreeMem()) so the parser is unit-testable on canned
# text without mocking system2().
.nvidiaSmiFreeMemRaw <- function() {
  tryCatch(
    {
      idx <- torch::cuda_current_device()
      out <- suppressWarnings(system2(
        "nvidia-smi",
        c("--query-gpu=memory.free", "--format=csv,noheader,nounits", paste0("-i=", idx)),
        stdout = TRUE, stderr = FALSE, timeout = 5
      ))
      status <- attr(out, "status")
      if (!is.null(status) && status != 0) {
        return(NA_character_)
      }
      if (length(out) == 0) {
        return(NA_character_)
      }
      out[1]
    },
    error = function(e) NA_character_
  )
}

# parse a single nvidia-smi "--format=csv,noheader,nounits" value (MiB) to
# bytes; NA on anything that doesn't parse as a non-negative number
.parseNvidiaSmiFreeMem <- function(raw) {
  if (is.na(raw) || !nzchar(trimws(raw))) {
    return(NA_real_)
  }
  mib <- suppressWarnings(as.numeric(trimws(raw)))
  if (is.na(mib) || mib < 0) {
    return(NA_real_)
  }
  mib * 1024^2
}

.cudaMemFreeBytes <- function() {
  .parseNvidiaSmiFreeMem(.nvidiaSmiFreeMemRaw())
}

# raw (unparsed) total physical memory report from sysctl (macOS), or
# NA_character_ if unavailable/errors.
.sysctlMemRaw <- function() {
  tryCatch(
    {
      out <- suppressWarnings(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE, timeout = 5))
      if (length(out) == 0) {
        return(NA_character_)
      }
      out[1]
    },
    error = function(e) NA_character_
  )
}

.parseSysctlMem <- function(raw) {
  if (is.na(raw) || !nzchar(trimws(raw))) {
    return(NA_real_)
  }
  bytes <- suppressWarnings(as.numeric(trimws(raw)))
  if (is.na(bytes) || bytes < 0) {
    return(NA_real_)
  }
  bytes
}

# MPS shares unified memory with the OS/other processes and exposes no
# per-process memory-query binding via torch's R package (Metal's
# recommendedMaxWorkingSetSize is not exposed by these bindings); approximate
# the usable budget as a conservative fraction of total system memory.
.mpsMemBudgetBytes <- function(safety.fraction = 0.5) {
  total <- .parseSysctlMem(.sysctlMemRaw())
  if (is.na(total)) {
    return(NA_real_)
  }
  total * safety.fraction
}

#' Determine the accelerator memory budget for blocked fitting
#'
#' Resolves (once per session, then cached -- see \code{\link{resetGPUCache}})
#' the accelerator memory budget (bytes) to use for memory-aware gene-blocked
#' fitting: free CUDA memory via \code{nvidia-smi}, a conservative fraction of
#' total system memory for MPS (unified memory has no per-process query), or a
#' fixed fallback if detection fails. A user-supplied \code{gpu.mem.budget}
#' bypasses detection entirely -- set it to \code{Inf} to disable blocking
#' outright. Returns \code{Inf} when no accelerator is in use. Exposed so
#' downstream packages (e.g. spiDE) doing their own GPU-memory-aware blocking
#' share the same budget/cache as SpaNorm's own fitting.
#'
#' @param gpu.mem.budget \code{NULL} (default, auto-detect) or a single
#'   positive number (bytes).
#' @return a numeric, the memory budget in bytes (\code{Inf} if no accelerator
#'   is in use).
#'
#' @examples
#' getGPUMemoryBudget()
#'
#' @export
getGPUMemoryBudget <- function(gpu.mem.budget = NULL) {
  if (!is.null(gpu.mem.budget)) {
    if (!is.numeric(gpu.mem.budget) || length(gpu.mem.budget) != 1 || is.na(gpu.mem.budget) || gpu.mem.budget <= 0) {
      stop("'gpu.mem.budget' should be a single positive number (bytes), or NULL to auto-detect")
    }
    return(gpu.mem.budget)
  }

  if (exists("gpu.mem.budget", envir = .spanorm_env, inherits = FALSE)) {
    return(get("gpu.mem.budget", envir = .spanorm_env, inherits = FALSE))
  }

  budget <- if (!checkGPU()) {
    Inf # no accelerator in use: blocking is irrelevant, never engaged
  } else {
    dev <- getBackendDevice()
    detected <- switch(dev,
      cuda = .cudaMemFreeBytes(),
      mps = .mpsMemBudgetBytes(),
      NA_real_
    )
    if (is.na(detected)) {
      message(sprintf(
        "Could not auto-detect available %s memory; using a conservative %.1f GiB budget. Override with 'gpu.mem.budget' (bytes).",
        dev, GPU_MEM_BUDGET_FALLBACK_BYTES / 1024^3
      ))
      GPU_MEM_BUDGET_FALLBACK_BYTES
    } else {
      detected
    }
  }

  assign("gpu.mem.budget", budget, envir = .spanorm_env)
  budget
}

# Number of contiguous gene-blocks needed to keep a genes x cells fit within
# `budget.bytes`. Always 1L for backend = "cpu" or when no accelerator is
# available -- this is the single place that guarantees the CPU path is
# never affected by blocking. Mirrors normaliseBlocked()'s
# `ceiling(total/block.size)` clamp (R/mainSpaNorm.R) but is bytes- and
# dtype-aware (float32 on MPS, float64 elsewhere).
#
# `tensor.multiplier` defaults to GPU_BLOCK_TENSOR_MULTIPLIER but is a real
# parameter (not hardcoded) because it estimates a *caller-specific*
# computation graph's peak co-resident tensor count -- this function itself
# stays a generic "bytes needed vs. budget" calculator, agnostic to which
# fitting routine is asking.
geneBlockCount <- function(n_genes, n_cells, n_cov, backend = c("auto", "cpu", "gpu"), budget.bytes = getGPUMemoryBudget(), tensor.multiplier = GPU_BLOCK_TENSOR_MULTIPLIER) {
  backend <- match.arg(backend)
  if (!(backend %in% c("gpu", "auto") && checkGPU())) {
    return(1L)
  }
  if (!is.finite(budget.bytes)) {
    return(1L)
  }

  # n_cov intentionally does not enter the formula: alpha/gmean/psi/W are
  # genes x ncov or ncells x ncov and stay small regardless of gene count --
  # only the genes x cells intermediates drive peak memory.
  needed.bytes <- as.numeric(n_genes) * as.numeric(n_cells) * gpuDtypeBytes() * tensor.multiplier
  nblocks <- if (needed.bytes > budget.bytes) as.integer(ceiling(needed.bytes / budget.bytes)) else 1L
  nblocks <- max(1L, min(as.integer(n_genes), nblocks))
  nblocks
}

# split n_genes into nblocks contiguous index blocks; also called by
# normaliseBlocked() (R/mainSpaNorm.R) so both files carve genes identically.
geneBlockIndices <- function(n_genes, nblocks) {
  split(seq_len(n_genes), ceiling(seq_len(n_genes) / ceiling(n_genes / nblocks)))
}

#' Is an object a torch tensor?
#'
#' @param x an object to test.
#' @return a logical.
#'
#' @examples
#' is_torch_tensor(matrix(1:4, 2, 2))
#'
#' @export
is_torch_tensor <- function(x) {
  inherits(x, "torch_tensor")
}

#' Convert a torch tensor (or R object) to a base R matrix
#'
#' Converts a torch tensor to a base R matrix (moving it to the CPU first if
#' needed); non-tensor input is coerced via \code{\link{as.matrix}}. Exposed
#' so downstream packages (e.g. spiDE) doing their own GPU-blocked computation
#' can convert results back to plain R without depending on torch internals.
#'
#' @param x a torch tensor or an R object coercible via \code{as.matrix}.
#' @return a base R matrix.
#'
#' @examples
#' toRMatrix(matrix(1:4, 2, 2))
#'
#' @export
toRMatrix <- function(x) {
  if (is_torch_tensor(x)) {
    # branch on the tensor's own rank: as.array() drops the dim attribute for
    # 0-D/1-D tensors, so we cannot rely on dim() of the converted R object
    rank <- length(dim(x))
    # move to CPU and materialise as an R array
    arr <- tryCatch(as.array(x$cpu()), error = function(e) NULL)
    if (is.null(arr)) {
      stop("Failed to convert torch tensor to R matrix")
    }
    if (rank <= 1L) {
      # scalar (rank 0) or vector (rank 1) -> column matrix
      return(matrix(as.numeric(arr), ncol = 1L))
    }
    return(as.matrix(arr))
  }

  # Non-tensor input: standard coercion
  return(as.matrix(x))
}

#' Convert a matrix to a torch tensor on the active backend
#'
#' Converts a base R matrix to a torch tensor on the resolved backend device
#' and dtype (\code{\link{getBackendDevice}}/\code{\link{getBackendDtype}}) if
#' an accelerator is available and requested; otherwise returns \code{mat}
#' unchanged (already a tensor, or \code{backend = "cpu"}, or no accelerator
#' present). Exposed so downstream packages (e.g. spiDE) share the exact same
#' device/dtype resolution SpaNorm's own fitting uses.
#'
#' @param mat a base R matrix (or an existing torch tensor, returned as-is).
#' @param ... ignored.
#' @param backend one of \code{"auto"} (default), \code{"cpu"} or \code{"gpu"}.
#' @return a torch tensor, or \code{mat} unchanged.
#'
#' @examples
#' toGPUMatrix(matrix(1:4, 2, 2), backend = "cpu")
#'
#' @export
toGPUMatrix <- function(mat, ..., backend = c("auto", "cpu", "gpu")) {
  backend = match.arg(backend)
  # convert matrix to GPU matrix if an accelerator is available
  if (backend %in% c("gpu", "auto") && checkGPU()) {
    if (is_torch_tensor(mat)) {
      return(mat)
    }
    if (is.matrix(mat)) {
      return(torch::torch_tensor(mat, dtype = getBackendDtype(), device = getBackendDevice()))
    }
  }

  mat
}

#' Convert a vector to a torch tensor on the active backend
#'
#' As \code{\link{toGPUMatrix}}, for a numeric vector; optionally broadcasts a
#' length-1 \code{vec} to length \code{n}.
#'
#' @param vec a numeric vector (or an existing torch tensor, returned as-is).
#' @param n if supplied, the target length; a length-1 \code{vec} is
#'   broadcast to length \code{n}.
#' @param backend one of \code{"auto"} (default), \code{"cpu"} or \code{"gpu"}.
#' @return a torch tensor, or \code{vec} unchanged.
#'
#' @examples
#' toGPUVector(1:4, backend = "cpu")
#'
#' @export
toGPUVector <- function(vec, n = NULL, backend = c("auto", "cpu", "gpu")) {
  backend = match.arg(backend)
  # convert vector to a GPU (torch) vector if an accelerator is available
  if (backend %in% c("gpu", "auto") && checkGPU()) {
    if (is_torch_tensor(vec)) {
      return(vec)
    }
    if (is.null(n) || length(vec) == n) {
      return(torch::torch_tensor(as.numeric(vec), dtype = getBackendDtype(), device = getBackendDevice()))
    } else if (length(vec) == 1) {
      val <- torch::torch_tensor(as.numeric(vec), dtype = getBackendDtype(), device = getBackendDevice())
      ones <- torch::torch_ones(n, dtype = getBackendDtype(), device = getBackendDevice())
      return(torch::torch_multiply(ones, val))
    } else {
      stop("Length of vec does not match n.")
    }
  }

  vec
}

#' Build a diagonal matrix from a vector, GPU-aware
#'
#' As \code{\link[base]{diag}}, but returns a torch tensor (dense, to match
#' the CPU fallback's dense output) when an accelerator is active.
#'
#' @param vec a numeric vector (or torch tensor).
#' @param backend one of \code{"auto"} (default), \code{"cpu"} or \code{"gpu"}.
#' @return a diagonal matrix (base matrix or torch tensor).
#'
#' @examples
#' diag_mat(1:3, backend = "cpu")
#'
#' @export
diag_mat <- function(vec, backend = c("auto", "cpu", "gpu")) {
  backend = match.arg(backend)
  if (backend %in% c("gpu", "auto") && checkGPU()) {
    v <- if (is_torch_tensor(vec)) vec else torch::torch_tensor(as.numeric(vec), dtype = getBackendDtype(), device = getBackendDevice())
    # Ensure v is 1D for torch_diag (handle scalar case)
    if (length(dim(v)) == 0) {
      # Scalar tensor - reshape to 1D with length 1
      v <- v$view(1L)
    }
    mat <- torch::torch_diag(v)
  } else {
    # fallback to base R (dense, to match the GPU path's dense output);
    # nrow set explicitly so a length-1 vector is not misread by diag()
    mat = diag(x = as.numeric(vec), nrow = length(vec))
  }

  return(mat)
}

#' Invert a symmetric positive-definite matrix
#'
#' Inverts a (typically small) symmetric positive-definite matrix using a
#' Cholesky factorisation, falling back to a general solve when the Cholesky
#' fails. When a torch backend is active and \code{mat} is a tensor, the
#' factorisation runs on the accelerator (or CPU for MPS) and the result is
#' returned on the input device. Exposed for downstream packages (e.g. spiDE)
#' that reuse SpaNorm's fitting machinery for Wald-type inference.
#'
#' @param mat a symmetric positive-definite matrix (base matrix or torch tensor).
#' @return the matrix inverse, in the same representation as \code{mat}.
#'
#' @examples
#' m <- crossprod(matrix(rnorm(30), 10, 3))
#' invert_mat(m)
#'
#' @export
invert_mat <- function(mat) {
  # invert matrix using the accelerator if available
  if (is_torch_tensor(mat) && checkGPU()) {
    # MPS cannot run Cholesky/inverse ops; the matrices here are tiny
    # (n_cov x n_cov), so factorise on CPU and move the result back.
    on_mps <- getBackendDevice() == "mps"
    work <- if (on_mps) mat$cpu() else mat

    res <- tryCatch({
      chol <- torch::linalg_cholesky(work)
      torch::torch_cholesky_inverse(chol, upper = FALSE)
    }, error = function(e) {
      torch::torch_inverse(work)
    })

    if (on_mps) {
      res <- res$to(device = mat$device, dtype = mat$dtype)
    }
    return(res)
  }

  # fallback to base R: Cholesky inverse when SPD, general solve otherwise
  # (mirrors the GPU path, which also falls back for non-SPD matrices)
  return(tryCatch(
    Matrix::chol2inv(Matrix::chol(mat)),
    error = function(e) solve(mat)
  ))
}

#' Invert a batch of symmetric positive-definite matrices
#'
#' Batched counterpart of \code{\link{invert_mat}}: inverts a
#' \code{(batch, n, n)} stack of symmetric positive-definite matrices in a
#' single Cholesky factorisation call rather than looping \code{invert_mat()}
#' over each slice. SpaNorm's own model never needs this (its per-gene fit
#' shares a single coefficient across all genes, so it only ever inverts one
#' matrix at a time); it is exposed for downstream packages (e.g. spiDE) that
#' need a genuinely per-gene/per-feature covariance -- a different
#' \code{n x n} matrix for every slice of the batch, varying with a per-gene
#' working weight -- where looping \code{invert_mat()}'s R-level
#' \code{tryCatch}/Cholesky call once per gene is the dominant cost.
#'
#' Unlike \code{invert_mat()}'s per-matrix \code{tryCatch}, a failure here
#' (e.g. one non-SPD slice) falls back to a general solve for the *whole*
#' batch, not just the failing slice -- correctness is unaffected (the
#' fallback solves every slice correctly regardless of conditioning), it is
#' simply not the exact same numerical routine applied per-slice in that
#' (rare) case.
#'
#' @param mat a \code{(batch, n, n)} torch tensor, or a base R array of the
#'   same shape (\code{dim(mat) == c(batch, n, n)}).
#' @return the batch of matrix inverses, in the same representation as
#'   \code{mat}.
#'
#' @examples
#' m <- array(0, c(4, 3, 3))
#' for (i in 1:4) m[i, , ] <- crossprod(matrix(rnorm(9), 3, 3)) + diag(3)
#' invert_mat_batched(m)
#'
#' @export
invert_mat_batched <- function(mat) {
  if (is_torch_tensor(mat) && checkGPU()) {
    # MPS cannot run batched Cholesky/inverse ops either; the matrices here
    # are tiny (n_cov x n_cov, batched over genes), so factorise on CPU.
    on_mps <- getBackendDevice() == "mps"
    work <- if (on_mps) mat$cpu() else mat

    res <- tryCatch({
      chol <- torch::linalg_cholesky(work) # batches natively over dim 1
      torch::torch_cholesky_inverse(chol, upper = FALSE)
    }, error = function(e) {
      torch::torch_inverse(work) # whole-batch fallback; see docs above
    })

    if (on_mps) {
      res <- res$to(device = mat$device, dtype = mat$dtype)
    }
    return(res)
  }

  # CPU fallback: per-slice loop using the same Cholesky-then-solve primitive
  # invert_mat() uses, reassembled into a (batch, n, n) array.
  n <- dim(mat)[2]
  out <- vapply(seq_len(dim(mat)[1]), function(i) {
    m <- mat[i, , ]
    tryCatch(Matrix::chol2inv(Matrix::chol(m)), error = function(e) solve(m))
  }, matrix(0, n, n))
  aperm(out, c(3, 1, 2))
}

#' Tensor-aware \code{tcrossprod}
#'
#' Computes \code{x \%*\% t(y)} (or \code{x \%*\% t(x)} if \code{y} is
#' \code{NULL}) using the accelerator if either argument is a torch tensor and
#' an accelerator is active; falls back to \code{\link[base]{tcrossprod}}
#' otherwise. Exposed so downstream packages (e.g. spiDE) can build their own
#' GPU-blocked linear algebra with the same dual (tensor-or-matrix) semantics
#' SpaNorm's own fitting uses.
#'
#' @param x a matrix or torch tensor.
#' @param y a matrix or torch tensor, or \code{NULL} (default).
#' @return \code{x \%*\% t(y)}, as a matrix or torch tensor matching the input.
#'
#' @examples
#' m <- matrix(rnorm(12), nrow = 3)
#' tcrossprod_gpu(m)
#'
#' @export
tcrossprod_gpu <- function(x, y = NULL) {
  # compute tcrossprod (x %*% t(y)) using the accelerator if available
  if ((is_torch_tensor(x) || (!is.null(y) && is_torch_tensor(y))) && checkGPU()) {
    xt <- if (is_torch_tensor(x)) x else torch::torch_tensor(as.matrix(x), dtype = getBackendDtype(), device = getBackendDevice())
    if (is.null(y)) {
      yt <- xt
    } else {
      yt <- if (is_torch_tensor(y)) y else torch::torch_tensor(as.matrix(y), dtype = getBackendDtype(), device = getBackendDevice())
    }
    res <- torch::torch_matmul(xt, torch::torch_t(yt))
    return(res)
  }

  # fallback to base R
  if (is.null(y)) {
    return(base::tcrossprod(x))
  } else {
    return(base::tcrossprod(x, y))
  }
}

# Common vector-matrix operation helper (GPU-aware with CPU fallback)
.vec_mat_op_gpu <- function(vec, mat, op = c("add", "multiply"), backend = c("auto", "cpu", "gpu")) {
  backend <- match.arg(backend)
  op <- match.arg(op)

  if (backend %in% c("gpu", "auto") && (is_torch_tensor(vec) || is_torch_tensor(mat)) && checkGPU()) {
    v <- if (is_torch_tensor(vec)) vec else torch::torch_tensor(as.numeric(vec), dtype = getBackendDtype(), device = getBackendDevice())
    m <- if (is_torch_tensor(mat)) mat else torch::torch_tensor(as.matrix(mat), dtype = getBackendDtype(), device = getBackendDevice())

    # Determine broadcasting orientation
    mshape <- dim(m)
    m_nr <- tryCatch(as.integer(mshape[[1]]), error = function(e) NULL)
    m_nc <- tryCatch(as.integer(mshape[[2]]), error = function(e) NULL)
    vshape <- dim(v)

    if (length(vshape) == 1) {
      vlen <- tryCatch(as.integer(vshape[[1]]), error = function(e) NULL)
      if (!is.null(m_nr) && !is.null(vlen) && vlen == m_nr) {
        v <- torch::torch_unsqueeze(v, 2L) # (nr, 1)
      } else if (!is.null(m_nc) && !is.null(vlen) && vlen == m_nc) {
        v <- torch::torch_unsqueeze(v, 1L) # (1, nc)
      } else if (!is.null(vlen) && vlen == 1L) {
        # scalar-like
        v <- v$view(1L)
      } else {
        # default to row-wise orientation
        v <- torch::torch_unsqueeze(v, 2L)
      }
    }

    if (op == "add") {
      return(torch::torch_add(m, v))
    } else {
      return(torch::torch_multiply(m, v))
    }
  }

  # CPU fallback: match the GPU broadcast orientation explicitly rather than
  # relying on R's column-major recycling (which only matches when the vector
  # length equals nrow(mat))
  mat <- as.matrix(mat)
  vec <- as.vector(vec)
  fun <- if (op == "add") `+` else `*`
  if (length(vec) == nrow(mat)) {
    return(sweep(mat, 1, vec, FUN = fun)) # per-row broadcast
  } else if (length(vec) == ncol(mat)) {
    return(sweep(mat, 2, vec, FUN = fun)) # per-column broadcast
  } else {
    # scalar or non-conforming length: defer to R recycling
    return(fun(mat, vec))
  }
}

# Cross-product using tcrossprod_gpu: computes t(x) %*% y
crossprod_gpu <- function(x, y = NULL) {
  # compute crossprod using the accelerator if available by reusing tcrossprod_gpu
  if ((is_torch_tensor(x) || (!is.null(y) && is_torch_tensor(y))) && checkGPU()) {
    x_t <- if (is_torch_tensor(x)) torch::torch_t(x) else base::t(as.matrix(x))
    if (is.null(y)) {
      y_t <- x_t
    } else {
      y_t <- if (is_torch_tensor(y)) torch::torch_t(y) else base::t(as.matrix(y))
    }
    return(tcrossprod_gpu(x_t, y_t))
  }

  # fallback to base R
  if (is.null(y)) {
    return(base::crossprod(x))
  } else {
    return(base::crossprod(x, y))
  }
}

#' Broadcast-add a vector to a matrix, GPU-aware
#'
#' Adds \code{vec} to \code{mat}, broadcasting per-row if
#' \code{length(vec) == nrow(mat)} or per-column if
#' \code{length(vec) == ncol(mat)} (matching the GPU path's broadcast
#' orientation explicitly, rather than relying on R's column-major
#' recycling). Runs on the accelerator when either argument is a torch tensor
#' and an accelerator is active.
#'
#' @param vec a numeric vector (or torch tensor).
#' @param mat a matrix (or torch tensor).
#' @param backend one of \code{"auto"} (default), \code{"cpu"} or \code{"gpu"}.
#' @return \code{mat + vec}, broadcast as described above.
#'
#' @examples
#' add_vec_mat_gpu(1:3, matrix(1, 3, 4), backend = "cpu")
#'
#' @export
add_vec_mat_gpu <- function(vec, mat, backend = c("auto", "cpu", "gpu")) {
  .vec_mat_op_gpu(vec, mat, op = "add", backend = backend)
}

#' Broadcast-multiply a vector and a matrix, GPU-aware
#'
#' As \code{\link{add_vec_mat_gpu}}, but multiplying rather than adding.
#'
#' @inheritParams add_vec_mat_gpu
#' @return \code{mat * vec}, broadcast as per \code{\link{add_vec_mat_gpu}}.
#'
#' @examples
#' mult_vec_mat_gpu(1:3, matrix(1, 3, 4), backend = "cpu")
#'
#' @export
mult_vec_mat_gpu <- function(vec, mat, backend = c("auto", "cpu", "gpu")) {
  .vec_mat_op_gpu(vec, mat, op = "multiply", backend = backend)
}

# Helper: convert object to 2D torch tensor with orientation
.to_torch_2d <- function(obj, is_left = TRUE) {
  if (is_torch_tensor(obj)) {
    shape_list <- dim(obj)
    rank <- length(shape_list)
    if (rank == 1) {
      len <- as.integer(shape_list[[1]])
      if (is_left) {
        return(obj$view(c(1L, len)))
      } else {
        return(obj$view(c(len, 1L)))
      }
    }
    return(obj)
  }

  if (is.matrix(obj)) {
    return(torch::torch_tensor(as.matrix(obj), dtype = getBackendDtype(), device = getBackendDevice()))
  }

  v <- as.numeric(obj)
  t <- torch::torch_tensor(v, dtype = getBackendDtype(), device = getBackendDevice())
  if (is_left) {
    return(t$view(c(1L, as.integer(length(v)))))
  } else {
    return(t$view(c(as.integer(length(v)), 1L)))
  }
}

#' Tensor-aware matrix multiplication
#'
#' Computes \code{x \%*\% y} using the accelerator if either argument is a
#' torch tensor and an accelerator is active (vectors are reshaped so the
#' result matches base R's \code{\%*\%} shape semantics); falls back to base
#' \code{\%*\%} otherwise. Named explicitly (rather than overriding
#' \code{\%*\%}) so base matrix multiplication semantics are preserved
#' everywhere else. Exposed so downstream packages (e.g. spiDE) can build
#' batched linear algebra (e.g. a shared design matrix against many
#' per-gene/per-block weight vectors) without depending on torch internals.
#'
#' @param x,y a matrix, numeric vector, or torch tensor.
#' @return \code{x \%*\% y}, as a matrix or torch tensor matching the input.
#'
#' @examples
#' matmul_gpu(matrix(1:6, 2, 3), matrix(1:6, 3, 2))
#'
#' @export
matmul_gpu <- function(x, y) {
  # Tensor-aware matrix multiplication. Falls back to base %*%.
  # Named explicitly (rather than overriding `%*%`) so base matrix
  # multiplication semantics are preserved everywhere else in the package.
  if ((is_torch_tensor(x) || is_torch_tensor(y)) && checkGPU()) {
    ax <- .to_torch_2d(x, is_left = TRUE)
    by <- .to_torch_2d(y, is_left = FALSE)
    res <- torch::torch_matmul(ax, by)

    # Mimic R shape behavior when vectors involved
    x_is_vec <- (!is.matrix(x)) && !is_torch_tensor(x) && is.atomic(x) && is.null(dim(x))
    y_is_vec <- (!is.matrix(y)) && !is_torch_tensor(y) && is.atomic(y) && is.null(dim(y))
    if (is_torch_tensor(x)) {
      x_is_vec <- tryCatch(length(dim(x)) == 1, error = function(e) FALSE)
    }
    if (is_torch_tensor(y)) {
      y_is_vec <- tryCatch(length(dim(y)) == 1, error = function(e) FALSE)
    }

    if (x_is_vec && y_is_vec) {
      # return scalar (inner product)
      res <- torch::torch_squeeze(res)
    } else if (x_is_vec && !y_is_vec) {
      # row vector result (1, ncol) -> squeeze dim 1
      res <- torch::torch_squeeze(res, dim = 1L)
    } else if (!x_is_vec && y_is_vec) {
      # column vector result (nrow, 1) -> squeeze dim 2
      res <- torch::torch_squeeze(res, dim = 2L)
    }
    return(res)
  }

  # CPU/base fallback
  base::`%*%`(x, y)
}

rowMeans_gpu <- function(mat) {
  # calculate row means using the accelerator if available
  if (is_torch_tensor(mat) && checkGPU()) {
    res <- torch::torch_mean(mat, dim = 2L)
    return(res)
  }

  # fallback to base R
  return(matrixStats::rowMeans2(mat))
}

colMeans_gpu <- function(mat) {
  # calculate column means using the accelerator if available
  if (is_torch_tensor(mat) && checkGPU()) {
    res <- torch::torch_mean(mat, dim = 1L)
    return(res)
  }

  # fallback to base R
  return(matrixStats::colMeans2(mat))
}

#' Row sums, GPU-aware
#'
#' As \code{\link[matrixStats]{rowSums2}}, but runs on the accelerator when
#' \code{mat} is a torch tensor and an accelerator is active.
#'
#' @param mat a matrix or torch tensor.
#' @return a numeric vector (or torch tensor) of row sums.
#'
#' @examples
#' rowSums_gpu(matrix(1:6, 2, 3))
#'
#' @export
rowSums_gpu <- function(mat) {
  # calculate row sums using the accelerator if available
  if (is_torch_tensor(mat) && checkGPU()) {
    res <- torch::torch_sum(mat, dim = 2L)
    return(res)
  }

  # fallback to base R
  return(matrixStats::rowSums2(mat))
}

colSums_gpu <- function(mat) {
  # calculate column sums using the accelerator if available
  if (is_torch_tensor(mat) && checkGPU()) {
    res <- torch::torch_sum(mat, dim = 1L)
    return(res)
  }

  # fallback to base R
  return(matrixStats::colSums2(mat))
}

# --- Negative-binomial (mean/size parameterisation) on the torch backend ------
# stats::dnbinom/pnbinom/qnbinom use size = 1/dispersion and mean mu. `size` is
# typically a per-gene (per-row) vector broadcast across cells (columns), matching
# how stats:: recycles it (column-major) against the genes x cells matrices.
#
# The torch paths target the domain SpaNorm actually uses and are validated
# against stats:: there: non-negative integer support, mu > 0, and (for the
# quantile) p in [0, 1). Outside that domain they are not drop-in replacements
# for stats:: -- e.g. mu = 0 gives NaN rather than a point mass at 0, negative q
# is treated as 0, and p >= 1 (whose true quantile is +Inf) is bounded by `cap`
# in qnbinom_gpu rather than returning Inf.

# coerce a matrix/scalar argument to a torch tensor on ref's device and dtype
.nb_to_tensor <- function(m, ref) {
  if (is_torch_tensor(m)) return(m$to(device = ref$device, dtype = ref$dtype))
  torch::torch_tensor(as.matrix(m), dtype = ref$dtype, device = ref$device)
}

# coerce `size` to a tensor; a per-gene vector (length == nrow) becomes (nrow, 1)
# so it broadcasts row-wise, exactly as stats:: recycles it
.nb_size_tensor <- function(size, ref) {
  s <- if (is_torch_tensor(size)) size$to(device = ref$device, dtype = ref$dtype)
       else torch::torch_tensor(as.numeric(size), dtype = ref$dtype, device = ref$device)
  nr <- dim(ref)[[1]]
  if (length(dim(s)) == 1 && dim(s)[[1]] == nr) s <- s$view(c(nr, 1L))
  s
}

# NB log-pmf at tensor `x`: lgamma(x+r) - lgamma(r) - lgamma(x+1)
#   + r*log(r/(r+mu)) + x*log(mu/(r+mu))   (r = size)
.nb_logpmf_torch <- function(x, mu, sz) {
  torch::torch_lgamma(x + sz) - torch::torch_lgamma(sz) - torch::torch_lgamma(x + 1) +
    sz * (torch::torch_log(sz) - torch::torch_log(sz + mu)) +
    x  * (torch::torch_log(mu) - torch::torch_log(sz + mu))
}

#' Negative binomial density (mean/size parameterisation), GPU-aware
#'
#' As \code{\link[stats]{dnbinom}} (with \code{mu} instead of \code{prob}),
#' but runs the log-pmf on the accelerator when \code{x} or \code{mu} is a
#' torch tensor and an accelerator is active. Targets the domain SpaNorm
#' actually uses (non-negative integer \code{x}, \code{mu > 0}) -- outside
#' that domain this is not a drop-in replacement for \code{stats::dnbinom}
#' (e.g. \code{mu = 0} gives \code{NaN} rather than a point mass at 0).
#' Exposed so downstream packages (e.g. spiDE) computing their own per-gene NB
#' log-likelihoods reuse the same GPU/CPU-dual-path NB math SpaNorm's own
#' fitting uses, rather than duplicating it.
#'
#' @param x quantiles (matrix or torch tensor).
#' @param mu the mean (matrix or torch tensor, same shape as \code{x}).
#' @param size the NB size parameter (\code{1/dispersion}); a length-\code{nrow(x)}
#'   vector broadcasts row-wise, matching how \code{stats::dnbinom} recycles it.
#' @param log a logical, return the log-density (default \code{FALSE}).
#' @return the (log-)density, as a matrix or torch tensor matching the input.
#'
#' @examples
#' dnbinom_gpu(matrix(0:3, 2, 2), mu = matrix(1, 2, 2), size = c(2, 2))
#'
#' @export
dnbinom_gpu <- function(x, mu, size, log = FALSE) {
  if ((is_torch_tensor(x) || is_torch_tensor(mu)) && checkGPU()) {
    ref  <- if (is_torch_tensor(mu)) mu else x
    mu.t <- .nb_to_tensor(mu, ref)
    x.t  <- .nb_to_tensor(x, ref)
    sz.t <- .nb_size_tensor(size, ref)
    lp   <- .nb_logpmf_torch(x.t, mu.t, sz.t)
    return(if (log) lp else torch::torch_exp(lp))
  }
  # CPU fallback (as.numeric, not as.vector, so a tensor `size` coerces correctly)
  dnbinom(as.matrix(x), mu = as.matrix(mu), size = as.numeric(size), log = log)
}

pnbinom_gpu <- function(q, mu, size, log = FALSE) {
  if ((is_torch_tensor(q) || is_torch_tensor(mu)) && checkGPU()) {
    # discrete CDF: cumulative sum of the pmf up to floor(q).
    ref  <- if (is_torch_tensor(mu)) mu else q
    mu.t <- .nb_to_tensor(mu, ref)
    sz.t <- .nb_size_tensor(size, ref)
    # broadcast floor(q) to mu's shape for the per-element comparison
    qf   <- torch::torch_clamp(torch::torch_floor(.nb_to_tensor(q, ref)), min = 0) +
            torch::torch_zeros_like(mu.t)
    K    <- as.integer(as.numeric(qf$max()$cpu()))
    res  <- torch::torch_zeros_like(mu.t)
    run  <- torch::torch_zeros_like(mu.t)
    for (k in 0:K) {
      run <- run + torch::torch_exp(.nb_logpmf_torch(torch::torch_full_like(mu.t, k), mu.t, sz.t))
      res <- torch::torch_where(qf == k, run, res)
    }
    return(res)
  }
  pnbinom(as.matrix(q), mu = as.matrix(mu), size = as.numeric(size))
}

qnbinom_gpu <- function(p, mu, size, log = FALSE) {
  if ((is_torch_tensor(p) || is_torch_tensor(mu)) && checkGPU()) {
    # discrete quantile: smallest integer k with CDF(k) >= p, found by
    # accumulating the pmf until every element has crossed its probability.
    ref  <- if (is_torch_tensor(mu)) mu else p
    mu.t <- .nb_to_tensor(mu, ref)
    sz.t <- .nb_size_tensor(size, ref)
    p.t  <- .nb_to_tensor(p, ref) + torch::torch_zeros_like(mu.t) # broadcast to mu shape
    res  <- torch::torch_zeros_like(mu.t)
    run  <- torch::torch_zeros_like(mu.t)
    done <- torch::torch_zeros_like(mu.t)$to(dtype = torch::torch_bool())
    k    <- 0L
    cap  <- 1000000L # guard against p == 1 (quantile is +Inf)
    repeat {
      run   <- run + torch::torch_exp(.nb_logpmf_torch(torch::torch_full_like(mu.t, k), mu.t, sz.t))
      newly <- (run >= p.t) & done$logical_not()
      res   <- torch::torch_where(newly, torch::torch_full_like(mu.t, k), res)
      done  <- done | newly
      if (as.logical(done$all()$cpu()) || k >= cap) break
      k <- k + 1L
    }
    return(res)
  }
  qnbinom(as.matrix(p), mu = as.matrix(mu), size = as.numeric(size))
}

# Winsorise each column of `alpha` to median +/- k*MAD (statistics taken per
# column, over rows). On the accelerator this stays on-device -- previously the
# IRLS loop copied `alpha` GPU->CPU->GPU every inner iteration just to run
# matrixStats here. The CPU branch reproduces that exact behaviour (matrixStats
# centre = colMedians, constant = 1.4826; torch's linearly-interpolated 0.5
# quantile matches R's type-7 median, so the two paths agree).
winsoriseCols <- function(alpha, k = DEFAULT_WINSOR) {
  .checkWinsor(k)
  if (is.infinite(k)) return(alpha) # no winsorisation (avoids k*mad == NaN when mad == 0)

  if (is_torch_tensor(alpha) && checkGPU()) {
    ncov <- as.integer(dim(alpha)[[2]])
    med <- torch::torch_quantile(alpha, 0.5, dim = 1L) # per-column median
    mad <- torch::torch_quantile(
      torch::torch_abs(alpha - med$view(c(1L, ncov))), 0.5, dim = 1L
    ) * 1.4826
    a.max <- (med + k * mad)$view(c(1L, ncov))
    a.min <- (med - k * mad)$view(c(1L, ncov))
    return(torch::torch_maximum(torch::torch_minimum(alpha, a.max), a.min))
  }

  # CPU: exact previous behaviour (matrixStats + transpose recycling)
  alpha <- as.matrix(alpha)
  a.med <- matrixStats::colMedians(alpha)
  a.mad <- matrixStats::colMads(alpha)
  a.max <- a.med + k * a.mad
  a.min <- a.med - k * a.mad
  alpha <- pmin(t(alpha), a.max) # clamp +ve outliers
  alpha <- t(pmax(alpha, a.min)) # clamp -ve outliers
  return(alpha)
}

# One-sided, row-wise GPU winsorisation of log-means: clamps each row (gene)
# to at most its own median + k*MAD (mirrors calculateMu()'s CPU-only
# `pmin(mu, rowMedians(mu) + k*rowMads(mu))` exactly). Unlike winsoriseCols()
# above this is (a) one-sided -- only the upper tail is clamped, matching
# calculateMu()'s intent of bounding exploding fitted means, not both tails of
# a coefficient -- and (b) row- rather than column-wise, since here rows are
# genes and columns are cells.
.winsoriseRowsGPU <- function(lmu, k) {
  if (is_torch_tensor(lmu) && checkGPU()) {
    nr <- as.integer(dim(lmu)[[1]])
    med <- torch::torch_quantile(lmu, 0.5, dim = 2L) # per-row median
    mad <- torch::torch_quantile(
      torch::torch_abs(lmu - med$view(c(nr, 1L))), 0.5, dim = 2L
    ) * 1.4826
    lmu.max <- (med + k * mad)$view(c(nr, 1L))
    return(torch::torch_minimum(lmu, lmu.max))
  }

  lmu <- as.matrix(lmu)
  lmu.max <- matrixStats::rowMedians(lmu) + k * matrixStats::rowMads(lmu)
  pmin(lmu, lmu.max)
}

#' Does an object hold any NA/NaN/Inf?
#'
#' \code{TRUE} if \code{x} (tensor or matrix) holds any \code{NA}/\code{NaN}/
#' \code{Inf}. On a tensor this is a scalar reduction copied back, not a full
#' GPU->CPU materialisation of \code{x}.
#'
#' @param x a matrix or torch tensor.
#' @return a logical.
#'
#' @examples
#' hasBadValues(matrix(c(1, NA, 3, 4), 2, 2))
#'
#' @export
hasBadValues <- function(x) {
  if (is_torch_tensor(x)) {
    return(as.logical((torch::torch_isnan(x) | torch::torch_isinf(x))$any()$cpu()))
  }
  any(!is.finite(as.matrix(x)))
}
