#' @importFrom stats dnbinom qnbinom pnbinom mad median quantile model.matrix
NULL

# session-level cache for the torch device; probing the torch runtime on every
# call is expensive and this is invoked many times per IRLS iteration
.spanorm_env <- new.env(parent = emptyenv())

# clear the cached device (e.g. for tests, or if the device set changes within
# a session)
resetGPUCache <- function() {
  for (key in c("gpu", "device", "gpu.mem.budget")) {
    if (exists(key, envir = .spanorm_env, inherits = FALSE)) {
      rm(list = key, envir = .spanorm_env)
    }
  }
  invisible(NULL)
}

# resolve (once per session) the torch device to use for accelerated ops:
# "cuda", "mps", or "cpu"
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

# per-device dtype: float64 for exact CPU parity on cuda/cpu; MPS cannot
# represent float64, so it uses float32 (its maximum precision)
getBackendDtype <- function() {
  if (getBackendDevice() == "mps") torch::torch_float32() else torch::torch_float64()
}

# bytes per element for the active device's dtype (float32 on MPS, float64
# elsewhere) -- used for sizing memory budgets against genes x cells tensors
.gpuDtypeBytes <- function() {
  if (getBackendDevice() == "mps") 4 else 8
}

checkGPU <- function() {
  # TRUE when an accelerated (non-CPU) torch device is available
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

# determine the accelerator memory budget (bytes) to use for blocked fitting.
# Resolved once per session and cached in the same environment used by
# getBackendDevice(); cleared by resetGPUCache(). A user-supplied
# gpu.mem.budget bypasses detection entirely -- set it to Inf to disable
# blocking outright. Returns Inf when no accelerator is in use.
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
  needed.bytes <- as.numeric(n_genes) * as.numeric(n_cells) * .gpuDtypeBytes() * tensor.multiplier
  nblocks <- if (needed.bytes > budget.bytes) as.integer(ceiling(needed.bytes / budget.bytes)) else 1L
  nblocks <- max(1L, min(as.integer(n_genes), nblocks))
  nblocks
}

# split n_genes into nblocks contiguous index blocks; also called by
# normaliseBlocked() (R/mainSpaNorm.R) so both files carve genes identically.
geneBlockIndices <- function(n_genes, nblocks) {
  split(seq_len(n_genes), ceiling(seq_len(n_genes) / ceiling(n_genes / nblocks)))
}

# helper to detect torch tensors
is_torch_tensor <- function(x) {
  inherits(x, "torch_tensor")
}

# Convert a torch tensor or R object to a base R matrix
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
      n <- dim(work)[[1]]
      I <- torch::torch_eye(as.integer(n), dtype = work$dtype, device = work$device)
      torch::torch_cholesky_solve(I, chol, upper = FALSE)
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

add_vec_mat_gpu <- function(vec, mat, backend = c("auto", "cpu", "gpu")) {
  .vec_mat_op_gpu(vec, mat, op = "add", backend = backend)
}

# Multiply vector and matrix similar to add_vec_mat_gpu
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

# TRUE if `x` (tensor or matrix) holds any NA/NaN/Inf. On a tensor this is a
# scalar reduction copied back, not a full GPU->CPU materialisation of `x`.
hasBadValues <- function(x) {
  if (is_torch_tensor(x)) {
    return(as.logical((torch::torch_isnan(x) | torch::torch_isinf(x))$any()$cpu()))
  }
  any(!is.finite(as.matrix(x)))
}
