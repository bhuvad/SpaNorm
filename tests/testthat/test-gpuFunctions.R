# skip_if_no_gpu()/gpu_tol()/tinyGpuBudget() live in helper-gpu.R (auto-sourced
# by testthat before any test file runs)

as_numeric_matrix <- function(x) {
  if (is_torch_tensor(x)) {
    return(toRMatrix(x))
  }
  as.matrix(x)
}

as_numeric_vector <- function(x) {
  res <- as_numeric_matrix(x)
  as.numeric(res)
}

test_that("GPU matrix and vector conversions round-trip", {
  skip_if_no_gpu()

  set.seed(1)
  m <- matrix(rnorm(9), nrow = 3)
  v <- rnorm(5)

  mgpu <- toGPUMatrix(m, backend = "gpu")
  vgpu <- toGPUVector(v, backend = "gpu")

  expect_true(is_torch_tensor(mgpu))
  expect_true(is_torch_tensor(vgpu))

  expect_equal(toRMatrix(mgpu), m, tolerance = gpu_tol())
  expect_equal(as_numeric_vector(vgpu), as.numeric(v), tolerance = gpu_tol())
})

test_that("GPU tcrossprod and crossprod match CPU", {
  skip_if_no_gpu()

  set.seed(2)
  m <- matrix(rnorm(12), nrow = 3)

  tc_gpu <- toRMatrix(tcrossprod_gpu(toGPUMatrix(m, backend = "gpu")))
  cp_gpu <- toRMatrix(crossprod_gpu(toGPUMatrix(m, backend = "gpu")))

  expect_equal(tc_gpu, tcrossprod(m), tolerance = gpu_tol())
  expect_equal(cp_gpu, crossprod(m), tolerance = gpu_tol())
})

test_that("Vector-matrix add and multiply broadcast match CPU", {
  skip_if_no_gpu()

  set.seed(3)
  m <- matrix(rnorm(15), nrow = 5)
  v_row <- rnorm(5) # length matches nrow
  v_col <- rnorm(3) # length matches ncol

  add_gpu <- toRMatrix(add_vec_mat_gpu(v_row, toGPUMatrix(m, backend = "gpu")))
  mult_gpu <- toRMatrix(mult_vec_mat_gpu(v_col, toGPUMatrix(m, backend = "gpu")))

  expect_equal(add_gpu, m + v_row, tolerance = gpu_tol())
  expect_equal(mult_gpu, sweep(m, 2, v_col, `*`), tolerance = gpu_tol())
})

test_that("diag_mat and invert_mat match CPU", {
  skip_if_no_gpu()

  set.seed(4)
  v <- rnorm(4)
  base_diag <- diag(v)
  gpu_diag <- toRMatrix(diag_mat(v, backend = "gpu"))
  expect_equal(gpu_diag, base_diag, tolerance = gpu_tol())

  # Test length-1 vector (edge case)
  v1 <- 2.5
  base_diag1 <- diag(v1, nrow = 1, ncol = 1)
  gpu_diag1 <- toRMatrix(diag_mat(v1, backend = "gpu"))
  expect_equal(gpu_diag1, base_diag1, tolerance = gpu_tol())

  a <- matrix(rnorm(16), nrow = 4)
  pd <- crossprod(a) + diag(4) * 0.1
  inv_gpu <- toRMatrix(invert_mat(toGPUMatrix(pd, backend = "gpu")))
  expect_equal(inv_gpu, solve(pd), tolerance = gpu_tol())
})

test_that("Row/col summaries match CPU", {
  skip_if_no_gpu()

  set.seed(5)
  m <- matrix(rnorm(18), nrow = 6)
  mgpu <- toGPUMatrix(m, backend = "gpu")

  expect_equal(as_numeric_vector(rowMeans_gpu(mgpu)), rowMeans(m), tolerance = gpu_tol())
  expect_equal(as_numeric_vector(colMeans_gpu(mgpu)), colMeans(m), tolerance = gpu_tol())
  expect_equal(as_numeric_vector(rowSums_gpu(mgpu)), rowSums(m), tolerance = gpu_tol())
  expect_equal(as_numeric_vector(colSums_gpu(mgpu)), colSums(m), tolerance = gpu_tol())
})

test_that("Tensor-aware matmul_gpu matches CPU", {
  skip_if_no_gpu()

  x <- 1:3
  y <- 4:6
  xy_gpu <- as_numeric_vector(matmul_gpu(toGPUVector(x, backend = "gpu"), toGPUVector(y, backend = "gpu")))
  expect_equal(xy_gpu, sum(x * y), tolerance = gpu_tol())

  mat <- matrix(1:6, nrow = 2)
  vec <- c(1, 2, 3)
  mat_gpu <- toGPUMatrix(mat, backend = "gpu")
  vec_gpu <- toGPUVector(vec, backend = "gpu")
  mv_gpu <- as_numeric_vector(matmul_gpu(mat_gpu, vec_gpu))
  expect_equal(mv_gpu, as.numeric(mat %*% vec), tolerance = gpu_tol())
})

test_that("GPU negative-binomial d/p/q match stats::", {
  skip_if_no_gpu()

  set.seed(20)
  G <- 6; N <- 5
  mu   <- matrix(runif(G * N, 0.5, 20), G, N)
  size <- runif(G, 0.4, 5)                 # per-gene (per-row)
  x    <- matrix(rpois(G * N, 8), G, N)     # counts
  p    <- matrix(runif(G * N, 0.02, 0.98), G, N)

  # per-gene stats:: reference (recycles size across columns)
  ref <- function(fun, val) {
    out <- matrix(0, G, N)
    for (g in seq_len(G)) out[g, ] <- fun(val[g, ], size = size[g], mu = mu[g, ])
    out
  }

  mu_gpu <- toGPUMatrix(mu, backend = "gpu")
  x_gpu  <- toGPUMatrix(x, backend = "gpu")
  p_gpu  <- toGPUMatrix(p, backend = "gpu")

  # density (log) — the hot-loop path
  d_gpu <- toRMatrix(dnbinom_gpu(x_gpu, mu = mu_gpu, size = size, log = TRUE))
  expect_equal(d_gpu, ref(function(v, ...) dnbinom(v, ..., log = TRUE), x), tolerance = gpu_tol() * 1e3)

  # CDF
  p_cdf <- toRMatrix(pnbinom_gpu(x_gpu, mu = mu_gpu, size = size))
  expect_equal(p_cdf, ref(pnbinom, x), tolerance = gpu_tol() * 1e3)

  # quantile (integer-valued; must match exactly)
  q_gpu <- toRMatrix(qnbinom_gpu(p_gpu, mu = mu_gpu, size = size))
  expect_equal(q_gpu, ref(qnbinom, p), tolerance = 0)

  # scalar p (as used by normaliseMedianBio: qnbinom_gpu(0.5, ...))
  q_med <- toRMatrix(qnbinom_gpu(0.5, mu = mu_gpu, size = size))
  ref_med <- matrix(0, G, N)
  for (g in seq_len(G)) ref_med[g, ] <- qnbinom(0.5, size = size[g], mu = mu[g, ])
  expect_equal(q_med, ref_med, tolerance = 0)
})

# ---- CPU fallbacks (run without a GPU) ----
# These guard the code paths that execute for every user who does not have a
# GPU. They must NOT skip_if_no_gpu().

test_that("CPU fallbacks match base R", {
  set.seed(10)
  m <- matrix(rnorm(15), nrow = 5) # 5 x 3
  v_row <- rnorm(5) # length == nrow(m)
  v_col <- rnorm(3) # length == ncol(m)

  # cross products and matrix multiply
  expect_equal(as.matrix(tcrossprod_gpu(m)), tcrossprod(m), tolerance = 1e-8)
  expect_equal(as.matrix(crossprod_gpu(m)), crossprod(m), tolerance = 1e-8)
  b <- matrix(rnorm(6), nrow = 3) # 3 x 2
  expect_equal(as.matrix(matmul_gpu(m, b)), m %*% b, tolerance = 1e-8)

  # vector-matrix broadcast must match sweep() for BOTH orientations
  expect_equal(add_vec_mat_gpu(v_row, m, backend = "cpu"), sweep(m, 1, v_row, `+`), tolerance = 1e-8)
  expect_equal(add_vec_mat_gpu(v_col, m, backend = "cpu"), sweep(m, 2, v_col, `+`), tolerance = 1e-8)
  expect_equal(mult_vec_mat_gpu(v_row, m, backend = "cpu"), sweep(m, 1, v_row, `*`), tolerance = 1e-8)
  expect_equal(mult_vec_mat_gpu(v_col, m, backend = "cpu"), sweep(m, 2, v_col, `*`), tolerance = 1e-8)

  # diagonal matrix
  expect_equal(as.matrix(diag_mat(v_col, backend = "cpu")), diag(v_col), tolerance = 1e-8)

  # row/col summaries
  expect_equal(as.numeric(rowMeans_gpu(m)), rowMeans(m), tolerance = 1e-8)
  expect_equal(as.numeric(colMeans_gpu(m)), colMeans(m), tolerance = 1e-8)
  expect_equal(as.numeric(rowSums_gpu(m)), rowSums(m), tolerance = 1e-8)
  expect_equal(as.numeric(colSums_gpu(m)), colSums(m), tolerance = 1e-8)
})

test_that("negative-binomial d/p/q CPU path matches stats::", {
  # plain R-matrix inputs always take the stats:: path (no tensor operands)
  set.seed(12)
  G <- 4; N <- 5
  mu   <- matrix(runif(G * N, 0.5, 15), G, N)
  size <- runif(G, 0.5, 4)
  x    <- matrix(rpois(G * N, 6), G, N)
  p    <- matrix(runif(G * N, 0.02, 0.98), G, N)
  ref <- function(fun, val, ...) {
    out <- matrix(0, G, N)
    for (g in seq_len(G)) out[g, ] <- fun(val[g, ], size = size[g], mu = mu[g, ], ...)
    out
  }
  expect_equal(dnbinom_gpu(x, mu = mu, size = size, log = TRUE),
               ref(dnbinom, x, log = TRUE), tolerance = 1e-8)
  expect_equal(pnbinom_gpu(x, mu = mu, size = size), ref(pnbinom, x), tolerance = 1e-8)
  expect_equal(qnbinom_gpu(p, mu = mu, size = size), ref(qnbinom, p), tolerance = 0)
})

test_that("winsoriseCols matches the exact matrixStats winsorisation", {
  set.seed(31)
  a <- matrix(rnorm(50 * 4), 50, 4)
  a[1, 1] <- 40; a[2, 3] <- -40 # column outliers

  # reference == the previous in-loop code (matrixStats + transpose recycling)
  med <- matrixStats::colMedians(a)
  md <- matrixStats::colMads(a)
  ref <- t(pmax(pmin(t(a), med + 4 * md), med - 4 * md))

  # CPU path must be bit-identical
  expect_equal(winsoriseCols(a, 4), ref, tolerance = 1e-12)

  skip_if_no_gpu()
  g <- toRMatrix(winsoriseCols(toGPUMatrix(a, backend = "gpu"), 4))
  expect_equal(g, ref, tolerance = gpu_tol())
})

test_that("winsoriseCols validates k and handles k = Inf without NaN", {
  # a column with zero MAD (all values identical) is a realistic case: e.g.
  # several genes converging to the same fitted coefficient
  m <- cbind(rep(5, 6), 1:6)

  # k <= 0 must error, not silently collapse every gene to one value
  expect_error(winsoriseCols(m, k = 0), "positive")
  expect_error(winsoriseCols(m, k = -1), "positive")

  # k = Inf means "no winsorisation": must not compute Inf * mad == NaN when
  # mad == 0, and must leave the input unchanged
  res <- winsoriseCols(m, k = Inf)
  expect_false(anyNA(res))
  expect_equal(res, m)
})

test_that("winsorisePsi and calculateMu validate winsor and handle Inf", {
  psi <- c(1, 2, 3)
  expect_error(winsorisePsi(psi, winsor = 0), "positive")
  expect_error(winsorisePsi(psi, winsor = -2), "positive")
  expect_equal(winsorisePsi(psi, winsor = Inf), psi)

  gmean <- 0
  alpha <- matrix(c(1, 1, 2, 2), 2, 2) # row 1 and row 2 identical -> mad == 0
  W <- diag(2)
  expect_error(calculateMu(gmean, alpha, W, winsor = 0), "positive")
  res <- calculateMu(gmean, alpha, W, winsor = Inf)
  expect_false(anyNA(res))
  expect_equal(res, exp(gmean + tcrossprod(alpha, W)))
})

test_that("hasBadValues detects NA/NaN/Inf (matrix and tensor)", {
  m <- matrix(as.numeric(1:6), 2)
  expect_false(hasBadValues(m))
  expect_true(hasBadValues(replace(m, 1, NA)))
  expect_true(hasBadValues(replace(m, 1, NaN)))
  expect_true(hasBadValues(replace(m, 1, Inf)))

  skip_if_no_gpu()
  expect_false(hasBadValues(toGPUMatrix(m, backend = "gpu")))
})

# ---- GPU memory budget detection and gene-block sizing -------------------

test_that(".parseNvidiaSmiFreeMem parses canned nvidia-smi output to bytes", {
  expect_equal(.parseNvidiaSmiFreeMem("1024"), 1024 * 1024^2)
  expect_equal(.parseNvidiaSmiFreeMem("0"), 0)
  expect_true(is.na(.parseNvidiaSmiFreeMem(NA_character_)))
  expect_true(is.na(.parseNvidiaSmiFreeMem("")))
  expect_true(is.na(.parseNvidiaSmiFreeMem("not a number")))
  expect_true(is.na(.parseNvidiaSmiFreeMem("-5")))
})

test_that(".parseSysctlMem parses canned sysctl output to bytes", {
  expect_equal(.parseSysctlMem("17179869184"), 17179869184)
  expect_true(is.na(.parseSysctlMem(NA_character_)))
  expect_true(is.na(.parseSysctlMem("")))
  expect_true(is.na(.parseSysctlMem("garbage")))
  expect_true(is.na(.parseSysctlMem("-1")))
})

test_that("getGPUMemoryBudget validates and honours a user override", {
  resetGPUCache()
  expect_error(getGPUMemoryBudget(gpu.mem.budget = -1), "positive")
  expect_error(getGPUMemoryBudget(gpu.mem.budget = 0), "positive")
  expect_error(getGPUMemoryBudget(gpu.mem.budget = NA_real_), "positive")
  expect_error(getGPUMemoryBudget(gpu.mem.budget = c(1, 2)), "positive")

  expect_equal(getGPUMemoryBudget(gpu.mem.budget = 12345), 12345)
  expect_equal(getGPUMemoryBudget(gpu.mem.budget = Inf), Inf) # disables blocking
  resetGPUCache()
})

test_that("getGPUMemoryBudget returns Inf and skips detection when no accelerator is in use", {
  resetGPUCache()
  testthat::local_mocked_bindings(checkGPU = function() FALSE)
  expect_equal(getGPUMemoryBudget(), Inf)
  resetGPUCache()
})

test_that("getGPUMemoryBudget uses nvidia-smi on a mocked CUDA device", {
  resetGPUCache()
  testthat::local_mocked_bindings(
    checkGPU = function() TRUE, getBackendDevice = function() "cuda",
    .nvidiaSmiFreeMemRaw = function() "8192"
  )
  expect_equal(getGPUMemoryBudget(), 8192 * 1024^2)
  resetGPUCache()
})

test_that("getGPUMemoryBudget uses sysctl (fraction of total RAM) on a mocked MPS device", {
  resetGPUCache()
  testthat::local_mocked_bindings(
    checkGPU = function() TRUE, getBackendDevice = function() "mps",
    .sysctlMemRaw = function() "17179869184"
  )
  expect_equal(getGPUMemoryBudget(), 17179869184 * 0.5)
  resetGPUCache()
})

test_that("getGPUMemoryBudget falls back to a fixed constant and messages when detection fails", {
  resetGPUCache()
  testthat::local_mocked_bindings(
    checkGPU = function() TRUE, getBackendDevice = function() "cuda",
    .nvidiaSmiFreeMemRaw = function() NA_character_
  )
  expect_message(budget <- getGPUMemoryBudget(), "Could not auto-detect")
  expect_equal(budget, GPU_MEM_BUDGET_FALLBACK_BYTES)
  resetGPUCache()
})

test_that("getGPUMemoryBudget caches the detected value for the session", {
  resetGPUCache()
  testthat::local_mocked_bindings(
    checkGPU = function() TRUE, getBackendDevice = function() "cuda",
    .nvidiaSmiFreeMemRaw = function() "1000"
  )
  first <- getGPUMemoryBudget()
  # a second call must return the cached value even if detection would now differ
  testthat::local_mocked_bindings(.nvidiaSmiFreeMemRaw = function() "999999")
  expect_equal(getGPUMemoryBudget(), first)
  resetGPUCache()
})

test_that("resetGPUCache clears the cached memory budget", {
  resetGPUCache()
  testthat::local_mocked_bindings(
    checkGPU = function() TRUE, getBackendDevice = function() "cuda",
    .nvidiaSmiFreeMemRaw = function() "1000"
  )
  getGPUMemoryBudget()
  expect_true(exists("gpu.mem.budget", envir = .spanorm_env, inherits = FALSE))
  resetGPUCache()
  expect_false(exists("gpu.mem.budget", envir = .spanorm_env, inherits = FALSE))
})

test_that("geneBlockCount never engages for backend = 'cpu', regardless of budget", {
  expect_equal(geneBlockCount(1e5, 1e5, 20, backend = "cpu", budget.bytes = 1), 1L)
})

test_that("geneBlockCount returns 1L when the problem already fits the budget", {
  expect_equal(geneBlockCount(100, 100, 10, backend = "auto", budget.bytes = Inf), 1L)
  expect_equal(geneBlockCount(100, 100, 10, backend = "auto", budget.bytes = 1e12), 1L)
})

test_that("geneBlockCount returns >1L and is clamped to n_genes under a tiny budget", {
  testthat::local_mocked_bindings(checkGPU = function() TRUE, getBackendDevice = function() "cuda")
  nb <- geneBlockCount(1000, 1000, 10, backend = "gpu", budget.bytes = 1000)
  expect_gt(nb, 1L)
  expect_lte(nb, 1000L)

  # pathologically tiny budget relative to a small gene count: still clamped
  nb.small <- geneBlockCount(5, 1e6, 10, backend = "gpu", budget.bytes = 1)
  expect_equal(nb.small, 5L)
})

test_that("geneBlockIndices splits genes into nblocks contiguous, exhaustive blocks", {
  blocks <- geneBlockIndices(10, 3)
  expect_equal(length(blocks), 3)
  expect_equal(sort(unname(unlist(blocks))), 1:10)
  expect_equal(unname(unlist(blocks)), 1:10) # contiguous and in order

  expect_equal(unname(unlist(geneBlockIndices(5, 1))), 1:5)
})

# ---- Blocked-fit helper parity against the whole-tensor computation -------
# these exercise the actual accelerator path (GPU/MPS) since the blocked
# helpers are only ever reached when checkGPU() is TRUE

test_that("blocked-fit helpers match their whole-tensor equivalents", {
  skip_if_no_gpu()

  set.seed(50)
  ng <- 17; ns <- 23; ncov <- 4
  alpha <- matrix(rnorm(ng * ncov, sd = 0.3), ng, ncov)
  gmean <- rnorm(ng, 1, 0.2)
  psi <- runif(ng, 0.1, 1)
  Y <- matrix(rpois(ng * ns, 8), ng, ns)
  W <- cbind(1, matrix(rnorm(ns * (ncov - 1), sd = 0.5), ns, ncov - 1))
  Wsub <- toGPUMatrix(W, backend = "gpu")
  step <- rep(1, ns)
  blocks <- geneBlockIndices(ng, 4) # force multiple, uneven-sized blocks

  # .blockedNBLoglik vs whole-tensor loglik
  ref.loglik <- as.numeric(toRMatrix(colSums_gpu(dnbinom_gpu(
    toGPUMatrix(Y, backend = "gpu"),
    mu = exp(add_vec_mat_gpu(gmean, tcrossprod_gpu(toGPUMatrix(alpha, backend = "gpu"), Wsub))),
    size = 1 / psi, log = TRUE
  ))))
  expect_equal(.blockedNBLoglik(alpha, gmean, Y, Wsub, psi, "gpu", blocks), ref.loglik, tolerance = gpu_tol())

  # .blockedOffsPsi vs whole-tensor offset slice
  psi.cols <- c(1, 3, 5, 7, 9)
  Wsub.psi <- Wsub[psi.cols, , drop = FALSE]
  ref.offs <- toRMatrix(tcrossprod_gpu(toGPUMatrix(alpha, backend = "gpu"), Wsub.psi))
  expect_equal(.blockedOffsPsi(alpha, Wsub.psi, blocks), ref.offs, tolerance = gpu_tol())

  # .irlsBlockedWtCell (Pass A) vs whole-tensor wt.cell + gmean-fold caches
  lmu.hat <- add_vec_mat_gpu(gmean, tcrossprod_gpu(toGPUMatrix(alpha, backend = "gpu"), Wsub))
  sig.inv <- 1 / add_vec_mat_gpu(psi, exp(-lmu.hat), backend = "gpu")
  ref.wt.cell <- as.numeric(toRMatrix(colSums_gpu(sig.inv))) / ng
  ref.wt.cell <- pmin(ref.wt.cell, quantile(ref.wt.cell, probs = 0.98))
  Zr <- (toGPUMatrix(Y, backend = "gpu") + 0.01) / (exp(lmu.hat) + 0.01) - 1
  Z <- lmu.hat + mult_vec_mat_gpu(step, Zr, backend = "gpu")
  ref.rowsum.Zsig <- as.numeric(toRMatrix(rowSums_gpu(Z * sig.inv)))
  ref.M <- toRMatrix(matmul_gpu(sig.inv, Wsub))
  ref.rowsum.sig <- as.numeric(toRMatrix(rowSums_gpu(sig.inv)))

  pass.a <- .irlsBlockedWtCell(alpha, gmean, Y, Wsub, psi, step, TRUE, "gpu", blocks)
  expect_equal(pass.a$wt.cell, ref.wt.cell, tolerance = gpu_tol())
  expect_equal(pass.a$rowsum.Zsig, ref.rowsum.Zsig, tolerance = gpu_tol())
  expect_equal(pass.a$M, ref.M, tolerance = gpu_tol())
  expect_equal(pass.a$rowsum.sig, ref.rowsum.sig, tolerance = gpu_tol())

  # is.spanorm = FALSE: only wt.cell is computed, the gmean-fold caches stay NULL
  pass.a.generic <- .irlsBlockedWtCell(alpha, gmean, Y, Wsub, psi, step, FALSE, "gpu", blocks)
  expect_equal(pass.a.generic$wt.cell, ref.wt.cell, tolerance = gpu_tol())
  expect_null(pass.a.generic$rowsum.Zsig)
  expect_null(pass.a.generic$M)
  expect_null(pass.a.generic$rowsum.sig)

  # .irlsBlockedAlphaNumerator vs whole-tensor b
  wt.cell <- runif(ns, 0.5, 2)
  Wsub.wt <- mult_vec_mat_gpu(wt.cell, Wsub, backend = "gpu")
  Z2 <- lmu.hat + mult_vec_mat_gpu(step, Zr, backend = "gpu")
  Zc <- add_vec_mat_gpu(-gmean, Z2, backend = "gpu")
  ref.b <- toRMatrix(matmul_gpu(Zc, Wsub.wt))
  expect_equal(.irlsBlockedAlphaNumerator(alpha, gmean, Y, Wsub, step, Wsub.wt, "gpu", blocks), ref.b, tolerance = gpu_tol())

  # .irlsBlockedGmeanUpdate vs whole-tensor gmean formula
  alpha.new <- alpha + 0.05 # arbitrary "updated" alpha
  Z.res <- Z2 - tcrossprod_gpu(toGPUMatrix(alpha.new, backend = "gpu"), Wsub)
  ref.gmean <- as.numeric(toRMatrix(rowSums_gpu(Z.res * sig.inv) / rowSums_gpu(sig.inv)))
  expect_equal(.irlsBlockedGmeanUpdate(alpha.new, ref.rowsum.Zsig, ref.M, ref.rowsum.sig), ref.gmean, tolerance = gpu_tol())
})

test_that("invert_mat CPU path handles SPD and non-SPD symmetric matrices", {
  set.seed(11)
  a <- matrix(rnorm(16), 4)
  spd <- crossprod(a) + diag(4) # symmetric positive definite
  expect_equal(as.matrix(invert_mat(spd)), solve(spd), tolerance = 1e-6)

  indef <- matrix(c(0, 1, 1, 0), 2) # symmetric indefinite -> Cholesky fails
  expect_equal(as.matrix(invert_mat(indef)), solve(indef), tolerance = 1e-6)
})

# a (batch, n, n) array of SPD matrices, and its per-slice reference inverse
.spd_batch <- function(batch, n, seed = 12) {
  set.seed(seed)
  arr <- array(0, c(batch, n, n))
  ref <- array(0, c(batch, n, n))
  for (i in seq_len(batch)) {
    a <- matrix(rnorm(n * n), n, n)
    m <- crossprod(a) + diag(n)
    arr[i, , ] <- m
    ref[i, , ] <- solve(m)
  }
  list(arr = arr, ref = ref)
}

test_that("invert_mat_batched CPU path matches per-slice solve()", {
  b <- .spd_batch(5, 3)
  out <- invert_mat_batched(b$arr)
  expect_equal(dim(out), c(5L, 3L, 3L))
  expect_equal(out, b$ref, tolerance = 1e-6)
})

test_that("invert_mat_batched GPU path matches per-slice invert_mat", {
  skip_if_no_gpu()

  b <- .spd_batch(5, 3)
  t <- toGPUMatrix(matrix(0, 1, 1), backend = "gpu") # force torch load
  arr.t <- torch::torch_tensor(b$arr, dtype = getBackendDtype(),
                               device = getBackendDevice())
  out <- as.array(invert_mat_batched(arr.t)$cpu())
  expect_equal(out, b$ref, tolerance = gpu_tol())
})

test_that("calculateMu GPU path matches the CPU path (incl. winsorisation)", {
  skip_if_no_gpu()

  set.seed(13)
  ng <- 8
  nc <- 20
  p <- 3
  alpha <- matrix(rnorm(ng * p), ng, p)
  W <- matrix(rnorm(nc * p), nc, p)
  gmean <- rnorm(ng)

  # default winsorisation
  cpu <- calculateMu(gmean, alpha, W, backend = "cpu")
  gpu <- toRMatrix(calculateMu(gmean, alpha, W, backend = "gpu"))
  expect_equal(gpu, cpu, tolerance = gpu_tol())

  # winsor = Inf (no clamp)
  cpu.inf <- calculateMu(gmean, alpha, W, winsor = Inf, backend = "cpu")
  gpu.inf <- toRMatrix(calculateMu(gmean, alpha, W, winsor = Inf, backend = "gpu"))
  expect_equal(gpu.inf, cpu.inf, tolerance = gpu_tol())

  # mad == 0 rows (identical rows) must not produce NaN on either path
  alpha0 <- matrix(c(1, 1, 2, 2), 2, 2)
  W0 <- diag(2)
  cpu0 <- calculateMu(0, alpha0, W0, backend = "cpu")
  gpu0 <- toRMatrix(calculateMu(0, alpha0, W0, backend = "gpu"))
  expect_false(hasBadValues(gpu0))
  expect_equal(gpu0, cpu0, tolerance = gpu_tol())
})
