skip_if_no_gpu <- function() {
  testthat::skip_if_not(checkGPU(), "torch GPU/MPS not available")
}

# accelerator tolerance: near machine precision on float64 devices (cuda/cpu),
# looser on MPS which can only run float32
gpu_tol <- function() {
  if (getBackendDevice() == "mps") 1e-6 else 1e-10
}

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

test_that("invert_mat CPU path handles SPD and non-SPD symmetric matrices", {
  set.seed(11)
  a <- matrix(rnorm(16), 4)
  spd <- crossprod(a) + diag(4) # symmetric positive definite
  expect_equal(as.matrix(invert_mat(spd)), solve(spd), tolerance = 1e-6)

  indef <- matrix(c(0, 1, 1, 0), 2) # symmetric indefinite -> Cholesky fails
  expect_equal(as.matrix(invert_mat(indef)), solve(indef), tolerance = 1e-6)
})
