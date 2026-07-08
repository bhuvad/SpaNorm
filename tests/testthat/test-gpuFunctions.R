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

test_that("invert_mat CPU path handles SPD and non-SPD symmetric matrices", {
  set.seed(11)
  a <- matrix(rnorm(16), 4)
  spd <- crossprod(a) + diag(4) # symmetric positive definite
  expect_equal(as.matrix(invert_mat(spd)), solve(spd), tolerance = 1e-6)

  indef <- matrix(c(0, 1, 1, 0), 2) # symmetric indefinite -> Cholesky fails
  expect_equal(as.matrix(invert_mat(indef)), solve(indef), tolerance = 1e-6)
})
