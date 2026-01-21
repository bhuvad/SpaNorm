skip_if_no_gpu <- function() {
  testthat::skip_if_not(checkGPU(), "TensorFlow GPU not available")
}

as_numeric_matrix <- function(x) {
  if (is_tf_tensor(x)) {
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

  expect_true(is_tf_tensor(mgpu))
  expect_true(is_tf_tensor(vgpu))

  expect_equal(toRMatrix(mgpu), m, tolerance = 1e-6)
  expect_equal(as_numeric_vector(vgpu), as.numeric(v), tolerance = 1e-6)
})

test_that("GPU tcrossprod and crossprod match CPU", {
  skip_if_no_gpu()

  set.seed(2)
  m <- matrix(rnorm(12), nrow = 3)

  tc_gpu <- toRMatrix(tcrossprod_gpu(toGPUMatrix(m, backend = "gpu")))
  cp_gpu <- toRMatrix(crossprod_gpu(toGPUMatrix(m, backend = "gpu")))

  expect_equal(tc_gpu, tcrossprod(m), tolerance = 1e-6)
  expect_equal(cp_gpu, crossprod(m), tolerance = 1e-6)
})

test_that("Vector-matrix add and multiply broadcast match CPU", {
  skip_if_no_gpu()

  set.seed(3)
  m <- matrix(rnorm(15), nrow = 5)
  v_row <- rnorm(5) # length matches nrow
  v_col <- rnorm(3) # length matches ncol

  add_gpu <- toRMatrix(add_vec_mat_gpu(v_row, toGPUMatrix(m, backend = "gpu")))
  mult_gpu <- toRMatrix(mult_vec_mat_gpu(v_col, toGPUMatrix(m, backend = "gpu")))

  expect_equal(add_gpu, m + v_row, tolerance = 1e-6)
  expect_equal(mult_gpu, sweep(m, 2, v_col, `*`), tolerance = 1e-6)
})

test_that("diag_mat and invert_mat match CPU", {
  skip_if_no_gpu()

  set.seed(4)
  v <- rnorm(4)
  base_diag <- diag(v)
  gpu_diag <- toRMatrix(diag_mat(v, backend = "gpu"))
  expect_equal(gpu_diag, base_diag, tolerance = 1e-6)

  # Test length-1 vector (edge case)
  v1 <- 2.5
  base_diag1 <- diag(v1, nrow = 1, ncol = 1)
  gpu_diag1 <- toRMatrix(diag_mat(v1, backend = "gpu"))
  expect_equal(gpu_diag1, base_diag1, tolerance = 1e-6)

  a <- matrix(rnorm(16), nrow = 4)
  pd <- crossprod(a) + diag(4) * 0.1
  inv_gpu <- toRMatrix(invert_mat(toGPUMatrix(pd, backend = "gpu")))
  expect_equal(inv_gpu, solve(pd), tolerance = 1e-6)
})

test_that("Row/col summaries match CPU", {
  skip_if_no_gpu()

  set.seed(5)
  m <- matrix(rnorm(18), nrow = 6)
  mgpu <- toGPUMatrix(m, backend = "gpu")

  expect_equal(as_numeric_vector(rowMeans_gpu(mgpu)), rowMeans(m), tolerance = 1e-6)
  expect_equal(as_numeric_vector(colMeans_gpu(mgpu)), colMeans(m), tolerance = 1e-6)
  expect_equal(as_numeric_vector(rowSums_gpu(mgpu)), rowSums(m), tolerance = 1e-6)
  expect_equal(as_numeric_vector(colSums_gpu(mgpu)), colSums(m), tolerance = 1e-6)
})

test_that("Tensor-aware %*% matches CPU", {
  skip_if_no_gpu()

  x <- 1:3
  y <- 4:6
  xy_gpu <- as_numeric_vector(toGPUVector(x, backend = "gpu") %*% toGPUVector(y, backend = "gpu"))
  expect_equal(xy_gpu, sum(x * y), tolerance = 1e-6)

  mat <- matrix(1:6, nrow = 2)
  vec <- c(1, 2, 3)
  mat_gpu <- toGPUMatrix(mat, backend = "gpu")
  vec_gpu <- toGPUVector(vec, backend = "gpu")
  mv_gpu <- as_numeric_vector(mat_gpu %*% vec_gpu)
  expect_equal(mv_gpu, as.numeric(mat %*% vec), tolerance = 1e-6)
})
