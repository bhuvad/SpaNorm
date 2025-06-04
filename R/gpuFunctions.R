#' @importFrom stats dnbinom qnbinom pnbinom mad median quantile model.matrix
NULL

checkGPU <- function() {
  # check if gpuR is installed and GPUs are available
  if (requireNamespace("gpuR", quietly = TRUE)) {
    if (gpuR::detectGPUs() > 0) {
      return(TRUE)
    }
  }
  return(FALSE)
}

toGPUMatrix <- function(mat) {
  # convert matrix to GPU matrix if GPUs are available
  if (checkGPU() && !inherits(mat, "gpuRmatrix")) {
    mat <- gpuR::gpuMatrix(mat, type = "float")
  }

  mat
}

diag_mat <- function(vec) {
  if (checkGPU()) {
    mat = gpuR::gpuMatrix(nrow = length(vec), ncol = length(vec), type = "float")
    for (i in seq_along(vec)) {
      mat[i, i] = vec[i]
    }
  } else {
    # fallback to base R
    mat = Matrix::sparseMatrix(
      i = seq_along(vec),
      j = seq_along(vec),
      x = vec
    )
  }

  return(mat)
}

invert_mat <- function(mat) {
  # invert matrix using GPU if available
  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    mat <- gpuR::solve(mat)
    return(mat)
  }

  # fallback to base R
  return(Matrix::chol2inv(Matrix::chol(mat)))
}

rowMeans_gpu <- function(mat) {
  # calculate row means using GPU if available
  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    return(as.vector(gpuR::rowMeans(mat)))
  }

  # fallback to base R
  return(matrixStats::rowMeans2(mat))
}

colMeans_gpu <- function(mat) {
  # calculate row means using GPU if available
  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    return(as.vector(gpuR::colMeans(mat)))
  }

  # fallback to base R
  return(matrixStats::colMeans2(mat))
}

rowSums_gpu <- function(mat) {
  # calculate row sums using GPU if available
  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    return(as.vector(gpuR::rowSums(mat)))
  }

  # fallback to base R
  return(matrixStats::rowSums2(mat))
}

colSums_gpu <- function(mat) {
  # calculate row sums using GPU if available
  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    return(as.vector(gpuR::colSums(mat)))
  }

  # fallback to base R
  return(matrixStats::colSums2(mat))
}

dnbinom_gpu <- function(x, mu, size, log = FALSE) {
  dnbinom(as.matrix(x), mu = as.matrix(mu), size = size, log = TRUE)
}

qnbinom_gpu <- function(p, mu, size, log = FALSE) {
  qnbinom(as.matrix(p), mu = as.matrix(mu), size = size)
}

pnbinom_gpu <- function(q, mu, size, log = FALSE) {
  pnbinom(as.matrix(q), mu = as.matrix(mu), size = size)
}

copy <- function(x) {
  # clone object
  if (checkGPU() && (inherits(x, "gpuRmatrix") || inherits(x, "gpuVector") || inherits(x, "vclMatrix"))) {
    return(gpuR::deepcopy(x))
  }

  # fallback to base R
  return(x)
}

mat_vec_prod <- function(mat, vec) {
  vec = as.vector(vec)

  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    mat = mat * toGPUMatrix(matrix(vec, nrow(mat), ncol(mat)))
    return(mat)
  } else {
    # fallback to base R
    mat = mat * vec
  }

  return(mat)
}
