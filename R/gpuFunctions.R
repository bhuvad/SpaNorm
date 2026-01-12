#' @importFrom stats dnbinom qnbinom pnbinom mad median quantile model.matrix
NULL

checkGPU <- function() {
  # check if gpuR is installed and GPUs are available
  if (requireNamespace("tensorflow", quietly = TRUE)) {
    res = tryCatch({
      if (length(tensorflow::tf$config$experimental$list_physical_devices("GPU")) > 0) {
        return(TRUE)
      } else {
        FALSE
      }
    }, error = \(e) {
      return(FALSE)
    })
    return(res)
  }
  return(FALSE)
}

# helper to detect TensorFlow tensors
is_tf_tensor <- function(x) {
  if (!requireNamespace("tensorflow", quietly = TRUE)) return(FALSE)
  out <- tryCatch(tensorflow::tf$is_tensor(x), error = function(e) FALSE)
  isTRUE(out)
}

toGPUMatrix <- function(mat, ..., backend = c("auto", "cpu", "gpu")) {
  backend = match.arg(backend)
  # convert matrix to GPU matrix if GPUs are available
  if (checkGPU() && backend %in% c("gpu", "auto")) {
    if (is_tf_tensor(mat)) {
      return(mat)
    }
    if (is.matrix(mat)) {
      return(tensorflow::tf$constant(mat, dtype = tensorflow::tf$float32))
    }
  }

  mat
}

toGPUVector <- function(vec, n = NULL, backend = c("auto", "cpu", "gpu")) {
  backend = match.arg(backend)
  # convert vector to GPU (TensorFlow) vector if GPUs are available
  if (checkGPU() && backend %in% c("gpu", "auto")) {
    if (is_tf_tensor(vec)) {
      return(vec)
    }
    if (is.null(n) || length(vec) == n) {
      return(tensorflow::tf$constant(as.numeric(vec), dtype = tensorflow::tf$float32))
    } else if (length(vec) == 1) {
      val <- tensorflow::tf$constant(as.numeric(vec), dtype = tensorflow::tf$float32)
      ones <- tensorflow::tf$ones(list(n), dtype = tensorflow::tf$float32)
      return(tensorflow::tf$math$multiply(ones, val))
    } else {
      stop("Length of vec does not match n.")
    }
  }

  vec
}

diag_mat <- function(vec, backend = c("auto", "cpu", "gpu")) {
  backend = match.arg(backend)
  if (checkGPU() && backend %in% c("gpu", "auto")) {
    v <- if (is_tf_tensor(vec)) vec else tensorflow::tf$constant(as.numeric(vec), dtype = tensorflow::tf$float32)
    mat <- tensorflow::tf$linalg$diag(v)
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
  if (checkGPU() && is_tf_tensor(mat)) {
    res <- tryCatch({
      chol <- tensorflow::tf$linalg$cholesky(mat)
      n <- mat$shape$as_list()[[1]]
      dtype <- tryCatch(mat$dtype, error = function(e) tensorflow::tf$float32)
      I <- tensorflow::tf$eye(as.integer(n), dtype = dtype)
      tensorflow::tf$linalg$cholesky_solve(chol, I)
    }, error = function(e) {
      tensorflow::tf$linalg$inv(mat)
    })
    return(res)
  }

  # fallback to base R (uses Cholesky)
  return(Matrix::chol2inv(Matrix::chol(mat)))
}

rowMeans_gpu <- function(mat) {
  # calculate row means using GPU if available
  if (checkGPU() && is_tf_tensor(mat)) {
    res <- tensorflow::tf$reduce_mean(mat, axis = 1L)
    return(res)
  }

  # fallback to base R
  return(matrixStats::rowMeans2(mat))
}

colMeans_gpu <- function(mat) {
  # calculate row means using GPU if available
  if (checkGPU() && is_tf_tensor(mat)) {
    res <- tensorflow::tf$reduce_mean(mat, axis = 0L)
    return(res)
  }

  # fallback to base R
  return(matrixStats::colMeans2(mat))
}

rowSums_gpu <- function(mat) {
  # calculate row sums using GPU if available
  if (checkGPU() && is_tf_tensor(mat)) {
    res <- tensorflow::tf$reduce_sum(mat, axis = 1L)
    return(res)
  }

  # fallback to base R
  return(matrixStats::rowSums2(mat))
}

colSums_gpu <- function(mat) {
  # calculate row sums using GPU if available
  if (checkGPU() && is_tf_tensor(mat)) {
    res <- tensorflow::tf$reduce_sum(mat, axis = 0L)
    return(res)
  }

  # fallback to base R
  return(matrixStats::colSums2(mat))
}

dnbinom_gpu <- function(x, mu, size, log = FALSE) {
  dnbinom(as.matrix(x), mu = as.matrix(mu), size = as.vector(size), log = TRUE)
}

qnbinom_gpu <- function(p, mu, size, log = FALSE) {
  qnbinom(as.matrix(p), mu = as.matrix(mu), size = as.vector(size))
}

pnbinom_gpu <- function(q, mu, size, log = FALSE) {
  pnbinom(as.matrix(q), mu = as.matrix(mu), size = as.vector(size))
}

copy <- function(x) {
  # clone object
  if (checkGPU() && is_tf_tensor(x)) {
    return(tensorflow::tf$identity(x))
  }

  # fallback to base R
  return(x)
}
