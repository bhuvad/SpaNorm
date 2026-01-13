#' @importFrom stats dnbinom qnbinom pnbinom mad median quantile model.matrix
NULL

checkGPU <- function() {
  # check if tensorflow is installed and GPUs are available
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

# Convert a TensorFlow tensor or R object to a base R matrix
toRMatrix <- function(x) {
  if (is_tf_tensor(x)) {
    arr <- NULL
    if (is.null(arr)) {
      # Fallback: try eager tensor numpy() conversion
      arr <- tryCatch(x$numpy(), error = function(e) NULL)
    }
    if (is.null(arr)) {
      stop("Failed to convert TensorFlow tensor to R matrix")
    }
    dims <- dim(arr)
    if (is.null(dims)) {
      # scalar
      return(matrix(as.numeric(arr), nrow = 1L, ncol = 1L))
    } else if (length(dims) == 1) {
      # vector -> column matrix
      return(matrix(arr, nrow = dims[1], ncol = 1L))
    } else {
      return(as.matrix(arr))
    }
  }

  # Non-tensor input: standard coercion
  return(as.matrix(x))
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

tcrossprod_gpu <- function(x, y = NULL) {
  # compute tcrossprod (x %*% t(y)) using GPU if available
  if (checkGPU() && (is_tf_tensor(x) || (!is.null(y) && is_tf_tensor(y)))) {
    xt <- if (is_tf_tensor(x)) x else tensorflow::tf$constant(as.matrix(x), dtype = tensorflow::tf$float32)
    if (is.null(y)) {
      yt <- xt
    } else {
      yt <- if (is_tf_tensor(y)) y else tensorflow::tf$constant(as.matrix(y), dtype = tensorflow::tf$float32)
    }
    res <- tensorflow::tf$matmul(xt, tensorflow::tf$transpose(yt))
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

  if (checkGPU() && backend %in% c("gpu", "auto") && (is_tf_tensor(vec) || is_tf_tensor(mat))) {
    v <- if (is_tf_tensor(vec)) vec else tensorflow::tf$constant(as.numeric(vec), dtype = tensorflow::tf$float32)
    m <- if (is_tf_tensor(mat)) mat else tensorflow::tf$constant(as.matrix(mat), dtype = tensorflow::tf$float32)

    # Determine broadcasting orientation
    mshape <- m$shape$as_list()
    m_nr <- tryCatch(as.integer(mshape[[1]]), error = function(e) NULL)
    m_nc <- tryCatch(as.integer(mshape[[2]]), error = function(e) NULL)
    vshape <- v$shape$as_list()

    if (length(vshape) == 1) {
      vlen <- tryCatch(as.integer(vshape[[1]]), error = function(e) NULL)
      if (!is.null(m_nr) && !is.null(vlen) && vlen == m_nr) {
        v <- tensorflow::tf$expand_dims(v, axis = 1L) # (nr, 1)
      } else if (!is.null(m_nc) && !is.null(vlen) && vlen == m_nc) {
        v <- tensorflow::tf$expand_dims(v, axis = 0L) # (1, nc)
      } else if (!is.null(vlen) && vlen == 1L) {
        # scalar-like
        v <- tensorflow::tf$reshape(v, shape = list(1L))
      } else {
        # default to row-wise orientation
        v <- tensorflow::tf$expand_dims(v, axis = 1L)
      }
    }

    if (op == "add") {
      return(tensorflow::tf$math$add(m, v))
    } else {
      return(tensorflow::tf$math$multiply(m, v))
    }
  }

  # CPU fallback uses R's standard recycling
  if (op == "add") {
    return(as.matrix(mat) + as.vector(vec))
  } else {
    return(as.matrix(mat) * as.vector(vec))
  }
}

# Cross-product using tcrossprod_gpu: computes t(x) %*% y
crossprod_gpu <- function(x, y = NULL) {
  # compute crossprod using GPU if available by reusing tcrossprod_gpu
  if (checkGPU() && (is_tf_tensor(x) || (!is.null(y) && is_tf_tensor(y)))) {
    x_t <- if (is_tf_tensor(x)) tensorflow::tf$transpose(x) else base::t(as.matrix(x))
    if (is.null(y)) {
      y_t <- x_t
    } else {
      y_t <- if (is_tf_tensor(y)) tensorflow::tf$transpose(y) else base::t(as.matrix(y))
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

# Helper: convert object to 2D TensorFlow tensor with orientation
.to_tf_2d <- function(obj, is_left = TRUE) {
  if (is_tf_tensor(obj)) {
    shape_list <- tryCatch(obj$shape$as_list(), error = function(e) list())
    rank <- length(shape_list)
    if (rank == 1) {
      len <- as.integer(shape_list[[1]])
      if (is_left) {
        return(tensorflow::tf$reshape(obj, shape = list(1L, len)))
      } else {
        return(tensorflow::tf$reshape(obj, shape = list(len, 1L)))
      }
    }
    return(obj)
  }

  if (is.matrix(obj)) {
    return(tensorflow::tf$constant(as.matrix(obj), dtype = tensorflow::tf$float32))
  }

  v <- as.numeric(obj)
  if (is_left) {
    return(tensorflow::tf$reshape(tensorflow::tf$constant(v, dtype = tensorflow::tf$float32), shape = list(1L, as.integer(length(v)))) )
  } else {
    return(tensorflow::tf$reshape(tensorflow::tf$constant(v, dtype = tensorflow::tf$float32), shape = list(as.integer(length(v)), 1L)) )
  }
}

"%*%" <- function(x, y) {
  # Tensor-aware matrix multiplication. Falls back to base %*%.
  if (checkGPU() && (is_tf_tensor(x) || is_tf_tensor(y))) {
    ax <- .to_tf_2d(x, is_left = TRUE)
    by <- .to_tf_2d(y, is_left = FALSE)
    res <- tensorflow::tf$matmul(ax, by)

    # Mimic R shape behavior when vectors involved
    x_is_vec <- (!is.matrix(x)) && !is_tf_tensor(x) && is.atomic(x) && is.null(dim(x))
    y_is_vec <- (!is.matrix(y)) && !is_tf_tensor(y) && is.atomic(y) && is.null(dim(y))
    if (is_tf_tensor(x)) {
      x_is_vec <- tryCatch(length(x$shape$as_list()) == 1, error = function(e) FALSE)
    }
    if (is_tf_tensor(y)) {
      y_is_vec <- tryCatch(length(y$shape$as_list()) == 1, error = function(e) FALSE)
    }

    if (x_is_vec && y_is_vec) {
      # return scalar (inner product)
      res <- tensorflow::tf$squeeze(res)
    } else if (x_is_vec && !y_is_vec) {
      # row vector result (1, ncol) -> squeeze axis 0
      res <- tensorflow::tf$squeeze(res, axis = list(0L))
    } else if (!x_is_vec && y_is_vec) {
      # column vector result (nrow, 1) -> squeeze axis 1
      res <- tensorflow::tf$squeeze(res, axis = list(1L))
    }
    return(res)
  }

  # CPU/base fallback
  base::`%*%`(x, y)
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
