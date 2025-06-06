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

toGPUMatrix <- function(mat, ...) {
  # convert matrix to GPU matrix if GPUs are available
  if (checkGPU() && !inherits(mat, "gpuRmatrix")) {
    if (is(mat, "vclVector") || is(mat, "gpuVector") || is(mat, "vector")) {
      mat <- matrix(mat, ...)
    }
    mat <- gpuR::vclMatrix(mat, type = "float")
  }

  mat
}

toGPUVector <- function(vec, n = NULL) {
  # convert matrix to GPU matrix if GPUs are available
  if (checkGPU() && !inherits(vec, "vclVector")) {
    if (is.null(n) || length(vec) == n) {
      vec <- gpuR::vclVector(as.vector(vec), type = "float")
    } else if (length(vec) == 1) {
      vec <- gpuR::vclVector(as.vector(vec), length = n, type = "float")
    } else {
      stop("Length of vec does not match n.")
    }
  }

  vec
}

diag_mat <- function(vec) {
  if (checkGPU()) {
    mat = gpuR::vclMatrix(nrow = length(vec), ncol = length(vec), type = "float")
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
    return(gpuR::rowMeans(mat))
  }

  # fallback to base R
  return(matrixStats::rowMeans2(mat))
}

colMeans_gpu <- function(mat) {
  # calculate row means using GPU if available
  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    return(gpuR::colMeans(mat))
  }

  # fallback to base R
  return(matrixStats::colMeans2(mat))
}

rowSums_gpu <- function(mat) {
  # calculate row sums using GPU if available
  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    return(gpuR::rowSums(mat))
  }

  # fallback to base R
  return(matrixStats::rowSums2(mat))
}

colSums_gpu <- function(mat) {
  # calculate row sums using GPU if available
  if (checkGPU() && inherits(mat, "gpuRmatrix")) {
    return(gpuR::colSums(mat))
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
  if (checkGPU() && (inherits(x, "gpuRmatrix") || inherits(x, "gpuVector") || inherits(x, "vclMatrix"))) {
    return(gpuR::deepcopy(x))
  }

  # fallback to base R
  return(x)
}

setMethod("+", signature(e1 = "gpuMatrix", e2 = "gpuVector"), \(e1, e2) {
  e1 + tcrossprod(e2, gpuR::gpuVector(1, ncol(e1)))
})
setMethod("+", signature(e1 = "gpuVector", e2 = "gpuMatrix"), \(e1, e2) e2 + e1)
setMethod("+", signature(e1 = "vclMatrix", e2 = "vclVector"), \(e1, e2) {
  e1 + tcrossprod(e2, gpuR::vclVector(1, ncol(e1)))
})
setMethod("+", signature(e1 = "vclVector", e2 = "vclMatrix"), \(e1, e2) e2 + e1)

setMethod("-", signature(e1 = "gpuMatrix", e2 = "gpuVector"), \(e1, e2) {
  e1 - tcrossprod(e2, gpuR::gpuVector(1, ncol(e1)))
})
setMethod("-", signature(e1 = "gpuVector", e2 = "gpuMatrix"), \(e1, e2) e2 - e1)
setMethod("-", signature(e1 = "vclMatrix", e2 = "vclVector"), \(e1, e2) {
  e1 - tcrossprod(e2, gpuR::vclVector(1, ncol(e1)))
})
setMethod("-", signature(e1 = "vclVector", e2 = "vclMatrix"), \(e1, e2) e2 - e1)

setMethod("*", signature(e1 = "gpuMatrix", e2 = "gpuVector"), \(e1, e2) {
  e1 * tcrossprod(e2, gpuR::gpuVector(1, ncol(e1)))
})
setMethod("*", signature(e1 = "gpuVector", e2 = "gpuMatrix"), \(e1, e2) e2 * e1)
setMethod("*", signature(e1 = "vclMatrix", e2 = "vclVector"), \(e1, e2) {
  e1 * tcrossprod(e2, gpuR::vclVector(1, ncol(e1)))
})
setMethod("*", signature(e1 = "vclVector", e2 = "vclMatrix"), \(e1, e2) e2 * e1)

setMethod("/", signature(e1 = "gpuMatrix", e2 = "gpuVector"), \(e1, e2) {
  e1 / tcrossprod(e2, gpuR::gpuVector(1, ncol(e1)))
})
setMethod("/", signature(e1 = "gpuVector", e2 = "gpuMatrix"), \(e1, e2) e2 / e1)
setMethod("/", signature(e1 = "vclMatrix", e2 = "vclVector"), \(e1, e2) {
  e1 / tcrossprod(e2, gpuR::vclVector(1, ncol(e1)))
})
setMethod("/", signature(e1 = "vclVector", e2 = "vclMatrix"), \(e1, e2) e2 / e1)
