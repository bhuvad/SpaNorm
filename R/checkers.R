checkSPE <- function(spe) {
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(nrow(spe) > 0 & ncol(spe) > 0)

  if (!"counts" %in% SummarizedExperiment::assayNames(spe)) {
    stop("assay 'counts' not found in 'spe'")
  } else {
    expr <- SummarizedExperiment::assay(spe, "counts")
    expr <- ifelse(is(expr, "sparseMatrix"), expr@x, as.numeric(expr))
    if (any(expr %% 1 != 0)) {
      stop("'counts' has non-integer values")
    }
    if (any(expr < 0)) {
      stop("'counts' has negative values")
    }
  }
}

checkSeurat <- function(spe) {
  stopifnot(is(spe, "Seurat"))
  stopifnot(nrow(spe) > 0)

  if (!"counts" %in% SeuratObject::Layers(spe)) {
    stop("layer 'counts' not found in 'spe'")
  } else {
    expr <- SeuratObject::LayerData(spe, "counts")
    expr <- ifelse(is(expr, "sparseMatrix"), expr@x, as.numeric(expr))
    if (any(expr %% 1 != 0)) {
      stop("'counts' has non-integer values")
    }
    if (any(expr < 0)) {
      stop("'counts' has negative values")
    }
  }
}

checkBatch <- function(batch, nobs) {
  # if null, do nothing
  if (is.null(batch)) {
    batch <- c()
  } else {
    if (any(is.na(batch))) {
      stop("'batch' cannot have missing values")
    }
  }

  # if matrix, check and return
  if (is.matrix(batch)) {
    # check dimensions
    if (nrow(batch) != nobs) {
      stop("number of rows in the 'batch' matrix do not match number of cells/spots")
    }

    # check type
    if (!is.numeric(batch)) {
      stop("'batch' should be a numeric matrix (consider using 'model.matrix' to define the design)")
    }

    # check for intercept
    isintercept <- grepl("intercept", colnames(batch), ignore.case = TRUE)
    isintercept <- isintercept | matrixStats::colAlls(batch == 1)
    if (any(isintercept)) {
      warning("'intercept' term detected and will be removed")
      batch <- batch[, !isintercept, drop = FALSE]
    }
  } else if (is(batch, "vector_OR_factor")) {
    # check dimensions
    if (length(batch) != nobs) {
      stop("length of 'batch' vector does not match number of cells/spots")
    }
    batch <- model.matrix(~batch)[, -1, drop = FALSE]
  }

  return(batch)
}
