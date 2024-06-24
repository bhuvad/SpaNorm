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

checkNBParams <- function(ngenes, ncells, W, gmean, alpha, psi) {
  if (ncells != nrow(W)) {
    stop("nrow of 'W' does not match the number of cells/spots")
  }
  if (ngenes != length(gmean)) {
    stop("length of 'gmean' does not match the number of genes/features")
  }
  if (ngenes != length(psi)) {
    stop("length of 'psi' does not match the number of genes/features")
  }
  if (ngenes != nrow(alpha)) {
    stop("nrow of 'alpha' does not match the number of genes/features")
  }
  if (ncol(alpha) != ncol(W)) {
    stop("ncol of 'alpha' and 'W' do not match")
  }
  TRUE
}

validSpaNormSPE <- function(spe, gene.model) {
  valid <- FALSE

  # check validity based on model
  if (gene.model == "nb") {
    valid <- validSpaNormNBSPE(spe)
  }

  return(valid)
}

validSpaNormNBSPE <- function(spe) {
  valid <- FALSE

  # is object present in metadata
  if ("SpaNorm" %in% names(S4Vectors::metadata(spe))) {
    # is object valid
    fit.spanorm <- S4Vectors::metadata(spe)$SpaNorm
    valid <- checkNBParams(
      nrow(spe),
      ncol(spe),
      fit.spanorm$W,
      fit.spanorm$gmean,
      fit.spanorm$alpha,
      fit.spanorm$psi
    )
    # check isbio vector
    valid <- valid & "isbio" %in% names(fit.spanorm) & length(fit.spanorm$isbio) == ncol(fit.spanorm$W)
    # check df.tps is valis
    valid <- valid & "df.tps" %in% names(fit.spanorm) & ncol(fit.spanorm$W) >= fit.spanorm$df.tps^2 * 2 + 1 & sum(fit.spanorm$isbio) >= fit.spanorm$df.tps

    if (!valid) {
      warning("an invalid SpaNorm fit exists and will be replaced")
    }
  }

  return(valid)
}