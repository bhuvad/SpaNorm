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
  if (length(unique(spe$sample_id)) > 1) {
    stop("multiple samples/images detected in SpatialExperiment object; only single-sample objects are supported")
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

  # Confirm only a single sample is present
  # 1) Multiple images typically indicate multiple samples in ST contexts
  if (!is.null(spe@images) && length(spe@images) > 1) {
    stop("multiple samples/images detected in Seurat object; only single-sample objects are supported")
  }
  # 2) Use common sample-identifying columns from meta.data
  sample_cols <- c("orig.ident", "sample", "Sample", "sample_id", "dataset")
  md <- spe@meta.data
  sample_col <- intersect(sample_cols, colnames(md))
  if (length(sample_col) > 0) {
    nsamp <- length(unique(md[[sample_col[1]]]))
    if (nsamp != 1) {
      stop(sprintf("multiple samples detected by '%s' (n = %d); only single-sample objects are supported", sample_col[1], nsamp))
    }
  }

  # validate spatial coordinates exist and align to cells/spots
  coords <- extractSeuratCoords(spe)
  if (is.null(coords)) {
    stop("spatial coordinates not found in Seurat object (images or meta.data)")
  }
}

# Extract spatial coordinates from a Seurat object when available
# Tries images slot (preferred), then meta.data common column names.
extractSeuratCoords <- function(spe) {
  coords <- NULL

  # try images slot
  if (!is.null(spe@images) && length(spe@images) > 0) {
    coord_df <- spe@images[[1]]@coordinates
    if (is.data.frame(coord_df) || is.matrix(coord_df)) {
      cn <- colnames(coord_df)
      if (all(c("x", "y") %in% cn)) {
        coords <- coord_df[, c("x", "y"), drop = FALSE]
      } else if (all(c("imagecol", "imagerow") %in% cn)) {
        coords <- coord_df[, c("imagecol", "imagerow"), drop = FALSE]
      } else if (ncol(coord_df) >= 2) {
        coords <- coord_df[, 1:2, drop = FALSE]
      }
    }
  }

  # try meta.data if not found
  if (is.null(coords)) {
    md <- spe@meta.data
    if (!is.null(md)) {
      cn <- colnames(md)
      if (all(c("x", "y") %in% cn)) {
        coords <- md[, c("x", "y"), drop = FALSE]
      } else if (all(c("imagecol", "imagerow") %in% cn)) {
        coords <- md[, c("imagecol", "imagerow"), drop = FALSE]
      }
    }
  }

  # if still not found, return NULL
  if (is.null(coords)) return(NULL)

  # ensure matrix
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  # ensure numeric and finite
  if (!is.numeric(coords)) coords <- apply(coords, 2, as.numeric)
  if (any(!is.finite(coords))) {
    stop("spatial coordinates contain non-finite values")
  }
  # dimensions should match number of cells/spots in counts layer
  ncells <- ncol(SeuratObject::LayerData(spe, "counts"))
  if (nrow(coords) != ncells) {
    stop("number of coordinate rows does not match number of cells/spots")
  }

  return(coords)
}

checkBatch <- function(batch, nobs) {
  # if null, do nothing
  if (is.null(batch)) {
    batch <- c()
  }
  
  # check missing
  if (any(is.na(batch))) {
    stop("'batch' cannot have missing values")
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
