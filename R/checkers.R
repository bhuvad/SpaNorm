checkSPE <- function(spe) {
  stopifnot(is(spe, "SpatialExperiment"))
  stopifnot(nrow(spe) > 0 & ncol(spe) > 0)

  if (!"counts" %in% SummarizedExperiment::assayNames(spe)) {
    stop("assay 'counts' not found in 'spe'")
  } else {
    expr <- SummarizedExperiment::assay(spe, "counts")
    expr <- if (is(expr, "sparseMatrix")) expr@x else as.numeric(expr)
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

checkSeurat <- function(spe, assay = NULL) {
  stopifnot(is(spe, "Seurat"))
  stopifnot(nrow(spe) > 0)

  # default to the object's active assay; SpaNorm models raw counts, so the
  # caller can point this at the raw-count assay explicitly (GitHub #19)
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(spe)
  if (!assay %in% SeuratObject::Assays(spe)) {
    stop(sprintf("assay '%s' not found in the Seurat object", assay))
  }

  if (!"counts" %in% SeuratObject::Layers(spe, assay = assay)) {
    stop(sprintf("layer 'counts' not found in assay '%s'", assay))
  } else {
    expr <- SeuratObject::LayerData(spe, layer = "counts", assay = assay)
    expr <- if (is(expr, "sparseMatrix")) expr@x else as.numeric(expr)
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
  coords <- extractSeuratCoords(spe, assay = assay)
  if (is.null(coords)) {
    stop("spatial coordinates not found in Seurat object (images or meta.data)")
  }
}

# Extract spatial coordinates from a Seurat object when available.
# Uses SeuratObject::GetTissueCoordinates(), which works across image classes
# (VisiumV1, VisiumV2/FOV, SlideSeq, ...). The previous version read the image
# @coordinates slot directly, which does not exist on VisiumV2 objects and threw
# `no slot of name "coordinates"` (GitHub #18). Falls back to coordinates stored
# as meta.data columns.
extractSeuratCoords <- function(spe, assay = NULL) {
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(spe)
  cell_ids <- colnames(SeuratObject::LayerData(spe, layer = "counts", assay = assay))
  coords <- NULL

  # preferred: pull coordinates from the image(s) via the accessor
  if (!is.null(spe@images) && length(spe@images) > 0) {
    coord_df <- tryCatch(SeuratObject::GetTissueCoordinates(spe),
                         error = function(e) NULL)
    coords <- .coordsFromDf(coord_df, cell_ids)
  }

  # fallback: coordinates stored as meta.data columns
  if (is.null(coords)) {
    coords <- .coordsFromDf(spe@meta.data, cell_ids)
  }

  # if still not found, return NULL
  if (is.null(coords)) return(NULL)

  # ensure finite and aligned to the cells/spots in the counts layer
  if (any(!is.finite(coords))) {
    stop("spatial coordinates contain non-finite values")
  }
  if (nrow(coords) != length(cell_ids)) {
    stop("number of coordinate rows does not match number of cells/spots")
  }

  return(coords)
}

# Pull an (x, y) coordinate matrix out of a data.frame/matrix of coordinates,
# aligning rows to the count matrix's cell order (via a 'cell' column, else
# rownames) when the identifiers are available. Returns NULL if no usable pair
# of coordinate columns is found.
.coordsFromDf <- function(df, cell_ids) {
  if (is.null(df) || !(is.data.frame(df) || is.matrix(df))) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # align rows to the counts cell order when identifiers are present.
  # GetTissueCoordinates() returns a 'cell' column (rownames may be integers).
  if (!is.null(cell_ids)) {
    if ("cell" %in% colnames(df) && all(cell_ids %in% as.character(df$cell))) {
      rownames(df) <- as.character(df$cell)
    }
    if (!is.null(rownames(df)) && all(cell_ids %in% rownames(df))) {
      df <- df[cell_ids, , drop = FALSE]
    }
  }

  cn <- colnames(df)
  coords <-
    if (all(c("x", "y") %in% cn)) {
      df[, c("x", "y"), drop = FALSE]
    } else if (all(c("imagecol", "imagerow") %in% cn)) {
      df[, c("imagecol", "imagerow"), drop = FALSE]
    } else {
      # first two numeric columns (skips any 'cell'/barcode id column)
      numeric_cols <- which(vapply(df, is.numeric, logical(1)))
      if (length(numeric_cols) < 2) return(NULL)
      df[, numeric_cols[1:2], drop = FALSE]
    }

  coords <- as.matrix(coords)
  if (!is.numeric(coords)) coords <- apply(coords, 2, as.numeric)
  coords
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
    cn <- colnames(batch)
    isintercept <- if (is.null(cn)) rep(FALSE, ncol(batch)) else grepl("intercept", cn, ignore.case = TRUE)
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
