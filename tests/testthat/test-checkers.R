test_that("checkSPE validates all count values, not just the first", {
  set.seed(1)
  counts <- matrix(rpois(5 * 4, 5), nrow = 5, ncol = 4)
  coords <- cbind(x = seq_len(4), y = rev(seq_len(4)))
  make_spe <- function(m) {
    SpatialExperiment::SpatialExperiment(
      assays = list(counts = m),
      spatialCoords = coords
    )
  }

  # valid counts pass
  spe_ok <- make_spe(counts)
  expect_silent(checkSPE(spe_ok))

  # a negative value away from position 1 must be rejected
  neg <- counts
  neg[3, 3] <- -5L
  expect_error(checkSPE(make_spe(neg)), "negative")

  # a non-integer value away from position 1 must be rejected
  frac <- counts
  frac[4, 2] <- 1.5
  expect_error(checkSPE(make_spe(frac)), "non-integer")

  # the same must hold for sparse counts (only stored values are inspected)
  sp <- Matrix::Matrix(counts, sparse = TRUE)
  sp[2, 4] <- -3
  expect_error(checkSPE(make_spe(sp)), "negative")
})

# Build a Seurat object carrying a VisiumV2-style FOV image (no @coordinates
# slot), with the image cells deliberately shuffled to exercise alignment.
make_fov_seurat <- function(ng = 8, ns = 15) {
  set.seed(1)
  counts <- matrix(rpois(ng * ns, 5), ng, ns,
                   dimnames = list(paste0("g", seq_len(ng)), paste0("c", seq_len(ns))))
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts, assay = "Spatial"))
  df <- data.frame(x = runif(ns, 0, 100), y = runif(ns, 0, 100), cell = colnames(counts))
  fov <- suppressWarnings(
    SeuratObject::CreateFOV(df[sample(ns), ], type = "centroids", key = "fov_", assay = "Spatial"))
  suppressWarnings(obj[["slice1"]] <- fov) # assignment warns on unordered cells
  obj
}

test_that("extractSeuratCoords works on VisiumV2/FOV images (GitHub #18)", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  obj <- make_fov_seurat()
  # the image class has no @coordinates slot -- the old direct slot access threw
  expect_false("coordinates" %in% methods::slotNames(obj@images[[1]]))

  coords <- extractSeuratCoords(obj)
  expect_true(is.matrix(coords))
  expect_equal(dim(coords), c(ncol(obj), 2L))
  expect_true(all(is.finite(coords)))
  # rows must align to the counts cell order despite the shuffled image
  tc <- SeuratObject::GetTissueCoordinates(obj)
  rownames(tc) <- as.character(tc$cell)
  expect_equal(unname(coords), unname(as.matrix(tc[colnames(obj), c("x", "y")])))

  expect_no_error(checkSeurat(obj))
})

test_that("checkSeurat honours the assay argument (GitHub #19)", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  obj <- make_fov_seurat()
  # a transformed default assay; raw counts remain in 'Spatial'
  obj[["SCT"]] <- suppressWarnings(SeuratObject::CreateAssayObject(
    counts = SeuratObject::LayerData(obj, layer = "counts", assay = "Spatial")))
  SeuratObject::DefaultAssay(obj) <- "SCT"

  expect_no_error(checkSeurat(obj, assay = "Spatial"))
  expect_error(checkSeurat(obj, assay = "NotAnAssay"), "not found")

  # a non-integer assay is rejected when explicitly selected
  bad <- SeuratObject::LayerData(obj, layer = "counts", assay = "Spatial")
  bad[1, 1] <- 0.5
  obj[["Frac"]] <- suppressWarnings(SeuratObject::CreateAssayObject(counts = bad))
  expect_error(checkSeurat(obj, assay = "Frac"), "non-integer")
})
