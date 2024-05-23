test_that("check SpatialExperiment works", {
  library(SpatialExperiment)
  data(HumanDLPFC)
  expect_no_error(checkSPE(HumanDLPFC))
  
  # check non integer
  counts(HumanDLPFC)[1, 1] = 0.1
  expect_error(checkSPE(HumanDLPFC), "non-integer")

  counts(HumanDLPFC)[1, 1] = -1
  expect_error(checkSPE(HumanDLPFC), "negative")

  expect_error(checkSPE(HumanDLPFC[0, 0]), "nrow.+ncol")
  expect_error(checkSPE(HumanDLPFC[0, ]), "nrow.+ncol")
  expect_error(checkSPE(HumanDLPFC[, 0]), "nrow.+ncol")
})

test_that("check Seurat works", {
  library(Seurat)
  data(HumanDLPFC)
  HumanDLPFC = suppressWarnings(Seurat::as.Seurat(HumanDLPFC))

  expect_no_error(checkSeurat(HumanDLPFC))
  
  # check non integer
  LayerData(HumanDLPFC, layer = "counts")[1, 1] = 0.1
  expect_error(checkSeurat(HumanDLPFC), "non-integer")

  LayerData(HumanDLPFC, layer = "counts")[1, 1] = -1
  expect_error(checkSeurat(HumanDLPFC), "negative")

  expect_error(checkSeurat(HumanDLPFC[0, 0]), "No cells")
  # expect_error(checkSeurat(HumanDLPFC[0, ]), "nrow.+ncol")
  expect_error(checkSeurat(HumanDLPFC[, 0]), "No cells")
})

test_that("check filterGenes works", {
  library(SpatialExperiment)
  library(Seurat)
  data(HumanDLPFC)
  
  # check out of range prop - SpatialExperiment
  expect_error(filterGenes(HumanDLPFC[, 0]), "ncol") # ensure checkSPE is used
  expect_error(filterGenes(HumanDLPFC, 1.1), "prop")
  expect_error(filterGenes(HumanDLPFC, -0.1), "prop")
  expect_type(filterGenes(HumanDLPFC), "logical")
  expect_equal(length(filterGenes(HumanDLPFC)), nrow(HumanDLPFC))
  expect_equal(sum(filterGenes(HumanDLPFC, 0)), nrow(HumanDLPFC))
  expect_equal(sum(filterGenes(HumanDLPFC, 1)), 5) # sum(apply(counts(HumanDLPFC) > 0, 1, \(x) all(x)))
  
  # check out of range prop - Seurat
  HumanDLPFC = suppressWarnings(Seurat::as.Seurat(HumanDLPFC))
  expect_error(filterGenes(HumanDLPFC[, 0]), "No cells") # ensure checkSPE is used
  expect_error(filterGenes(HumanDLPFC, 1.1), "prop")
  expect_error(filterGenes(HumanDLPFC, -0.1), "prop")
  expect_type(filterGenes(HumanDLPFC), "logical")
  expect_equal(length(filterGenes(HumanDLPFC)), nrow(HumanDLPFC))
  expect_equal(sum(filterGenes(HumanDLPFC, 0)), nrow(HumanDLPFC))
  expect_equal(sum(filterGenes(HumanDLPFC, 1)), 5) # sum(apply(counts(HumanDLPFC) > 0, 1, \(x) all(x)))
})
