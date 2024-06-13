test_that("fastSizeFactors works", {
  suppressPackageStartupMessages(library(SpatialExperiment))
  data(HumanDLPFC)

  expect_error(filterGenes(HumanDLPFC[, 0]), "ncol") # ensure checkSPE is used
  expect_equal(length(fastSizeFactors(HumanDLPFC)$sizeFactor), ncol(HumanDLPFC))
})
