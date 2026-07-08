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
