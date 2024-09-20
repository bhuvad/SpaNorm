test_that("batch checking works", {
  batch = rep(LETTERS[1:3], each = 2)

  # NULL case
  expect_equal(checkBatch(NULL, 5), c())
  # char/factor vector
  expect_error(checkBatch(batch, 5), "length")
  expect_equal(nrow(checkBatch(batch, 6)), 6)
  expect_equal(ncol(checkBatch(batch, 6)), 2)
  # numeric vector
  expect_equal(nrow(checkBatch(1:6, 6)), 6)
  expect_equal(ncol(checkBatch(1:6, 6)), 1)
  # check numeric factor
  expect_equal(ncol(checkBatch(as.factor(rep(1:3, 2)), 6)), 2)
  # check NAs
  batch[1] = NA
  expect_error(ncol(checkBatch(batch, 6)), "missing")

  # character matrix
  batch = rep(LETTERS[1:3], each = 2)
  expect_error(checkBatch(matrix(batch), 6), "numeric matrix")
  # numeric matrix
  batch = model.matrix(~batch)
  expect_error(checkBatch(batch, 5), "number of rows")
  # numeric matrix with intercept
  expect_warning(checkBatch(batch, 6), "intercept")
  colnames(batch)[1] = "column1"
  expect_warning(checkBatch(batch, 6), "intercept")
  # check NAs
  batch[1, 1] = NA
  expect_error(ncol(checkBatch(batch, 6)), "missing")
  
  batch = batch[, -1]
  # numeric matrix without intercept
  expect_equal(nrow(checkBatch(batch, 6)), 6)
  expect_equal(ncol(checkBatch(batch, 6)), 2)
})
