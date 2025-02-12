test_that("tps function works", {
  expect_error(bs.tps(1:4, 1:4, 0), "greater")
  expect_error(bs.tps(1:4, 1:4, -1), "greater")
  expect_error(bs.tps(1:4, 1:4, 1.2), "integer")
  expect_error(ncol(bs.tps(1:4, 1:6, 3)), "length")
  expect_equal(ncol(bs.tps(1:4, 1:4, 3)), 3 ^ 2)
  expect_equal(ncol(bs.tps(c(1:4, 1:2), 1:6, 3)), 3 * 2)
  expect_equal(bs.tps(1:4, 1:4, 1)[, 1], scale(splines::ns(1:4, 1)^2, scale = FALSE)[, 1])
})

test_that("fitSpaNorm internal function checks works", {
  set.seed(36)
  emat = matrix(rpois(6 * 20, 10), 6, 20)
  coords = cbind(1:20, 1:20)

  # lambda.a
  expect_error(fitSpaNorm(emat, coords, sample.p = 1, gene.model = "nb", lambda.a = 0), "greater")
  expect_error(fitSpaNorm(emat, coords, sample.p = 1, gene.model = "nb", lambda.a = -1), "greater")

  # sample.p
  expect_error(fitSpaNorm(emat, coords, sample.p = 0, gene.model = "nb"), "0,1")
  expect_error(fitSpaNorm(emat, coords, sample.p = 1.1, gene.model = "nb"), "0,1")
  expect_error(fitSpaNorm(emat, coords, sample.p = -1, gene.model = "nb"), "0,1")
  # expect_warning(fitSpaNorm(emat, coords, sample.p = 1), "consider")
  # expect_error(suppressMessages(suppressWarnings(fitSpaNorm(emat, coords, sample.p = 1 / 1e5, gene.model = "nb"))), "too small")
})

test_that("sampleRandom works", {
  coords = cbind(1:6, 1:6)
  expect_type(sampleRandom(coords, 3), "logical")
  expect_length(sampleRandom(coords, 3), 6)
  expect_equal(sum(sampleRandom(coords, 3)), 3)
})

test_that("sampleRandom works", {
  coords = cbind(1:6, 1:6)
  expect_type(sampleRandom(coords, 3), "logical")
  expect_length(sampleRandom(coords, 3), 6)
  expect_equal(sum(sampleRandom(coords, 3)), 3)
})

test_that("getAdjustmentFun works", {
  expect_type(getAdjustmentFun("nb", "auto"), "closure")
  expect_type(getAdjustmentFun("nb", "logpac"), "closure")
  expect_type(getAdjustmentFun("nb", "pearson"), "closure")
  expect_type(getAdjustmentFun("nb", "medbio"), "closure")
  expect_type(getAdjustmentFun("nb", "meanbio"), "closure")
  expect_error(getAdjustmentFun("nb", "abc"), "invalid")

  expect_error(getAdjustmentFun("abc", "abc"), "gene model")
})

test_that("matchDftps handles single and paired df.tps values correctly", {
  # Test single value df1 matches
  expect_true(matchDftps(6, c(6, 4, 3, 2)))
  expect_true(matchDftps(6, c(6, 5, 3, 2)))  # Should match if first two values max is 6
  expect_true(matchDftps(4, c(4, 3, 2, 1)))  # Should match with smaller values
  
  # Test single value df1 non-matches
  expect_false(matchDftps(6, c(5, 4, 3, 2)))  # Max of first pair should match df1
  expect_false(matchDftps(6, c(6, 5, 4, 4)))  # Second pair max should match ceiling(df1/2)
  expect_false(matchDftps(5, c(6, 4, 3, 2)))  # df1 should match max of first pair
  
  # Test paired value df1 matches
  expect_true(matchDftps(c(6, 3), c(6, 4, 3, 2)))
  expect_true(matchDftps(c(6, 4), c(6, 5, 4, 3)))
  expect_true(matchDftps(c(4, 2), c(4, 3, 2, 1)))
  
  # Test paired value df1 non-matches
  expect_false(matchDftps(c(6, 3), c(5, 4, 3, 2)))  # First value should match max of first pair
  expect_false(matchDftps(c(6, 3), c(6, 5, 2, 1)))  # Second value should match max of second pair
  expect_false(matchDftps(c(5, 3), c(6, 4, 3, 2)))  # Values should match respective maxes
  
  # Test edge cases
  expect_error(matchDftps(c(6), c(6)))  # df2 should be length 4
  expect_error(matchDftps(c(6, 3, 1), c(6, 4, 3, 2)))  # df1 should be length 1 or 2
})
