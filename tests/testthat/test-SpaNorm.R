test_that("tps function works", {
  expect_error(bs.tps(1:4, 1:4, 0), "greater")
  expect_error(bs.tps(1:4, 1:4, -1), "greater")
  expect_error(bs.tps(1:4, 1:4, 1.2), "integer")
  expect_equal(ncol(bs.tps(1:4, 1:4, 3)), 3 ^ 2)
  expect_equal(bs.tps(1:4, 1:4, 1)[, 1], scale(splines::ns(1:4, 1)^2, scale = FALSE)[, 1])
})

test_that("SpaNorm internal function checks works", {
  set.seed(36)
  emat = matrix(rpois(6 * 4e3, 10), 6, 4e3)
  coords = cbind(1:4e3, 1:4e3)

  # lambda.a
  expect_error(fitSpaNorm(emat, coords, "nb", lambda.a = 0), "greater")
  expect_error(fitSpaNorm(emat, coords, "nb", lambda.a = -1), "greater")

  # sample.p
  expect_error(fitSpaNorm(emat, coords, sample.p = 0), "0,1")
  expect_error(fitSpaNorm(emat, coords, sample.p = 1.1), "0,1")
  expect_error(fitSpaNorm(emat, coords, sample.p = -1), "0,1")
  # expect_warning(fitSpaNorm(emat, coords, sample.p = 1), "consider")
  expect_error(suppressMessages(suppressWarnings(fitSpaNorm(emat, coords, sample.p = 1 / 1e5))), "too small")

  # step.factor
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, step.factor = 0), "0,1")
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, step.factor = 1), "0,1")
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, step.factor = 1.1), "0,1")
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, step.factor = -1.1), "0,1")

  # maxit
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, inner.maxit = 0), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, inner.maxit = -1), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, inner.maxit = 1.1), "integer")
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, outer.maxit = 0), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, outer.maxit = -1), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), 0.0001, outer.maxit = 1.1), "integer")
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
