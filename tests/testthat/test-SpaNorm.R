test_that("tps function works", {
  expect_error(bs.tps(1:4, 1:4, 0), "greater")
  expect_error(bs.tps(1:4, 1:4, -1), "greater")
  expect_error(bs.tps(1:4, 1:4, 1.2), "integer")
  expect_equal(ncol(bs.tps(1:4, 1:4, 3)), 3 ^ 2)
  expect_equal(bs.tps(1:4, 1:4, 1)[, 1], scale(splines::ns(1:4, 1)^2, scale = FALSE)[, 1])
})

test_that("fitSpaNorm internal function checks works", {
  set.seed(36)
  emat = matrix(rpois(6 * 4e3, 10), 6, 4e3)
  coords = cbind(1:4e3, 1:4e3)

  # lambda.a
  expect_error(fitSpaNorm(emat, coords, sample.p = 1, gene.model = "nb", lambda.a = 0), "greater")
  expect_error(fitSpaNorm(emat, coords, sample.p = 1, gene.model = "nb", lambda.a = -1), "greater")

  # sample.p
  expect_error(fitSpaNorm(emat, coords, sample.p = 0, gene.model = "nb"), "0,1")
  expect_error(fitSpaNorm(emat, coords, sample.p = 1.1, gene.model = "nb"), "0,1")
  expect_error(fitSpaNorm(emat, coords, sample.p = -1, gene.model = "nb"), "0,1")
  # expect_warning(fitSpaNorm(emat, coords, sample.p = 1), "consider")
  expect_error(suppressMessages(suppressWarnings(fitSpaNorm(emat, coords, sample.p = 1 / 1e5, gene.model = "nb"))), "too small")
})

test_that("fitSpaNormNB internal function checks works", {
  set.seed(36)
  emat = matrix(rpois(6 * 4e3, 10), 6, 4e3)
  coords = cbind(1:4e3, 1:4e3)

  # maxit
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = 0), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = -1), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = 1.1), "integer")
})

test_that("fitNBGivenPsi internal function checks works", {
  set.seed(36)
  Y = matrix(0, 6, 10)
  gmean = psi = rep(0, 6)
  W = matrix(0, 10, 2)
  alpha = matrix(0, 6, 2)

  # maxit
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = 0), "greater")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = -1), "greater")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = 1.1), "integer")

  # step.factor
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 0), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 1), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 1.1), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = -1.1), "0,1")
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
