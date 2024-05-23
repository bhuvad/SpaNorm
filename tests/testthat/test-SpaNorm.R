test_that("tps function works", {
  expect_error(bs.tps(1:4, 1:4, 0), "greater")
  expect_error(bs.tps(1:4, 1:4, -1), "greater")
  expect_error(bs.tps(1:4, 1:4, 1.2), "integer")
  expect_equal(ncol(bs.tps(1:4, 1:4, 3)), 3 ^ 2)
  expect_equal(bs.tps(1:4, 1:4, 1)[, 1], scale(splines::ns(1:4, 1)^2, scale = FALSE)[, 1])
})

test_that("SpaNorm checks works", {
  set.seed(36)
  emat = matrix(rpois(6 * 4e3, 10), 6, 4e3)
  coords = cbind(1:4e3, 1:4e3)

  # step.factor
  expect_error(SpaNorm_intl(emat, coords, step.factor = 0), "0,1")
  expect_error(SpaNorm_intl(emat, coords, step.factor = 1), "0,1")
  expect_error(SpaNorm_intl(emat, coords, step.factor = 1.1), "0,1")
  expect_error(SpaNorm_intl(emat, coords, step.factor = -1.1), "0,1")

  # lambda.a
  expect_error(SpaNorm_intl(emat, coords, lambda.a = 0), "greater")
  expect_error(SpaNorm_intl(emat, coords, lambda.a = -1), "greater")

  # maxit
  expect_error(SpaNorm_intl(emat, coords, inner.maxit = 0), "greater")
  expect_error(SpaNorm_intl(emat, coords, inner.maxit = -1), "greater")
  expect_error(SpaNorm_intl(emat, coords, inner.maxit = 1.1), "integer")
  expect_error(SpaNorm_intl(emat, coords, outer.maxit = 0), "greater")
  expect_error(SpaNorm_intl(emat, coords, outer.maxit = -1), "greater")
  expect_error(SpaNorm_intl(emat, coords, outer.maxit = 1.1), "integer")

  # sample.p
  expect_error(SpaNorm_intl(emat, coords, sample.p = 0), "0,1")
  expect_error(SpaNorm_intl(emat, coords, sample.p = 1.1), "0,1")
  expect_error(SpaNorm_intl(emat, coords, sample.p = -1), "0,1")
  # expect_warning(SpaNorm_intl(emat, coords, sample.p = 1), "consider")
  expect_error(suppressWarnings(SpaNorm_intl(emat, coords, sample.p = 1 / 1e5)), "too small")
})
