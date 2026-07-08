test_that("devianceResiduals matches a dense base-R reference", {
  set.seed(7)
  ng <- 8; nc <- 12; ncov <- 3
  W <- matrix(rnorm(nc * ncov), nc, ncov)
  alpha <- matrix(rnorm(ng * ncov) * 0.1, ng, ncov)
  gmean <- rnorm(ng)
  psi <- runif(ng, 0.1, 0.5)
  Ycounts <- matrix(rpois(ng * nc, 5), ng, nc)
  Ycounts[sample.int(length(Ycounts), 10)] <- 0 # ensure some structural zeros
  Y <- Matrix::Matrix(Ycounts, sparse = TRUE)

  ft <- list(gmean = gmean, alpha = alpha, W = W, psi = psi)

  # dense base-R reference for NB deviance residuals (independent of the
  # sparse/dense representation used internally)
  mu <- calculateMu(ft$gmean, ft$alpha, ft$W)
  pp <- ft$psi
  pp.max <- exp(median(log(pp)) + 3 * mad(log(pp)))
  pp <- pmin(pp, pp.max)
  Yd <- as.matrix(Y)
  d <- Yd * log(Yd / mu)
  d[is.nan(d)] <- 0
  d <- d - (Yd + 1 / pp) * log((1 + Yd * pp) / (1 + mu * pp))
  d <- pmax(d, 0)
  ref <- sign(Yd - mu) * sqrt(2 * d)

  got <- as.matrix(devianceResiduals(Y, ft))
  expect_equal(got, ref, tolerance = 1e-8)
  expect_false(any(is.na(got)))
})
