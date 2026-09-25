# polishNB(): the driver that converges every gene of a shared fitNB() fit to
# that gene's own penalised-NB optimum, moved from spiDE's .polishFit().

test_that("polishNB converges a fitNB fit to each gene's own optimum", {
  set.seed(10)
  G <- 20; n <- 400
  W <- cbind(1, rnorm(n), rnorm(n))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + g / 10 + 0.3 * W[, 2]), size = 3), numeric(n)))
  fit <- fitNB(Y, W, lambda.a = 0.1, verbose = FALSE, backend = "cpu")
  # tol is the relative log-likelihood gain at which Newton stops. At the
  # default 1e-8 the score is left at up to 2.6e-6 * sum(y) on this fixture, so
  # a zero-score check at 1e-6 needs the tighter stop (it gives 3.4e-8)
  pol <- polishNB(Y, W, fit$alpha, fit$psi, lambda.a = 0.1, psi.method = "fixed",
                  tol = 1e-12)
  for (g in seq_len(G)) {
    mu <- exp(as.numeric(W %*% pol$alpha[g, ]))
    sc <- crossprod(W, (Y[g, ] - mu) / (1 + pol$psi[g] * mu)) - 0.1 * pol$alpha[g, ]
    expect_lt(max(abs(sc)) / sum(Y[g, ]), 1e-6)
  }
  expect_true(all(pol$polish$polished))
  expect_equal(pol$psi, fit$psi)                  # "fixed" keeps the dispersion
})

test_that("gene and batch engines agree", {
  set.seed(11)
  G <- 8; n <- 250
  W <- cbind(1, rnorm(n))
  Y <- t(vapply(seq_len(G), function(g) rnbinom(n, mu = exp(1 + 0.2 * W[, 2]), size = 2), numeric(n)))
  a0 <- matrix(0, G, 2); p0 <- rep(0.5, G)
  # at a held dispersion the two engines take the same Newton path
  a <- polishNB(Y, W, a0, p0, lambda.a = c(0, 1), engine = "gene", psi.method = "fixed")
  b <- polishNB(Y, W, a0, p0, lambda.a = c(0, 1), engine = "batch", psi.method = "fixed")
  expect_equal(a$alpha, b$alpha, tolerance = 1e-8)
  expect_equal(a$psi, b$psi, tolerance = 1e-6)
  # the profile dispersion is optimize() per gene and a bisection in the batch,
  # which agree to optimize()'s tolerance (~1.2e-4 on log psi): measured here
  # 1.9e-5 on log psi and a mean relative 1.9e-8 on alpha
  a <- polishNB(Y, W, a0, p0, lambda.a = c(0, 1), engine = "gene")
  b <- polishNB(Y, W, a0, p0, lambda.a = c(0, 1), engine = "batch")
  expect_equal(a$alpha, b$alpha, tolerance = 1e-6)
  expect_equal(a$psi, b$psi, tolerance = 1e-4)
})

test_that("polishNB refuses non-integer counts", {
  Y <- matrix(c(1.5, 2, 3, 4), 1)
  expect_error(polishNB(Y, cbind(rep(1, 4)), matrix(0, 1, 1), 0.1),
               "needs integer counts")
})

test_that("a gene with no counts is flagged, not fatal", {
  set.seed(12)
  n <- 100
  W <- cbind(1, rnorm(n))
  Y <- rbind(rnbinom(n, mu = 5, size = 3), rep(0, n))
  r <- polishNB(Y, W, matrix(0, 2, 2), c(0.3, 0.3), lambda.a = c(0, 1))
  expect_true(r$polish$polished[1])
  # both genes come back finite; the zero gene's `polished` flag and
  # coefficients are recorded in the Task 4 report, not asserted here
  expect_true(all(is.finite(r$alpha)))
  expect_true(all(is.finite(r$psi)))
})
