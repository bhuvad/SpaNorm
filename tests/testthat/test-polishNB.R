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

# A random-slope design in the shape spiDE hands the polish: cell-type
# intercepts, a covariate, nested (sample x cell type) indicators and a slope
# per sample. `group` is the per-sample grouping of the whole random block (the
# per-gene solver's absorption); `nested` marks the indicators alone, the 1x1
# blocks the shared-factor batched solver can absorb.
.slopeDesign <- function(G = 5, seed = 20) {
  set.seed(seed)
  n <- 240
  smp <- rep(1:4, each = 60)
  ct <- rep(c("A", "B"), times = 120)
  x <- rnorm(n)
  Zn <- do.call(cbind, lapply(1:4, function(s)
    cbind(as.numeric(smp == s & ct == "A"), as.numeric(smp == s & ct == "B"))))
  Zs <- do.call(cbind, lapply(1:4, function(s) x * (smp == s)))
  W <- cbind(A = as.numeric(ct == "A"), B = as.numeric(ct == "B"), x = x, Zn, Zs)
  p <- ncol(W)
  u <- rnorm(8, 0, 0.3)
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + 0.2 * g * (ct == "A") + 0.3 * x +
                          u[(smp - 1) * 2 + (ct == "B") + 1]), size = 4), numeric(n)))
  list(Y = Y, W = W, A0 = matrix(0, G, p), psi = rep(0.3, G),
       pen = c(0, 0, 0.01, rep(2, 8), rep(5, 4)),
       start = c(TRUE, TRUE, rep(FALSE, p - 2)),
       nested = c(rep(FALSE, 3), rep(TRUE, 8), rep(FALSE, 4)),
       group = c(rep(NA, 3), rep(1:4, each = 2), 1:4))
}

test_that("absorb.batch reaches the shared-factor batched solver, and a grouping alone goes dense", {
  d <- .slopeDesign()
  # the shared-factor solver absorbs 1x1 blocks only, which is why the
  # grouping cannot be handed to it
  expect_error(.newtonSolverBatch(d$W, d$pen, d$group), "1x1 blocks only")
  cpu <- polishNB(d$Y, d$W, d$A0, d$psi, lambda.a = d$pen, absorb = d$group,
                  start.cols = d$start, backend = "cpu")
  expect_true(all(cpu$polish$polished))

  # Reach the device path on the CPU: report a GPU and make the transfer the
  # identity (every batched kernel runs on a base matrix too), and record what
  # .polishBatch() is handed.
  real <- .polishBatch
  seen <- list()
  local_mocked_bindings(
    checkGPU = function(...) TRUE,
    toGPUMatrix = function(x, ...) x,
    .requireFloat64 = function(...) invisible(TRUE),
    .polishBatch = function(..., shared.factor = FALSE, nested = NULL) {
      seen[[length(seen) + 1L]] <<- list(shared = shared.factor, nested = nested)
      real(..., shared.factor = shared.factor, nested = nested)
    }
  )
  dev <- polishNB(d$Y, d$W, d$A0, d$psi, lambda.a = d$pen, absorb = d$group,
                  absorb.batch = d$nested, start.cols = d$start, backend = "gpu")
  expect_gt(length(seen), 0L)
  for (s in seen) {
    expect_true(s$shared)
    expect_identical(s$nested, d$nested)
  }
  expect_true(all(dev$polish$polished))
  # a shared factorisation refreshes on a different schedule from the per-gene
  # one, so the two reach the same optimum rather than along one path
  expect_equal(dev$alpha, cpu$alpha, tolerance = 1e-5)
  expect_equal(dev$psi, cpu$psi, tolerance = 1e-5)

  # without absorb.batch the grouping falls back to the dense batched solver
  seen <- list()
  dns <- polishNB(d$Y, d$W, d$A0, d$psi, lambda.a = d$pen, absorb = d$group,
                  start.cols = d$start, backend = "gpu")
  expect_gt(length(seen), 0L)
  for (s in seen) {
    expect_true(s$shared)
    expect_identical(s$nested, rep(FALSE, ncol(d$W)))
  }
  expect_true(all(dns$polish$polished))
  expect_equal(dns$alpha, cpu$alpha, tolerance = 1e-5)
  expect_equal(dns$psi, cpu$psi, tolerance = 1e-5)

  # and a logical absorb with no absorb.batch is passed through unchanged
  seen <- list()
  polishNB(d$Y, d$W, d$A0, d$psi, lambda.a = d$pen, absorb = d$nested,
           start.cols = d$start, backend = "gpu")
  for (s in seen) expect_identical(s$nested, d$nested)
})

test_that("absorb.batch must be a logical over the columns of W", {
  d <- .slopeDesign(G = 2)
  expect_error(polishNB(d$Y, d$W, d$A0, d$psi, lambda.a = d$pen,
                        absorb = d$group, absorb.batch = d$group),
               "absorb.batch")
  expect_error(polishNB(d$Y, d$W, d$A0, d$psi, lambda.a = d$pen,
                        absorb = d$group, absorb.batch = d$nested[-1]),
               "absorb.batch")
})

# the profile-ML dispersion of one gene at a fixed mean: the oracle for
# nbProfilePsi(). optimize() minimises the objective, which is flat at its
# optimum, so it places the argmin only to ~sqrt(.Machine$double.eps) relative
# (measured 1e-8 to 4e-8 here) whatever `tol` asks for; nbProfilePsi() bisects
# the score to ~1e-14. Hence 1e-6 against it below.
.mlPsi <- function(y, mu, range = c(1e-3, 1e3)) {
  exp(stats::optimize(function(lp) {
    -sum(stats::dnbinom(y, size = exp(-lp), mu = mu, log = TRUE))
  }, log(range), tol = 1e-12)$minimum)
}

test_that("nbProfilePsi returns the profile-ML dispersion, keeping the input on a bound", {
  set.seed(30)
  n <- 400
  W <- cbind(1, rnorm(n))
  A <- rbind(c(1, 0.3), c(0.5, -0.2), c(log(5), 0))
  Y <- rbind(rnbinom(n, mu = exp(W %*% A[1, ]), size = 2),
             rnbinom(n, mu = exp(W %*% A[2, ]), size = 5),
             rep(5, n))           # no spread at all: the ML runs to the lower bound
  psi_in <- c(0.7, 0.7, 0.123)
  got <- nbProfilePsi(Y, W, A, psi_in)
  for (g in 1:2) {
    expect_equal(got[g], .mlPsi(Y[g, ], exp(as.numeric(W %*% A[g, ]))),
                 tolerance = 1e-6)
  }
  expect_identical(got[3], 0.123)
  # gene blocking is exact
  expect_identical(nbProfilePsi(Y, W, A, psi_in, block.size = 1), got)
})

test_that("nbProfilePsi profiles at the mean the offset gives", {
  set.seed(31)
  n <- 400
  W <- cbind(1, rnorm(n))
  off <- log(runif(n, 0.3, 3))
  A <- rbind(c(1, 0.3), c(0.2, 0.5))
  Y <- t(vapply(1:2, function(g)
    rnbinom(n, mu = exp(as.numeric(W %*% A[g, ]) + off), size = 3), numeric(n)))
  with <- nbProfilePsi(Y, W, A, c(0.5, 0.5), offset = off)
  without <- nbProfilePsi(Y, W, A, c(0.5, 0.5))
  for (g in 1:2) {
    expect_equal(with[g], .mlPsi(Y[g, ], exp(as.numeric(W %*% A[g, ]) + off)),
                 tolerance = 1e-6)
  }
  expect_gt(min(abs(log(with) - log(without))), 0.01)
  # a genes x cells offset gives each gene its own row, also across blocks
  O <- rbind(off, off + 0.5)
  m <- nbProfilePsi(Y, W, A, c(0.5, 0.5), offset = O, block.size = 1)
  expect_equal(m[1], with[1])
  expect_equal(m[2], .mlPsi(Y[2, ], exp(as.numeric(W %*% A[2, ]) + O[2, ])),
               tolerance = 1e-6)
})
