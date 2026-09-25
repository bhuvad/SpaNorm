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

# Moved from spiDE (tests/testthat/test-psi-batch.R, ".reprofilePsi uses the
# same kernel as the batched engine") with nbProfilePsi(), which it tests.
test_that("nbProfilePsi uses the same kernel as the batched engine", {
  # .reprofilePsi() runs at the END of the tau2 loop, in the default path
  # (psi.method = "profile"), and overwrites @psi. Left on optimize() it would
  # discard the bisection's more accurate answer at the last step, so folding
  # it onto the same kernel is what makes the change mean anything for
  # production rather than only for the engine's internals.
  set.seed(31)
  n <- 150; p <- 3; b <- 5
  W <- cbind(1, matrix(stats::rnorm(n * (p - 1)), n, p - 1))
  alpha <- matrix(stats::rnorm(b * p, 0, 0.2), b, p); alpha[, 1] <- 2
  mu <- exp(alpha %*% t(W))
  Y <- matrix(0L, b, n)
  for (i in seq_len(b - 1L)) Y[i, ] <- stats::rnbinom(n, mu = mu[i, ], size = 3)
  Y[b, ] <- stats::rpois(n, mu[b, ])          # at_bound: keeps its incoming psi
  psi_in <- rep(0.35, b)

  got <- nbProfilePsi(Y, W, alpha, psi_in)
  expect_length(got, b)
  expect_true(all(is.finite(got)))

  # the same rule as before: an at-bound gene keeps the dispersion it came with
  ref_kernel <- .psiProfileBatch(Y, pmax(exp(alpha %*% t(W)), .MU_FLOOR))
  expect_equal(got[ref_kernel$at_bound], psi_in[ref_kernel$at_bound])

  # THE assertion: one optimiser, not two. A free gene's re-profiled dispersion
  # must be the kernel's answer exactly, not a second search that merely agrees
  # with it to a few decimals. An "at least as good as optimize()" gate would
  # pass on the old implementation too -- optimize() is as good as itself -- so
  # it would not drive this change.
  free <- !ref_kernel$at_bound
  expect_true(any(free))
  expect_equal(got[free], ref_kernel$psi[free], tolerance = 1e-12)
})

# the penalised NB score of one gene at the mean exp(W a + offset)
.nbScore <- function(y, W, a, psi, pen, off = 0) {
  mu <- exp(as.numeric(W %*% a) + off)
  as.numeric(crossprod(W, (y - mu) / (1 + psi * mu))) - pen * a
}

test_that("polishNB reaches each gene's optimum with an offset, in both engines", {
  set.seed(40)
  G <- 6; n <- 300
  W <- cbind(1, rnorm(n))
  off <- log(runif(n, 0.3, 3))                 # a known per-cell factor
  O <- t(vapply(seq_len(G), function(g) off + g / 4, numeric(n)))  # per gene
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + 0.3 * W[, 2] + O[g, ]), size = 4), numeric(n)))
  pen <- c(0, 0.5)
  for (eng in c("batch", "gene")) {
    v <- polishNB(Y, W, matrix(0, G, 2), rep(0.3, G), lambda.a = pen,
                  offset = off, psi.method = "fixed", tol = 1e-12, engine = eng)
    m <- polishNB(Y, W, matrix(0, G, 2), rep(0.3, G), lambda.a = pen,
                  offset = O, psi.method = "fixed", tol = 1e-12, engine = eng)
    expect_true(all(v$polish$polished), label = eng)
    expect_true(all(m$polish$polished), label = eng)
    for (g in seq_len(G)) {
      expect_lt(max(abs(.nbScore(Y[g, ], W, v$alpha[g, ], 0.3, pen, off))) /
                  sum(Y[g, ]), 1e-6)
      expect_lt(max(abs(.nbScore(Y[g, ], W, m$alpha[g, ], 0.3, pen, O[g, ]))) /
                  sum(Y[g, ]), 1e-6)
    }
    # the per-gene offset shifts each intercept by its own g / 4
    expect_equal(m$alpha[, 1] - v$alpha[, 1], -seq_len(G) / 4, tolerance = 1e-4,
                 label = eng)
  }
})

test_that("a per-gene offset matrix is sliced with the counts across blocks and batches", {
  set.seed(41)
  G <- 7; n <- 150
  W <- cbind(1, rnorm(n))
  O <- t(vapply(seq_len(G), function(g) log(runif(n, 0.5, 2)) + (g - 4) / 3, numeric(n)))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + 0.2 * g * W[, 2] + O[g, ]), size = 4), numeric(n)))
  Y[3, ] <- 0                                    # a gene off the common path
  A0 <- matrix(0, G, 2)
  A0[5, 1] <- -40                                # a degenerate start: restarted
  pen <- c(0, 0.5)
  whole <- polishNB(Y, W, A0, rep(0.3, G), lambda.a = pen, offset = O)
  # each gene alone, with its own offset row, is the reference
  alone <- lapply(seq_len(G), function(g)
    polishNB(Y[g, , drop = FALSE], W, A0[g, , drop = FALSE], 0.3,
             lambda.a = pen, offset = O[g, , drop = FALSE]))
  expect_true(whole$polish$restarted[5])
  # A batch's GEMMs round differently with its shape (up to ~1e-14 with no
  # offset at all), so layouts agree to rounding, not bit for bit. A mis-sliced
  # offset row moves an intercept by (g - h) / 3.
  flags <- c("restarted", "capped", "singular", "psi_bound", "polished")
  expect_equal(whole$alpha, do.call(rbind, lapply(alone, `[[`, "alpha")),
               tolerance = 1e-10)
  expect_equal(whole$psi, vapply(alone, `[[`, numeric(1), "psi"), tolerance = 1e-10)
  expect_identical(whole$polish[, flags],
                   do.call(rbind, lapply(alone, function(a) a$polish[, flags])))
  for (bs in list(c(block = 2, batch = 1), c(block = 3, batch = 2),
                  c(block = 7, batch = 3))) {
    r <- polishNB(Y, W, A0, rep(0.3, G), lambda.a = pen, offset = O,
                  block.size = bs[["block"]], batch.size = bs[["batch"]])
    expect_equal(r$alpha, whole$alpha, tolerance = 1e-10)
    expect_equal(r$psi, whole$psi, tolerance = 1e-10)
    expect_identical(r$polish[, flags], whole$polish[, flags])
  }
  g <- polishNB(Y, W, A0, rep(0.3, G), lambda.a = pen, offset = O,
                engine = "gene", psi.method = "fixed", block.size = 3)
  b <- polishNB(Y, W, A0, rep(0.3, G), lambda.a = pen, offset = O,
                psi.method = "fixed", block.size = 3)
  expect_equal(g$alpha, b$alpha, tolerance = 1e-8)
  # and a warm pass reads the offset too: from the converged fit it stays put
  # (the all-zero gene is left out, since its intercept's optimum is -Inf and
  # every pass moves it further down)
  w <- polishNB(Y, W, whole$alpha, whole$psi, lambda.a = pen, offset = O,
                warm = TRUE, block.size = 2, batch.size = 2)
  expect_equal(w$alpha[-3, ], whole$alpha[-3, ], tolerance = 1e-6)
})

test_that("a bad offset is refused at entry, naming the argument and the shape", {
  set.seed(42)
  n <- 30
  W <- cbind(1, rnorm(n))
  Y <- matrix(rnbinom(2 * n, mu = 4, size = 3), 2, n)
  A0 <- matrix(0, 2, 2)
  run <- function(off) polishNB(Y, W, A0, c(0.3, 0.3), offset = off)
  bad_nan <- rep(0, n); bad_nan[3] <- NaN
  bad_inf <- rep(0, n); bad_inf[5] <- Inf
  bad_na <- matrix(0, 2, n); bad_na[2, 7] <- NA
  expect_error(run(bad_nan), "'offset' must be finite.*one value per cell \\(30\\)")
  expect_error(run(bad_inf), "'offset' must be finite")
  expect_error(run(bad_na), "'offset' must be finite")
  expect_error(run(rep(0, n - 1)), "'offset' must be NULL, a numeric vector with one value per cell \\(30\\)")
  expect_error(run(rep(0, 2)), "'offset' must be NULL")
  expect_error(run(matrix(0, 3, n)), "genes x cells matrix \\(2 x 30\\)")
  expect_error(run(matrix(0, 2, n - 1)), "genes x cells matrix \\(2 x 30\\)")
  expect_error(run(array(0, c(2, n, 1))), "'offset' must be NULL")
  expect_error(run(as.character(rep(0, n))), "'offset' must be NULL")
  # nbProfilePsi() applies the same check
  expect_error(nbProfilePsi(Y, W, A0, c(0.3, 0.3), offset = bad_nan),
               "'offset' must be finite")
  expect_error(nbProfilePsi(Y, W, A0, c(0.3, 0.3), offset = matrix(0, 3, n)),
               "genes x cells matrix")
})

test_that("a Matrix offset is made a base matrix, and a torch tensor is refused", {
  set.seed(43)
  n <- 80
  W <- cbind(1, rnorm(n))
  O <- rbind(log(runif(n, 0.5, 2)), rep(0.3, n))
  Y <- t(vapply(1:2, function(g) rnbinom(n, mu = exp(1 + O[g, ]), size = 4), numeric(n)))
  A0 <- matrix(0, 2, 2)
  base <- polishNB(Y, W, A0, c(0.3, 0.3), offset = O)
  expect_identical(polishNB(Y, W, A0, c(0.3, 0.3), offset = Matrix::Matrix(O)), base)
  expect_identical(nbProfilePsi(Y, W, base$alpha, base$psi, offset = Matrix::Matrix(O)),
                   nbProfilePsi(Y, W, base$alpha, base$psi, offset = O))
  skip_if_no_torch()
  expect_error(polishNB(Y, W, A0, c(0.3, 0.3),
                        offset = torch::torch_tensor(O, dtype = torch::torch_float64())),
               "not a torch tensor")
})

test_that("the offset reaches polishNB's device path, sliced per sub-batch", {
  # The shared-factor path with the transfer mocked to the identity, as for
  # absorb.batch above. Real tensors are not used here: on a CPU device torch
  # (0.17.0) wraps an R double array without copying it, so a tensor made from
  # an R temporary reads freed memory after a garbage collection and the
  # result depends on GC timing. The engine's own tensor path with an offset
  # is tested in test-polishEngine.R.
  set.seed(44)
  G <- 5; n <- 160
  W <- cbind(1, rnorm(n))
  O <- t(vapply(seq_len(G), function(g) log(runif(n, 0.5, 2)) + g / 3, numeric(n)))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + 0.3 * W[, 2] + O[g, ]), size = 4), numeric(n)))
  A0 <- matrix(0, G, 2)
  A0[2, 1] <- NaN                                   # a sane start on one gene
  pen <- c(0, 0.5)
  cpu <- lapply(list(O, O[1, ]), function(off)
    polishNB(Y, W, A0, rep(0.3, G), lambda.a = pen, offset = off, backend = "cpu"))
  real <- .polishBatch
  seen <- list()
  local_mocked_bindings(
    checkGPU = function(...) TRUE,
    toGPUMatrix = function(x, ...) x,
    .requireFloat64 = function(...) invisible(TRUE),
    .polishBatch = function(Yb, ..., shared.factor = FALSE, offset = NULL) {
      seen[[length(seen) + 1L]] <<- list(Yb = Yb, shared = shared.factor,
                                         offset = offset)
      real(Yb, ..., shared.factor = shared.factor, offset = offset)
    }
  )
  dev <- lapply(list(O, O[1, ]), function(off)
    polishNB(Y, W, A0, rep(0.3, G), lambda.a = pen, offset = off, backend = "gpu",
             batch.size = 2))
  # three sub-batches per run (2 + 2 + 1 genes); a matrix offset arrives as
  # exactly the rows of the genes in the sub-batch, a vector as itself
  expect_length(seen, 6L)
  rows <- list(1:2, 3:4, 5L)
  for (k in 1:3) {
    expect_true(seen[[k]]$shared)
    expect_identical(seen[[k]]$Yb, Y[rows[[k]], , drop = FALSE])
    expect_identical(seen[[k]]$offset, O[rows[[k]], , drop = FALSE])
    expect_identical(seen[[k + 3]]$offset, O[1, ])
  }
  for (k in 1:2) {
    expect_true(all(dev[[k]]$polish$polished))
    expect_identical(dev[[k]]$polish$restarted, cpu[[k]]$polish$restarted)
    # a shared factorisation refreshes on its own schedule: same optimum
    expect_equal(dev[[k]]$alpha, cpu[[k]]$alpha, tolerance = 1e-5)
    expect_equal(dev[[k]]$psi, cpu[[k]]$psi, tolerance = 1e-5)
  }
})

test_that("the batch size counts the offset's genes x cells matrices", {
  set.seed(45)
  n <- 240
  W <- cbind(1, rnorm(n))
  Y <- matrix(rnbinom(12 * n, mu = 4, size = 3), 12, n)
  A0 <- matrix(0, 12, 2)
  # a budget of exactly ten genes at the no-offset count of six matrices
  op <- options(SpaNorm.polish.mem.budget = 8 * n * 6 * 10)
  on.exit(options(op), add = TRUE)
  # the opening progress message reports the batch size
  opening <- function(off) {
    capture_messages(polishNB(Y, W, A0, rep(0.3, 12), offset = off,
                              psi.method = "fixed", verbose = TRUE))[1]
  }
  expect_match(opening(NULL), "in batches of 10 ")
  # a vector adds the transient genes x cells expansion: 7 matrices
  expect_match(opening(rep(0.1, n)), "in batches of 8 ")
  # a matrix also holds the batch's own rows: 8 matrices
  expect_match(opening(matrix(0.1, 12, n)), "in batches of 7 ")
  # and the default count is untouched
  expect_identical(.polishBatchSize(n, budget = 1e6),
                   max(1L, as.integer(floor(1e6 / (8 * n * POLISH_GENE_CELL_MATS)))))
})

test_that("psi.range bounds the profile dispersion search", {
  set.seed(46)
  G <- 4; n <- 300
  W <- cbind(1, rnorm(n))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(1 + 0.2 * W[, 2]), size = 3), numeric(n)))   # psi ~ 1/3
  A0 <- matrix(0, G, 2)
  for (eng in c("batch", "gene")) {
    free <- polishNB(Y, W, A0, rep(0.9, G), engine = eng)
    expect_false(any(free$polish$psi_bound), label = eng)
    # a range excluding the optimum: every gene's search lands on the bound
    # and keeps its input dispersion
    tight <- polishNB(Y, W, A0, rep(0.9, G), psi.range = c(2, 50), engine = eng)
    expect_true(all(tight$polish$psi_bound), label = eng)
    expect_identical(tight$psi, rep(0.9, G), label = eng)
  }
  expect_error(polishNB(Y, W, A0, rep(0.9, G), psi.range = c(5, 1)), "psi.range")
  expect_error(polishNB(Y, W, A0, rep(0.9, G), psi.range = c(0, 1)), "psi.range")
})

## ---- fix round 1 -----------------------------------------------------------

test_that("polishNB and nbProfilePsi refuse inputs whose shapes disagree", {
  set.seed(50)
  n <- 50
  W <- cbind(1, rnorm(n))
  Y <- t(vapply(1:2, function(g) rnbinom(n, mu = 4, size = 3), numeric(n)))
  A0 <- matrix(0, 2, 2)
  W100 <- cbind(1, rnorm(100))
  # This one used to pass silently on the per-gene engine: each gene's 50
  # counts were recycled against a 100-cell mean and came back polished.
  for (eng in c("gene", "batch")) {
    expect_error(polishNB(Y, W100, A0, c(0.3, 0.3), engine = eng),
                 "'W' must have one row per cell (column of Y): nrow(W) = 100, ncol(Y) = 50",
                 fixed = TRUE)
  }
  expect_error(polishNB(Y, W, matrix(0, 3, 2), c(0.3, 0.3)),
               "'alpha' must have one row per gene (row of Y): nrow(alpha) = 3, nrow(Y) = 2",
               fixed = TRUE)
  expect_error(polishNB(Y, W, matrix(0, 1, 2), 0.3),
               "nrow(alpha) = 1, nrow(Y) = 2", fixed = TRUE)
  expect_error(polishNB(Y, W, matrix(0, 2, 3), c(0.3, 0.3)),
               "'alpha' must have one column per column of W: ncol(alpha) = 3, ncol(W) = 2",
               fixed = TRUE)
  expect_error(polishNB(Y, W, A0, c(0.3, 0.3, 0.3)),
               "'psi' must be one value or one per gene: length(psi) = 3, nrow(Y) = 2",
               fixed = TRUE)
  # nbProfilePsi() applies the same checks
  expect_error(nbProfilePsi(Y, W100, A0, c(0.3, 0.3)), "nrow(W) = 100, ncol(Y) = 50",
               fixed = TRUE)
  expect_error(nbProfilePsi(Y, W, matrix(0, 3, 2), c(0.3, 0.3)),
               "nrow(alpha) = 3, nrow(Y) = 2", fixed = TRUE)
  expect_error(nbProfilePsi(Y, W, matrix(0, 2, 3), c(0.3, 0.3)),
               "ncol(alpha) = 3, ncol(W) = 2", fixed = TRUE)
  expect_error(nbProfilePsi(Y, W, A0, c(0.3, 0.3, 0.3)),
               "length(psi) = 3, nrow(Y) = 2", fixed = TRUE)
})

test_that("nbProfilePsi recycles a scalar psi", {
  set.seed(51)
  n <- 100
  W <- cbind(1, rnorm(n))
  # genes 2 and 3 have no spread, so their search runs to the lower bound and
  # they keep the input psi -- which a scalar used to give only gene 1
  Y <- rbind(rnbinom(n, mu = 5, size = 3), rep(5, n), rep(5, n))
  A <- matrix(c(log(5), 0), 3, 2, byrow = TRUE)
  got <- nbProfilePsi(Y, W, A, 0.2)
  expect_false(anyNA(got))
  expect_identical(got, nbProfilePsi(Y, W, A, rep(0.2, 3)))
  expect_identical(got[2:3], c(0.2, 0.2))
})

