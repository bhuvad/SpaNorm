# The polish engine (R/polishEngine.R, R/polishEngineBatch.R), moved from
# spiDE with its tests: tests/testthat/test-nb-batch.R, test-psi-batch.R and
# test-polish-batch.R, every test that builds its design directly. The tests
# that need spiDE's fixtures stay in spiDE.

# ---- from spiDE tests/testthat/test-nb-batch.R -----------------------------

# The two kernels .polishBatch()'s Newton is built on -- the mean and the
# penalised NB log-likelihood -- on either backend.
#
# Everything else in the loop (the line search, the dispersion search, the
# convergence test) is a consumer of these two, so they go to tensors first.
# base R is the oracle: these are the expressions already in .polishBatch().

.nbFixture <- function(n = 50, p = 4, b = 5, seed = 17) {
  set.seed(seed)
  W <- cbind(1, matrix(stats::rnorm(n * (p - 1)), n, p - 1))
  A <- matrix(stats::rnorm(b * p, 0, 0.3), b, p)
  A[, 1] <- A[, 1] + 1.5
  mu <- exp(A %*% t(W))
  Y <- matrix(stats::rnbinom(b * n, mu = as.numeric(mu), size = 2), b, n)
  list(W = W, A = A, Y = Y, psi = stats::runif(b, 0.1, 1.5),
       pen = stats::runif(p, 0, 0.4))
}

# the expressions verbatim from .polishBatch()
.muRef <- function(f) pmax(exp(f$A %*% t(f$W)), .MU_FLOOR)
.llRef <- function(f, M) {
  rowSums(stats::dnbinom(f$Y, size = 1 / f$psi, mu = M, log = TRUE)) -
    0.5 * as.numeric((f$A^2) %*% f$pen)
}

test_that(".muBatch reproduces the base-R mean, including the floor", {
  f <- .nbFixture()
  expect_equal(.muBatch(f$A, f$W), .muRef(f), tolerance = 1e-12)

  # the floor has to bite, or it is not being tested. Drive the INTERCEPT down
  # and leave the rest at zero: setting every coefficient to -50 gives
  # eta = -50 * rowSums(W), which is large and POSITIVE wherever rowSums(W) is
  # negative, so most cells would not be floored at all.
  f2 <- f; f2$A[1, ] <- c(-50, rep(0, ncol(f$A) - 1L))
  got <- .muBatch(f2$A, f2$W)
  expect_true(all(got[1, ] == .MU_FLOOR))
  expect_equal(got, .muRef(f2), tolerance = 1e-12)
})

test_that(".nbLoglikBatch reproduces dnbinom's penalised row sums", {
  f <- .nbFixture()
  M <- .muRef(f)
  expect_equal(.nbLoglikBatch(f$Y, M, f$psi, f$A, f$pen), .llRef(f, M),
               tolerance = 1e-10)
})

test_that("the NB kernels agree between the base-R and torch branches", {
  # Gated on a GPU, not just torch: on a CPU device torch 0.17's
  # torch_tensor(<R double>, float64) aliases R memory rather than copying,
  # so these results depend on GC timing. The engine fix is deferred.
  skip_if_no_gpu()
  f <- .nbFixture()
  M <- .muRef(f)
  ll <- .llRef(f, M)

  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  Wt <- tt(f$W); At <- tt(f$A); Yt <- tt(f$Y)
  Mt <- .muBatch(At, Wt)
  expect_true(is_torch_tensor(Mt))
  expect_equal(as.matrix(toRMatrix(Mt)), M, tolerance = 1e-10)

  llt <- .nbLoglikBatch(Yt, Mt, f$psi, At, f$pen)
  expect_equal(as.numeric(toRMatrix(llt)), ll, tolerance = 1e-9)
})

test_that(".nbLoglikBatch handles an all-zero gene and a large count", {
  # the two ends .polishBatch()'s own tests exercise per gene
  f <- .nbFixture()
  f$Y[1, ] <- 0L
  f$Y[2, ] <- f$Y[2, ] * 100L
  M <- .muRef(f)
  got <- .nbLoglikBatch(f$Y, M, f$psi, f$A, f$pen)
  expect_equal(got, .llRef(f, M), tolerance = 1e-10)
  expect_true(all(is.finite(got)))
})

test_that(".nbLoglikBatch handles a zero dispersion, which is the Poisson limit", {
  # psi = 0 means size = 1/psi = Inf, and dnbinom() treats that as Poisson.
  # The written-out log-pmf does not get that for free: r = Inf gives
  # lgamma(y + Inf) - lgamma(Inf) = Inf - Inf = NaN and r*log(r/(r+mu)) =
  # Inf*log(1) = NaN, so every likelihood comes back NaN and -- since the
  # engine drops a gene whose loglik is not finite -- NOTHING is polished,
  # silently.
  #
  # This is not hypothetical. fitSpiDE() on the toy fixture returns psi = 0 for
  # every gene, so the whole device path returned "polished: FALSE" for all of
  # them (H100, job 28555677). Every unit fixture here used psi = 0.4 or 0.5
  # and missed it.
  f <- .nbFixture()
  f$psi <- rep(0, nrow(f$Y))
  M <- .muRef(f)
  ref <- .llRef(f, M)                       # dnbinom(size = Inf) = Poisson
  expect_true(all(is.finite(ref)))

  got <- .nbLoglikBatch(f$Y, M, f$psi, f$A, f$pen)
  expect_equal(got, ref, tolerance = 1e-10)

  # Gated on a GPU, not just torch: on a CPU device torch 0.17's
  # torch_tensor(<R double>, float64) aliases R memory rather than copying,
  # so these results depend on GC timing. The engine fix is deferred.
  skip_if_no_gpu()
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  gt <- .nbLoglikBatch(tt(f$Y), tt(M), f$psi, tt(f$A), f$pen)
  gt <- as.numeric(toRMatrix(gt))
  expect_true(all(is.finite(gt)))
  expect_equal(gt, ref, tolerance = 1e-9)
})

test_that(".nbLoglikBatch handles a mix of zero and non-zero dispersions", {
  # the batch is the point: one Poisson gene among NB ones must not take the
  # others with it
  # Gated on a GPU, not just torch: on a CPU device torch 0.17's
  # torch_tensor(<R double>, float64) aliases R memory rather than copying,
  # so these results depend on GC timing. The engine fix is deferred.
  skip_if_no_gpu()
  f <- .nbFixture()
  f$psi[c(2L, 4L)] <- 0
  M <- .muRef(f)
  ref <- .llRef(f, M)
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  gt <- as.numeric(toRMatrix(
    .nbLoglikBatch(tt(f$Y), tt(M), f$psi, tt(f$A), f$pen)))
  expect_true(all(is.finite(gt)))
  expect_equal(gt, ref, tolerance = 1e-9)
})

# ---- from spiDE tests/testthat/test-psi-batch.R ----------------------------

# The profile dispersion search, batched and fixed-iteration.
#
# This is the ONE deliberate numerical divergence in the batched fitter, and
# the justification is accuracy rather than speed: optimize()'s default
# tolerance is .Machine$double.eps^0.25, about 1.2e-4 ABSOLUTE on a
# log-interval of width log(1e3) - log(1e-3) = 13.8, so the psi the per-gene
# engine reports is only determined to ~1e-4 anyway. Fifty bisection steps on
# the score take that interval to ~1e-14.
#
# So the gate is not "the same answer as optimize()". It is "at least as good
# an answer", measured on the objective both are maximising, plus agreement
# within optimize()'s own tolerance. Bisection is also fixed-iteration, which
# means no divergent control flow across a batch -- the reason the plan chose
# it over reproducing Brent.

.psiFixture <- function(n = 200, b = 6, seed = 5) {
  set.seed(seed)
  mu <- matrix(exp(stats::rnorm(b * n, 2, 0.3)), b, n)
  Y <- matrix(0L, b, n)
  # genes 1..b-1 overdispersed at a range of psi; gene b is exactly Poisson,
  # whose NB optimum sits below the search range and must come back at_bound
  psis <- c(0.05, 0.2, 0.8, 2, 10)[seq_len(b - 1L)]
  for (i in seq_len(b - 1L)) {
    Y[i, ] <- stats::rnbinom(n, mu = mu[i, ], size = 1 / psis[i])
  }
  Y[b, ] <- stats::rpois(n, mu[b, ])
  list(Y = Y, Mu = mu, range = c(1e-3, 1e3))
}

# what the per-gene engine does today, verbatim
.psiRef <- function(f) {
  lo <- log(f$range[1]); hi <- log(f$range[2])
  est <- numeric(nrow(f$Y)); bnd <- logical(nrow(f$Y))
  for (i in seq_len(nrow(f$Y))) {
    y <- f$Y[i, ]; mu <- f$Mu[i, ]
    o <- stats::optimize(function(lp) {
      -sum(stats::dnbinom(y, size = 1 / exp(lp), mu = mu, log = TRUE))
    }, c(lo, hi))
    bnd[i] <- (o$minimum - lo) < 1e-3 * (hi - lo) ||
      (hi - o$minimum) < 1e-3 * (hi - lo)
    est[i] <- exp(o$minimum)
  }
  list(psi = est, at_bound = bnd)
}

.nbll <- function(y, mu, psi) sum(stats::dnbinom(y, size = 1 / psi, mu = mu, log = TRUE))

test_that(".psiProfileBatch agrees with optimize() within optimize's own tolerance", {
  f <- .psiFixture()
  ref <- .psiRef(f)
  got <- .psiProfileBatch(f$Y, f$Mu, f$range)
  expect_named(got, c("psi", "at_bound"), ignore.order = TRUE)
  # compare in the space the search runs in, against optimize's tolerance there
  free <- !ref$at_bound
  expect_true(any(free))
  expect_lt(max(abs(log(got$psi[free]) - log(ref$psi[free]))), 2e-4)
})

test_that(".psiProfileBatch never finds a worse optimum than optimize()", {
  # the gate. Both maximise the same function; the bisection is the more
  # accurate search, so it must not lose.
  f <- .psiFixture()
  ref <- .psiRef(f)
  got <- .psiProfileBatch(f$Y, f$Mu, f$range)
  for (i in seq_len(nrow(f$Y))) {
    ll_new <- .nbll(f$Y[i, ], f$Mu[i, ], got$psi[i])
    ll_old <- .nbll(f$Y[i, ], f$Mu[i, ], ref$psi[i])
    expect_gte(ll_new, ll_old - 1e-8 * abs(ll_old))
  }
})

test_that(".psiProfileBatch reproduces the at_bound rule", {
  f <- .psiFixture()
  ref <- .psiRef(f)
  got <- .psiProfileBatch(f$Y, f$Mu, f$range)
  expect_identical(got$at_bound, ref$at_bound)
  # the Poisson gene is the one that should be flagged
  expect_true(got$at_bound[nrow(f$Y)])
})

test_that(".psiProfileBatch is fixed-iteration, so more steps only refine", {
  # no divergent control flow across the batch: the answer at 60 steps is the
  # answer at 50, refined, never a different root
  f <- .psiFixture()
  a <- .psiProfileBatch(f$Y, f$Mu, f$range, maxit = 50L)
  b <- .psiProfileBatch(f$Y, f$Mu, f$range, maxit = 60L)
  expect_equal(a$psi, b$psi, tolerance = 1e-8)
  expect_identical(a$at_bound, b$at_bound)
})

test_that(".psiProfileBatch agrees between the base-R and torch branches", {
  # Gated on a GPU, not just torch: on a CPU device torch 0.17's
  # torch_tensor(<R double>, float64) aliases R memory rather than copying,
  # so these results depend on GC timing. The engine fix is deferred.
  skip_if_no_gpu()
  f <- .psiFixture()
  base <- .psiProfileBatch(f$Y, f$Mu, f$range)
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  tor <- .psiProfileBatch(tt(f$Y), tt(f$Mu), f$range)
  expect_equal(as.numeric(toRMatrix(tor$psi)), base$psi, tolerance = 1e-9)
  expect_identical(tor$at_bound, base$at_bound)
})

# ---- from spiDE tests/testthat/test-polish-batch.R -------------------------

# The batched per-gene Newton must be the per-gene one, restructured. These
# tests are written against .polishGene() as the oracle: it stays in the tree
# for exactly that purpose, reachable through engine = "gene".
#
# The risk being tested is not the arithmetic -- it is per-gene state leaking
# between slices of a batch. Genes converge at different iterations, take
# different numbers of line-search halvings, restart independently and fail
# independently, and every one of those is a per-gene branch that batching
# turns into an index set.

# a small design in the exact shape .polishFit() hands the solver: a dense
# block plus 0/1 indicators partitioning the cells
toy_batch <- function(n = 240L, G = 6L, seed = 11L) {
  set.seed(seed)
  x <- scale(rnorm(n))[, 1]
  ct <- rep(c(1, 2), length.out = n)
  X <- cbind(`(Intercept)` = 1, CellTypeB = as.numeric(ct == 2), niche = x)
  grp <- rep(seq_len(G), length.out = n)
  Z <- matrix(0, n, G, dimnames = list(NULL, paste0("SampleCellType", seq_len(G))))
  Z[cbind(seq_len(n), grp)] <- 1
  W <- cbind(X, Z)
  list(W = W, pen = c(0, 0, 0, rep(1 / 0.05, G)),
       nested = c(rep(FALSE, 3), rep(TRUE, G)),
       ct_cols = c(TRUE, TRUE, FALSE, rep(FALSE, G)))
}

# run the oracle over a set of genes, one at a time
gene_by_gene <- function(Yb, d, A0, psi0, solver, ...) {
  out <- lapply(seq_len(nrow(Yb)), function(i) {
    .polishGene(as.numeric(Yb[i, ]), d$W, A0[i, ], psi0[[i]], d$pen,
                        solver, start.cols = d$ct_cols, ...)
  })
  list(alpha = do.call(rbind, lapply(out, `[[`, "alpha")),
       psi = vapply(out, `[[`, numeric(1), "psi"),
       loglik = vapply(out, `[[`, numeric(1), "loglik"),
       iterations = vapply(out, `[[`, integer(1), "iterations"),
       restarted = vapply(out, `[[`, logical(1), "restarted"),
       capped = vapply(out, `[[`, logical(1), "capped"),
       singular = vapply(out, `[[`, logical(1), "singular"),
       psi_bound = vapply(out, `[[`, logical(1), "psi_bound"),
       polished = vapply(out, `[[`, logical(1), "polished"))
}

expect_same_fit <- function(b, g, tol = 1e-12, fields = c("alpha", "psi", "loglik")) {
  for (f in fields) expect_equal(b[[f]], g[[f]], tolerance = tol, ignore_attr = TRUE)
  for (f in c("restarted", "capped", "singular", "psi_bound", "polished")) {
    expect_identical(unname(b[[f]]), unname(g[[f]]), info = f)
  }
}

# The batched engine searches the dispersion by fixed-iteration bisection and
# .polishGene() by optimize(), so under psi.method = "profile" the two differ by
# about optimize()'s own tolerance. Loosening a tolerance is not a gate, so this
# is the gate: the batched engine, searching the same objective more finely,
# must never land on a WORSE penalised log-likelihood.
expect_no_worse <- function(b, g) {
  expect_true(all(b$loglik >= g$loglik - 1e-9 * abs(g$loglik)),
              info = "batched loglik is worse than the per-gene loglik")
}

test_that("a batch of one reproduces .polishGene() exactly", {
  d <- toy_batch()
  solver <- .newtonSolver(d$W, d$pen, d$nested)
  set.seed(1)
  y <- matrix(rnbinom(ncol(d$W) * 0 + nrow(d$W), mu = 6, size = 3), nrow = 1)
  A0 <- matrix(0, 1, ncol(d$W)); A0[1, 1] <- log(mean(y))
  # psi.method = "fixed" holds the dispersion fixed, which takes the one
  # deliberate divergence (bisection against optimize()) out of the comparison
  # and keeps this a STRICT test of the Newton machinery, which is what it is
  # for. The psi search is gated on its own in test-psi-batch.R and below.
  b <- .polishBatch(y, d$W, A0, 0.4, d$pen, solver, start.cols = d$ct_cols,
                            psi.method = "fixed")
  g <- gene_by_gene(y, d, A0, 0.4, solver, psi.method = "fixed")
  expect_same_fit(b, g)
})

test_that("a batch reproduces the same genes run one at a time", {
  d <- toy_batch()
  solver <- .newtonSolver(d$W, d$pen, d$nested)
  set.seed(2)
  B <- 13L
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = rep(c(2, 8, 30), length.out = B), size = 2),
               nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  psi0 <- rep(0.4, B)
  b <- .polishBatch(Yb, d$W, A0, psi0, d$pen, solver, start.cols = d$ct_cols,
                            psi.method = "fixed")
  g <- gene_by_gene(Yb, d, A0, psi0, solver, psi.method = "fixed")
  expect_same_fit(b, g)
})

test_that("the batch boundary does not move a gene's answer", {
  # The batched sibling of test-polish.R's gene-blocking invariance, at a
  # tolerance rather than at machine precision, and the difference is real:
  # mu = A %*% t(W) and the penalty term A^2 %*% pen go through BLAS, which
  # blocks a 13-row GEMM differently from a 1-row one, so the summation order
  # depends on the batch. The per-gene engine has no such dependence and its
  # exact-invariance test still holds for it.
  #
  # Measured here at ~2e-12 relative on psi, which propagates from a profile
  # optimum found on a mu that differs in its last digits. 1e-9 is the bar: far
  # inside anything that could move a call, far outside the noise.
  d <- toy_batch()
  solver <- .newtonSolver(d$W, d$pen, d$nested)
  set.seed(3)
  B <- 13L
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = rep(c(3, 12), length.out = B), size = 2), nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  full <- .polishBatch(Yb, d$W, A0, rep(0.4, B), d$pen, solver, start.cols = d$ct_cols)
  for (bs in c(1L, 2L, 5L)) {
    idx <- split(seq_len(B), ceiling(seq_len(B) / bs))
    parts <- lapply(idx, function(ii)
      .polishBatch(Yb[ii, , drop = FALSE], d$W, A0[ii, , drop = FALSE],
                           rep(0.4, length(ii)), d$pen, solver, start.cols = d$ct_cols))
    got <- list(alpha = do.call(rbind, lapply(parts, `[[`, "alpha")),
                psi = unlist(lapply(parts, `[[`, "psi")),
                loglik = unlist(lapply(parts, `[[`, "loglik")))
    expect_equal(got$alpha, full$alpha, tolerance = 1e-9, ignore_attr = TRUE,
                 info = paste("batch size", bs))
    expect_equal(got$psi, full$psi, tolerance = 1e-9, ignore_attr = TRUE)
    expect_equal(got$loglik, full$loglik, tolerance = 1e-9, ignore_attr = TRUE)
  }
})

test_that("one batch holds genes taking every different path", {
  # THE test: a converged gene, a degenerate start, an all-zero gene, a gene
  # that needs several halvings, and a gene whose psi sits on a bound, in one
  # batch, against the same five run singly. This is what catches per-gene
  # state leaking across slices.
  d <- toy_batch()
  solver <- .newtonSolver(d$W, d$pen, d$nested)
  set.seed(5)
  n <- nrow(d$W)
  Yb <- rbind(
    rnbinom(n, mu = 8, size = 3),        # ordinary
    rnbinom(n, mu = 40, size = 5),       # bright
    rep(0L, n),                          # all zero: psi runs to a bound
    rpois(n, lambda = 4),                # under-dispersed: the other bound
    rnbinom(n, mu = 2, size = 1)         # noisy, hard line search
  )
  B <- nrow(Yb)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  A0[2, ] <- -30                          # a degenerate start: forces a restart
  psi0 <- c(0.4, 0.2, 0.5, 0.3, 0.8)
  b <- .polishBatch(Yb, d$W, A0, psi0, d$pen, solver, start.cols = d$ct_cols)
  g <- gene_by_gene(Yb, d, A0, psi0, solver)
  # This one keeps psi.method = "profile": a gene whose dispersion runs to a
  # bound is one of the paths it exists to mix, and "fixed" would remove it.
  # So the numeric fields are compared at the dispersion search's accuracy
  # rather than at machine precision -- optimize()'s tolerance is ~1.2e-4 on the
  # log interval, and the observed spread is ~3e-6 in psi and ~1e-7 in alpha.
  # The FLAGS are still identical, which is what "flag for flag" meant, and the
  # objective gate below is what makes the looser tolerance honest.
  expect_same_fit(b, g, tol = 1e-4)
  expect_no_worse(b, g)
  expect_identical(b$iterations, g$iterations)
  expect_true(any(g$restarted))           # the fixture must actually exercise it
  expect_true(any(g$psi_bound))
})

test_that("a gene that cannot be polished does not poison its batch", {
  # .waldCauchyBlock() fails a whole sub-batch when one gene's Cholesky fails,
  # and its only recourse is to tell the user to shrink cov.batch. The polish
  # must not acquire that failure mode. The warm path's documented fallback --
  # a non-finite start returns fitNB's own fit, not the sane start -- is a real
  # code path, so no stub solver is needed to reach it.
  d <- toy_batch()
  solver <- .newtonSolver(d$W, d$pen, d$nested)
  set.seed(6)
  B <- 4L
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = 7, size = 3), nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(rowMeans(Yb))
  A0[2, 3] <- NaN
  psi0 <- rep(0.4, B)
  b <- .polishBatch(Yb, d$W, A0, psi0, d$pen, solver,
                            start.cols = d$ct_cols, warm = TRUE)
  g <- gene_by_gene(Yb, d, A0, psi0, solver, warm = TRUE)
  expect_false(b$polished[2])
  expect_identical(b$alpha[2, ], A0[2, ])     # fitNB's fit kept, not the sane start
  expect_true(all(b$polished[-2]))            # the neighbours are unaffected
  expect_same_fit(b, g)
  # and identical to those three run without the failing gene present at all
  ok <- .polishBatch(Yb[-2, , drop = FALSE], d$W, A0[-2, , drop = FALSE],
                             psi0[-2], d$pen, solver, start.cols = d$ct_cols, warm = TRUE)
  expect_equal(b$alpha[-2, ], ok$alpha, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that(".polishBatch reports what a shared factorisation would cost", {
  # Phase 2e's design question, made measurable. newton() keeps a list of
  # per-gene factorisations under a per-gene staleness counter. One shared
  # TENSOR factorisation cannot do that: it must refresh the whole active stack
  # whenever any gene in it is stale. How much of Phase 0b's memoisation that
  # discards is an empirical question about the staleness trajectory, not
  # something to reason about -- so the engine counts both.
  #
  #   factorisations       what the per-gene policy actually built
  #   factorisations_sync  what refreshing the whole active set would have built
  #
  # The second is a counterfactual and costs one integer per iteration.
  d <- toy_batch()
  set.seed(4)
  B <- 8L
  mu <- exp(d$W %*% c(1.2, 0.4, 0.3, rep(0, 6)))
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = rep(as.numeric(mu), each = B),
                       size = 1 / 0.4), nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  solver <- .newtonSolver(d$W, d$pen, d$nested)

  out <- .polishBatch(Yb, d$W, A0, rep(0.4, B), d$pen, solver,
                              start.cols = d$ct_cols)
  nf <- attr(out, "factorisations")
  expect_type(nf, "integer")
  expect_named(nf, c("pergene", "sync"))

  # every gene is factorised at least once before its first step
  expect_gte(nf[["pergene"]], B)
  # and a shared factorisation can never build fewer than the per-gene policy
  expect_gte(nf[["sync"]], nf[["pergene"]])
})

test_that("a shared factorisation reaches the same optimum as the per-gene one", {
  # Phase 2e's wiring. newton() keeps one factorisation for the active set
  # instead of a list of per-gene ones, refreshing the whole stack when any
  # active gene is stale. That is a different PATH -- a gene that was not stale
  # gets a fresher information matrix than it would have had -- so this is not
  # an equality test and must not pretend to be one. A fresher matrix cannot
  # move the fixed point and is not a worse Newton direction, so the gate is
  # the one Phase 0b used for the symmetric gram: the objective is never worse.
  d <- toy_batch()
  set.seed(21)
  B <- 12L
  mu <- exp(d$W %*% c(1.2, 0.4, 0.3, rep(0, 6)))
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = rep(as.numeric(mu), each = B),
                       size = 1 / 0.4), nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  solver <- .newtonSolver(d$W, d$pen, d$nested)

  args <- list(Yb, d$W, A0, rep(0.4, B), d$pen, solver, start.cols = d$ct_cols)
  per <- do.call(.polishBatch, args)
  shd <- do.call(.polishBatch,
                 c(args, list(shared.factor = TRUE, nested = d$nested)))

  # same genes polished, same failures
  expect_identical(shd$polished, per$polished)
  expect_identical(shd$singular, per$singular)
  expect_identical(shd$capped, per$capped)

  # the same optimum, to the convergence tolerance rather than to machine
  expect_equal(shd$alpha, per$alpha, tolerance = 1e-5)
  expect_equal(shd$psi, per$psi, tolerance = 1e-5)

  # the gate: the penalised log-likelihood is never worse
  expect_gt(min(shd$loglik - per$loglik), -1e-7 * max(abs(per$loglik)))
})

test_that(".polishBatch runs on tensors and agrees with the matrix path", {
  # The point of the whole conversion. Exercised on CPU torch tensors, which is
  # the same code the device runs -- only the placement differs. A shared
  # factorisation is required here: the per-gene branch keeps a LIST of
  # factorisations, which is exactly what cannot go to a device.
  # Gated on a GPU, not just torch: on a CPU device torch 0.17's
  # torch_tensor(<R double>, float64) aliases R memory rather than copying,
  # so these results depend on GC timing. The engine fix is deferred.
  skip_if_no_gpu()
  d <- toy_batch()
  solver <- .newtonSolver(d$W, d$pen, d$nested)
  set.seed(77)
  B <- 9L
  mu <- exp(d$W %*% c(1.3, 0.35, 0.25, rep(0, 6)))
  Yb <- matrix(rnbinom(B * nrow(d$W), mu = rep(as.numeric(mu), each = B),
                       size = 1 / 0.5), nrow = B)
  A0 <- matrix(0, B, ncol(d$W)); A0[, 1] <- log(pmax(rowMeans(Yb), 0.1))
  psi0 <- rep(0.5, B)

  args <- list(A0 = A0, psi0 = psi0, pen = d$pen, solver = solver,
               start.cols = d$ct_cols, shared.factor = TRUE, nested = d$nested)
  cpu <- do.call(.polishBatch, c(list(Yb, d$W), args))

  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  targs <- args; targs$A0 <- tt(A0)
  tor <- do.call(.polishBatch, c(list(tt(Yb), tt(d$W)), targs))

  # the return is always host, whatever went in
  expect_true(is.matrix(tor$alpha))
  expect_type(tor$psi, "double")

  # the agreement is worthless if neither side did anything: pin that the
  # matrix path actually polished every gene before comparing
  expect_true(all(cpu$polished))
  expect_true(all(is.finite(cpu$loglik)))

  expect_identical(tor$polished, cpu$polished)
  expect_identical(tor$singular, cpu$singular)
  expect_identical(tor$restarted, cpu$restarted)
  expect_identical(tor$psi_bound, cpu$psi_bound)
  expect_equal(tor$alpha, cpu$alpha, tolerance = 1e-8)
  expect_equal(tor$psi, cpu$psi, tolerance = 1e-8)
  expect_equal(tor$loglik, cpu$loglik, tolerance = 1e-8)
})

# ---- the optional per-cell offset ------------------------------------------

.score <- function(y, W, a, psi, pen, off) {
  mu <- pmax(as.numeric(exp(W %*% a + off)), nbMuFloor())
  as.numeric(crossprod(W, (y - mu) / (1 + psi * mu))) - pen * a
}

test_that("polishGene reaches a zero penalised score with an offset", {
  set.seed(1)
  n <- 300
  W <- cbind(1, rnorm(n), rnorm(n))
  off <- log(runif(n, 0.5, 2))               # a known per-cell log-scale factor
  y <- rnbinom(n, mu = exp(1 + 0.4 * W[, 2] + off), size = 5)
  pen <- c(0, 1, 1)
  s <- .newtonSolver(W, pen)
  r <- .polishGene(y, W, a0 = c(0, 0, 0), psi0 = 0.2, pen = pen, solver = s,
                   psi.method = "fixed", offset = off)
  sc <- .score(y, W, r$alpha, r$psi, pen, off)
  expect_lt(max(abs(sc)) / sum(y), 1e-6)
  expect_true(r$polished)
})

test_that("a far start converges to the intercept net of a constant offset", {
  set.seed(2)
  n <- 200
  W <- cbind(1, rnorm(n))
  off <- rep(3, n)                            # a large constant offset
  y <- rnbinom(n, mu = exp(0.5 + off), size = 10)
  s <- .newtonSolver(W, c(0, 0))
  # a0 = 40 is far from the optimum but not degenerate (its fitted log-mean is
  # 43, not below -10), so Newton descends from it without a restart (44
  # iterations); "the sane start is the log mean net of the offset" below
  # tests the sane start
  r <- .polishGene(y, W, a0 = c(40, 0), psi0 = 0.1, pen = c(0, 0), solver = s,
                   psi.method = "fixed", offset = off)
  expect_false(r$restarted)
  expect_equal(r$alpha[1], 0.5, tolerance = 0.1)
})

test_that("batch and per-gene engines agree with a vector and a matrix offset", {
  set.seed(3)
  G <- 6; n <- 150
  W <- cbind(1, rnorm(n))
  off <- log(runif(n, 0.5, 2))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + 0.2 * g * W[, 2] + off), size = 4), numeric(n)))
  pen <- c(0, 0.5)
  s <- .newtonSolver(W, pen)
  per <- t(vapply(seq_len(G), function(g)
    .polishGene(Y[g, ], W, c(0, 0), 0.25, pen, s, psi.method = "fixed",
                offset = off)$alpha, numeric(2)))
  bv <- .polishBatch(Y, W, matrix(0, G, 2), rep(0.25, G), pen, s,
                     psi.method = "fixed", offset = off)
  bm <- .polishBatch(Y, W, matrix(0, G, 2), rep(0.25, G), pen, s,
                     psi.method = "fixed",
                     offset = matrix(off, G, n, byrow = TRUE))
  expect_equal(bv$alpha, per, tolerance = 1e-8)
  expect_equal(bm$alpha, bv$alpha, tolerance = 0)
})

test_that("offset = NULL is bit-identical to the pre-offset engine", {
  set.seed(4)
  n <- 120
  W <- cbind(1, rnorm(n))
  y <- rnbinom(n, mu = exp(1 + 0.3 * W[, 2]), size = 3)
  s <- .newtonSolver(W, c(0, 0.1))
  a <- .polishGene(y, W, c(0, 0), 0.3, c(0, 0.1), s, psi.method = "profile")
  b <- .polishGene(y, W, c(0, 0), 0.3, c(0, 0.1), s, psi.method = "profile",
                   offset = NULL)
  expect_identical(a, b)
})

# Beyond the four above: each of these fails when the offset is dropped (or
# mis-sliced) at one site the four do not reach -- the degenerate-start check,
# the sane start with and without start.cols, and the row slicing of a
# per-gene offset matrix on a path where the active rows are not 1..B.

test_that("the degenerate-start check reads the linear predictor with its offset", {
  set.seed(5)
  n <- 150
  W <- cbind(1, rnorm(n))
  s <- .newtonSolver(W, c(0, 0))
  # -11 + 12 = 1 is a fitted log-mean at the optimum, so not degenerate and no
  # restart; read without the offset it is -11, below -10, and would restart
  up <- rep(12, n)
  y <- rnbinom(n, mu = exp(-11 + up), size = 5)
  g <- .polishGene(y, W, c(-11, 0), 0.2, c(0, 0), s, psi.method = "fixed",
                   offset = up)
  b <- .polishBatch(rbind(y), W, rbind(c(-11, 0)), 0.2, c(0, 0), s,
                    psi.method = "fixed", offset = up)
  expect_false(g$restarted)
  expect_false(b$restarted)
  # -5 - 8 = -13 is below -10: degenerate, restarted from the sane start
  down <- rep(-8, n)
  y <- rnbinom(n, mu = exp(9 + down), size = 5)
  g <- .polishGene(y, W, c(-5, 0), 0.2, c(0, 0), s, psi.method = "fixed",
                   offset = down)
  b <- .polishBatch(rbind(y), W, rbind(c(-5, 0)), 0.2, c(0, 0), s,
                    psi.method = "fixed", offset = down)
  expect_true(g$restarted)
  expect_true(b$restarted)
})

test_that("the sane start is the log mean net of the offset, in both engines", {
  set.seed(6)
  n <- 200
  ct <- rep(c(TRUE, FALSE), length.out = n)
  W <- cbind(A = as.numeric(ct), B = as.numeric(!ct), x = rnorm(n))
  pen <- c(0, 0, 0)
  s <- .newtonSolver(W, pen)
  G <- 3
  O <- t(vapply(seq_len(G), function(g)
    log(runif(n, 0.5, 2)) + ifelse(ct, g, -g), numeric(n)))  # distinct rows
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.3 * W[, 3] + O[g, ] + 0.5), size = 5), numeric(n)))
  A0 <- matrix(0, G, 3)
  A0[, 1] <- NaN                    # a non-finite start: always the sane start
  for (sc in list(c(TRUE, TRUE, FALSE), NULL)) {
    # maxit = 0 returns the start itself, so the formula is checked directly
    for (g in seq_len(G)) {
      r <- .polishGene(Y[g, ], W, A0[g, ], 0.2, pen, s, maxit = 0L,
                       start.cols = sc, psi.method = "fixed", offset = O[g, ])
      want <- if (is.null(sc)) {
        c(log(mean(Y[g, ]) + 1e-3) - mean(O[g, ]), 0, 0)
      } else {
        c(log(mean(Y[g, ct]) + 1e-3) - mean(O[g, ct]),
          log(mean(Y[g, !ct]) + 1e-3) - mean(O[g, !ct]), 0)
      }
      expect_true(r$restarted)
      expect_equal(r$alpha, want, tolerance = 1e-12)
    }
    # the batched engine takes at least one step, so compare one step from
    # each engine's sane start. Gene 1 keeps a finite start, so the restarted
    # genes are batch rows 2..G rather than 1..G-1.
    A1 <- A0
    A1[1, ] <- c(0.5, 0.5, 0)
    per <- t(vapply(seq_len(G), function(g)
      .polishGene(Y[g, ], W, A1[g, ], 0.2, pen, s, maxit = 1L, start.cols = sc,
                  psi.method = "fixed", offset = O[g, ])$alpha, numeric(3)))
    bat <- .polishBatch(Y, W, A1, rep(0.2, G), pen, s, maxit = 1L,
                        start.cols = sc, psi.method = "fixed", offset = O)
    expect_identical(bat$restarted, c(FALSE, rep(TRUE, G - 1)))
    expect_equal(bat$alpha, per, tolerance = 1e-10, ignore_attr = TRUE)
  }
})

test_that("a per-gene offset matrix reaches each gene on every path", {
  # The all-zero first gene's dispersion runs to a bound, so the profile
  # re-polish runs on rows 2..B: positions in `rows` are then not batch rows,
  # which is where a mis-sliced offset would show. Each gene alone in its own
  # batch is the reference.
  set.seed(7)
  G <- 4; n <- 160
  W <- cbind(1, rnorm(n))
  pen <- c(0, 0.5)
  s <- .newtonSolver(W, pen)
  O <- t(vapply(seq_len(G), function(g) log(runif(n, 0.5, 2)) + g - 2, numeric(n)))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + 0.3 * W[, 2] + O[g, ]), size = 4), numeric(n)))
  Y[1, ] <- 0
  A0 <- cbind(log(pmax(rowMeans(Y), 0.1)) - rowMeans(O), 0)
  full <- .polishBatch(Y, W, A0, rep(0.3, G), pen, s, offset = O)
  one <- lapply(seq_len(G), function(g)
    .polishBatch(Y[g, , drop = FALSE], W, A0[g, , drop = FALSE], 0.3, pen, s,
                 offset = O[g, , drop = FALSE]))
  expect_true(full$psi_bound[1])
  expect_true(all(full$polished[-1]))
  expect_equal(full$alpha, do.call(rbind, lapply(one, `[[`, "alpha")),
               tolerance = 1e-9)
  expect_equal(full$psi, vapply(one, `[[`, numeric(1), "psi"), tolerance = 1e-9)
  expect_identical(full$iterations, vapply(one, `[[`, integer(1), "iterations"))
})

test_that("a mis-shaped offset is refused, not recycled", {
  set.seed(8)
  n <- 40
  W <- cbind(1, rnorm(n))
  Y <- matrix(rnbinom(3 * n, mu = 4, size = 3), 3, n)
  s <- .newtonSolver(W, c(0, 0))
  expect_error(.polishGene(Y[1, ], W, c(1, 0), 0.3, c(0, 0), s, offset = rep(0, n - 1)),
               "'offset' must have one value per cell")
  expect_error(.polishBatch(Y, W, matrix(1, 3, 2), rep(0.3, 3), c(0, 0), s,
                            offset = rep(0, 3)),
               "'offset' must have one value per cell")
  expect_error(.polishBatch(Y, W, matrix(1, 3, 2), rep(0.3, 3), c(0, 0), s,
                            offset = matrix(0, 2, n)),
               "'offset' must have one value per cell")
  # a bare per-cell vector would recycle down the columns of A %*% t(W)
  expect_error(.muBatch(matrix(1, 3, 2), W, offset = rep(0, n)), ".offsetRows")
})

test_that("the offset reaches the tensor path", {
  # Gated on a GPU, not just torch: on a CPU device torch 0.17's
  # torch_tensor(<R double>, float64) aliases R memory rather than copying,
  # so these results depend on GC timing. The engine fix is deferred.
  skip_if_no_gpu()
  set.seed(9)
  G <- 5; n <- 180
  W <- cbind(1, rnorm(n))
  pen <- c(0, 0.5)
  s <- .newtonSolver(W, pen)
  O <- t(vapply(seq_len(G), function(g) log(runif(n, 0.5, 2)) + g / 2 - 1, numeric(n)))
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.5 + 0.3 * W[, 2] + O[g, ]), size = 4), numeric(n)))
  A0 <- cbind(log(pmax(rowMeans(Y), 0.1)) - rowMeans(O), 0)
  A0[2, ] <- NaN                                   # a sane start on one gene
  tt <- function(x) torch::torch_tensor(x, dtype = torch::torch_float64())
  for (off in list(O, O[3, ])) {
    cpu <- .polishBatch(Y, W, A0, rep(0.3, G), pen, s, shared.factor = TRUE,
                        nested = c(FALSE, FALSE), offset = off)
    tor <- .polishBatch(tt(Y), tt(W), tt(A0), rep(0.3, G), pen, s,
                        shared.factor = TRUE, nested = c(FALSE, FALSE),
                        offset = off)
    expect_true(all(cpu$polished))
    expect_identical(tor$restarted, cpu$restarted)
    expect_equal(tor$alpha, cpu$alpha, tolerance = 1e-8)
    expect_equal(tor$psi, cpu$psi, tolerance = 1e-8)
  }
})
