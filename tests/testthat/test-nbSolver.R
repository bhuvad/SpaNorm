# Tests for the penalised NB Newton solvers and batched grams, moved from
# spiDE: test-solver-batch.R, test-solver-slope-absorb.R, test-absorb-batch.R.
# Only the tests that build `W` directly moved here; any test that went
# through fitSpiDE()/spiDE()/toySpiDE/.toyClustered stayed in spiDE (none of
# the three source files had one -- see task-2-report.md).
#
# skip_if_no_torch()/gpu_tol() live in helper-gpu.R (auto-sourced by testthat
# before any test file runs). The torch-agreement tests below use
# skip_if_no_torch() in place of spiDE's skip_if_not_installed("torch"): a
# torch install without a usable libtorch/lantern backend must skip too.

## ---- from test-solver-batch.R --------------------------------------------
# A batched factor/solve for the nested design, which is what Phase 2e needs to
# put .polishBatch()'s Newton on a device.
#
# newton() currently keeps a LIST of per-gene factorisations and calls
# solver$solve(fac[[k]], S[k, ]) one gene at a time. A list of R objects cannot
# go to a device and a per-gene solve is a kernel launch per gene per
# iteration, so the state has to become one batched object.
#
# .newtonSolver() is the oracle throughout: same Schur absorption, same
# right-hand side, same back-substitution, one gene at a time.

.solverFixture <- function(n = 80, px = 4, G = 5, b = 4, seed = 3,
                           nested = TRUE) {
  set.seed(seed)
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  if (!nested) {
    return(list(W = X, px = px, G = 0L,
                nested = rep(FALSE, px), pen = rep(0, px),
                wt = matrix(stats::runif(b * n, 0.2, 2), b, n)))
  }
  grp <- c(seq_len(G), sample(seq_len(G), n - G, replace = TRUE))
  Z <- matrix(0, n, G)
  Z[cbind(seq_len(n), grp)] <- 1
  # pen_x = 0 is the production shape: the ridge is on the random-effect
  # columns only, which is also what makes a singular dense block reachable
  list(W = cbind(X, Z), px = px, G = G,
       nested = c(rep(FALSE, px), rep(TRUE, G)),
       pen = c(rep(0, px), rep(0.3, G)),
       wt = matrix(stats::runif(b * n, 0.2, 2), b, n))
}

# a per-gene score to solve against
.scores <- function(f) {
  set.seed(99)
  matrix(stats::rnorm(nrow(f$wt) * ncol(f$W)), nrow(f$wt), ncol(f$W))
}

test_that(".absorbBatch can return the parts, not only the Schur complement", {
  # the solve needs B and cvec as well as S, and they must be the same B and
  # cvec .newtonSolver() built
  f <- .solverFixture()
  got <- .absorbBatch(f$W, f$pen, f$nested, f$wt, parts = TRUE)
  expect_named(got, c("S", "B", "cvec", "xi", "zi"), ignore.order = TRUE)
  expect_equal(dim(got$S), c(nrow(f$wt), f$px, f$px))
  expect_equal(dim(got$B), c(nrow(f$wt), f$G, f$px))
  expect_equal(dim(got$cvec), c(nrow(f$wt), f$G))

  sol <- .newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))) {
    ref <- sol$factor(f$wt[g, ])
    expect_equal(got$S[g, , ], ref$S, tolerance = 1e-10)
    # .newtonSolver()'s B is px x G; the batched stack carries its transpose.
    # unname(): rowsum() puts the group levels on B's columns and cvec's names,
    # which a (batch, G, px) stack cannot carry and nothing downstream reads --
    # both index positionally.
    expect_equal(t(got$B[g, , ]), unname(ref$B), tolerance = 1e-10)
    expect_equal(got$cvec[g, ], unname(ref$cvec), tolerance = 1e-10)
  }
  # and the default return is unchanged
  expect_equal(.absorbBatch(f$W, f$pen, f$nested, f$wt), got$S)
})

test_that(".newtonSolverBatch's step matches .newtonSolver's, gene by gene", {
  f <- .solverFixture()
  Sc <- .scores(f)
  sb <- .newtonSolverBatch(f$W, f$pen, f$nested)
  st <- sb$factor(f$wt)
  D <- sb$solve(st, Sc)
  expect_equal(dim(D), dim(Sc))
  expect_true(all(st$ok))

  sol <- .newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))) {
    expect_equal(D[g, ], sol$solve(f$wt[g, ], Sc[g, ]), tolerance = 1e-9)
  }
})

test_that(".newtonSolverBatch works with no nested block at all", {
  f <- .solverFixture(nested = FALSE)
  Sc <- .scores(f)
  sb <- .newtonSolverBatch(f$W, f$pen, f$nested)
  D <- sb$solve(sb$factor(f$wt), Sc)
  sol <- .newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))) {
    expect_equal(D[g, ], sol$solve(f$wt[g, ], Sc[g, ]), tolerance = 1e-9)
  }
})

test_that(".newtonSolverBatch's xcov matches .newtonSolver's", {
  f <- .solverFixture()
  sb <- .newtonSolverBatch(f$W, f$pen, f$nested)
  V <- sb$xcov(sb$factor(f$wt))
  sol <- .newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))) {
    expect_equal(V[g, , ], sol$xcov(f$wt[g, ]), tolerance = 1e-9)
  }
})

test_that("a singular gene does not poison its batch", {
  # THE trap this whole batched path has to avoid, and the one
  # .waldCauchyBlock() fell into (inference.R: one singular gene kills the
  # sub-batch's Cholesky and the only recourse is telling the user to shrink
  # cov.batch). With pen_x = 0, a gene carrying no weight has a zero dense
  # block: singular, and its neighbours must be unaffected.
  f <- .solverFixture()
  f$wt[2, ] <- 0
  Sc <- .scores(f)
  sb <- .newtonSolverBatch(f$W, f$pen, f$nested)
  st <- sb$factor(f$wt)
  expect_false(st$ok[2])
  expect_true(all(st$ok[-2]))

  D <- sb$solve(st, Sc)
  expect_true(all(is.na(D[2, ])))
  sol <- .newtonSolver(f$W, f$pen, f$nested)
  for (g in seq_len(nrow(f$wt))[-2]) {
    expect_equal(D[g, ], sol$solve(f$wt[g, ], Sc[g, ]), tolerance = 1e-9)
  }
})

test_that(".newtonSolverBatch agrees between the base-R and torch branches", {
  skip_if_no_torch()
  f <- .solverFixture()
  Sc <- .scores(f)
  base_D <- .newtonSolverBatch(f$W, f$pen, f$nested) |>
    (\(sb) sb$solve(sb$factor(f$wt), Sc))()

  Wt <- torch::torch_tensor(f$W, dtype = torch::torch_float64())
  wtt <- torch::torch_tensor(f$wt, dtype = torch::torch_float64())
  Sct <- torch::torch_tensor(Sc, dtype = torch::torch_float64())
  sbt <- .newtonSolverBatch(Wt, f$pen, f$nested)
  tor_D <- sbt$solve(sbt$factor(wtt), Sct)
  expect_true(is_torch_tensor(tor_D))
  expect_equal(as.array(tor_D), base_D, tolerance = gpu_tol())
})

## ---- from test-solver-slope-absorb.R --------------------------------------
# Absorbing a random block that is BLOCK-diagonal by sample, not merely diagonal.
#
# .newtonSolver()'s Schur absorption was written for the nested (sample x cell
# type) intercept: 0/1 indicators partitioning the cells, so C = Z' diag(w) Z is
# diagonal and C^-1 is a reciprocal. Random SLOPES are indicator x covariate, so
# a sample's slope columns are not orthogonal to each other or to that sample's
# intercept -- C is block-diagonal by sample with a dense block per sample, and
# the scalar path cannot absorb it.
#
# It is still absorbable, and by the same identity: every random column belongs
# to exactly one sample, so cells of sample s load on no other sample's columns.
# The only change is that C^-1 is a per-sample dense solve rather than a
# reciprocal.
#
# The oracle throughout is the DENSE solve on the full design. Absorption is an
# exact re-arrangement, so "exact" is the standard, not "close".

# A design whose random block is block-diagonal by sample and NOT a 0/1
# partition: per-sample intercept plus per-sample slopes on K continuous bases.
.slopeFixture <- function(n = 120, px = 3, S = 4, K = 2, seed = 11) {
  set.seed(seed)
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  smp <- factor(c(seq_len(S), sample(seq_len(S), n - S, replace = TRUE)))
  Zint <- stats::model.matrix(~ 0 + smp)
  bases <- matrix(stats::rnorm(n * K), n, K)          # the CellType:niche bases
  Zsl <- do.call(cbind, lapply(seq_len(K), function(k) Zint * bases[, k]))
  Z <- cbind(Zint, Zsl)
  # block id per column: every random column belongs to exactly one sample
  blk <- c(as.integer(smp[!duplicated(smp)][order(unique(as.integer(smp)))]),
           rep(NA_integer_, 0))
  blk <- c(seq_len(S), rep(seq_len(S), times = K))
  list(W = cbind(X, Z), px = px, S = S, K = K,
       group = c(rep(NA_integer_, px), blk),
       pen = c(rep(0, px), rep(0.4, ncol(Z))),
       w = stats::runif(n, 0.2, 2),
       score = stats::rnorm(px + ncol(Z)))
}

# the dense oracle: (X'WX + diag(pen))^-1 s, no absorption anywhere
.denseSolve <- function(W, pen, w, s) {
  info <- crossprod(W * sqrt(w))
  diag(info) <- diag(info) + pen
  as.numeric(solve(info, s))
}

test_that(".newtonSolver absorbs a per-sample block-diagonal random block exactly", {
  f <- .slopeFixture()
  # sanity: this really is NOT the case the scalar path handles -- a sample's
  # slope column is not orthogonal to that sample's intercept
  zi <- which(!is.na(f$group))
  C <- crossprod(f$W[, zi, drop = FALSE] * sqrt(f$w))
  expect_gt(max(abs(C[upper.tri(C)])), 1e-6)

  sol <- .newtonSolver(f$W, f$pen, f$group)
  got <- sol$solve(f$w, f$score)
  expect_equal(got, .denseSolve(f$W, f$pen, f$w, f$score), tolerance = 1e-8)
})

test_that(".newtonSolver's xcov on a block-diagonal random block is the fixed-effect covariance", {
  f <- .slopeFixture()
  sol <- .newtonSolver(f$W, f$pen, f$group)
  got <- sol$xcov(f$w)
  info <- crossprod(f$W * sqrt(f$w)); diag(info) <- diag(info) + f$pen
  # the Schur complement's inverse IS the fixed-effect block of the full inverse
  expect_equal(got, solve(info)[seq_len(f$px), seq_len(f$px)], tolerance = 1e-8)
})

test_that("a logical `nested` still absorbs the indicator block exactly", {
  # the existing call sites pass a logical vector; that path must not move
  set.seed(5)
  n <- 90; px <- 3; G <- 5
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  g <- c(seq_len(G), sample(seq_len(G), n - G, replace = TRUE))
  Z <- matrix(0, n, G); Z[cbind(seq_len(n), g)] <- 1
  W <- cbind(X, Z); pen <- c(rep(0, px), rep(0.3, G))
  nested <- c(rep(FALSE, px), rep(TRUE, G))
  w <- stats::runif(n, 0.2, 2); s <- stats::rnorm(px + G)

  sol <- .newtonSolver(W, pen, nested)
  expect_equal(sol$solve(w, s), .denseSolve(W, pen, w, s), tolerance = 1e-8)
})

test_that("a block grouping that mixes indicators and slopes absorbs exactly", {
  # the production shape: SampleInt + SampleSlope + SampleCellTypeInt, all of
  # which belong to one sample, so all of them go in that sample's block
  set.seed(7)
  n <- 150; px <- 3; S <- 3; K <- 2; nct <- 2
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  smp <- factor(sample(seq_len(S), n, replace = TRUE))
  ct <- factor(sample(seq_len(nct), n, replace = TRUE))
  Zint <- stats::model.matrix(~ 0 + smp)
  bases <- matrix(stats::rnorm(n * K), n, K)
  Zsl <- do.call(cbind, lapply(seq_len(K), function(k) Zint * bases[, k]))
  grp <- interaction(smp, ct, drop = TRUE)
  Zct <- stats::model.matrix(~ 0 + grp)
  ct_sample <- as.integer(sub("\\..*$", "", levels(grp)))

  W <- cbind(X, Zint, Zsl, Zct)
  pen <- c(rep(0, px), rep(0.5, ncol(Zint) + ncol(Zsl) + ncol(Zct)))
  group <- c(rep(NA_integer_, px), seq_len(S), rep(seq_len(S), times = K), ct_sample)
  w <- stats::runif(n, 0.2, 2); s <- stats::rnorm(ncol(W))

  sol <- .newtonSolver(W, pen, group)
  expect_equal(sol$solve(w, s), .denseSolve(W, pen, w, s), tolerance = 1e-8)
})

test_that("the batched solver refuses a multi-column block rather than guessing", {
  # .absorbBatch()/.newtonSolverBatch() carry their own copy of the absorption
  # and only implement the diagonal (1x1 block) case -- C^-1 is a reciprocal
  # there, and a per-block Cholesky has no batched equivalent written yet. A
  # grouping they cannot honour must stop with a message that says so, because
  # the failure mode otherwise is a wrong Newton step, not an error.
  f <- .slopeFixture()
  expect_error(.newtonSolverBatch(f$W, f$pen, f$group),
               "1x1|one column|block", ignore.case = TRUE)
})

test_that("the batched solver still accepts the logical indicator case", {
  set.seed(5)
  n <- 60; px <- 3; G <- 4
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  g <- c(seq_len(G), sample(seq_len(G), n - G, replace = TRUE))
  Z <- matrix(0, n, G); Z[cbind(seq_len(n), g)] <- 1
  W <- cbind(X, Z); pen <- c(rep(0, px), rep(0.3, G))
  nested <- c(rep(FALSE, px), rep(TRUE, G))
  expect_no_error(.newtonSolverBatch(W, pen, nested))
})

## ---- from test-absorb-batch.R ----------------------------------------------
# The nested (sample x cell type) Schur absorption, batched and
# backend-agnostic. .newtonSolver()'s per-gene `parts()` is the oracle: it is
# the implementation every converged fit and the CPU inference path already
# use, and test-polish.R (spiDE) pins it against the dense inverse.
#
# The point of the batched form is the GPU inference path, which today skips
# the absorption entirely and builds a dense p x p gram -- 1,107 columns where
# 398 would do on the cohort's design, 7.7x the flops. Everything here is
# testable on CPU torch tensors; only device placement needs an accelerator.

# a design with a genuine nested indicator block: px dense columns, then G
# 0/1 columns that partition the cells
.absorbFixture <- function(n = 60, px = 4, G = 5, b = 3, seed = 7,
                           empty_group = FALSE) {
  set.seed(seed)
  X <- cbind(1, matrix(stats::rnorm(n * (px - 1)), n, px - 1))
  grp <- if (empty_group) {
    # group 2 gets no cells: rowsum() drops it, and the caller indexes cvec
    # positionally, so a dropped group silently misaligns every group after it
    sample(setdiff(seq_len(G), 2L), n, replace = TRUE)
  } else {
    c(seq_len(G), sample(seq_len(G), n - G, replace = TRUE))
  }
  Z <- matrix(0, n, G)
  Z[cbind(seq_len(n), grp)] <- 1
  W <- cbind(X, Z)
  list(W = W, X = X, grp = grp, G = G, px = px,
       nested = c(rep(FALSE, px), rep(TRUE, G)),
       pen = c(rep(0, px), rep(0.3, G)),
       wt = matrix(stats::runif(b * n, 0.2, 2), b, n))
}

.refS <- function(f) {
  sol <- .newtonSolver(f$W, f$pen, f$nested)
  out <- array(0, c(nrow(f$wt), f$px, f$px))
  for (g in seq_len(nrow(f$wt))) out[g, , ] <- sol$factor(f$wt[g, ])$S
  out
}

# S = A - B C^-1 B' written out, independent of .newtonSolver(). Needed for the
# empty-group case: .newtonSolver()'s rowsum() DROPS a group with no cells, so
# cvec comes back short and pen_z recycles against it. That is unreachable from
# the package's own designs -- .buildRandomEffects() (spiDE) builds the nested
# block with interaction(drop = TRUE), so every column holds at least one cell
# -- but .absorbBatch() gets the robustness free from .segmentSum() and should
# not quietly lose it.
.refS_direct <- function(f) {
  Z <- f$W[, f$nested, drop = FALSE]
  out <- array(0, c(nrow(f$wt), f$px, f$px))
  for (g in seq_len(nrow(f$wt))) {
    w <- f$wt[g, ]
    A <- crossprod(f$X * w, f$X) + diag(f$pen[!f$nested], f$px)
    Cd <- colSums(Z * w) + f$pen[f$nested]
    B <- crossprod(f$X * w, Z)
    out[g, , ] <- A - B %*% (t(B) / Cd)
  }
  out
}

test_that(".absorbBatch reproduces .newtonSolver()'s Schur complement per gene", {
  f <- .absorbFixture()
  got <- .absorbBatch(f$W, f$pen, f$nested, f$wt)
  expect_equal(dim(got), c(nrow(f$wt), f$px, f$px))
  expect_equal(got, .refS(f), tolerance = 1e-10)
  # and against the definition, so the two references corroborate each other
  expect_equal(got, .refS_direct(f), tolerance = 1e-10)
})

test_that(".absorbBatch keeps a group that holds no cells", {
  # compared against the written-out definition, NOT .newtonSolver(): the
  # oracle misaligns here, and that is a property of rowsum() rather than of
  # the absorption (see .refS_direct above)
  f <- .absorbFixture(empty_group = TRUE)
  expect_equal(.absorbBatch(f$W, f$pen, f$nested, f$wt),
               .refS_direct(f), tolerance = 1e-10)
})

test_that(".absorbBatch is invariant to the cell tile", {
  # the tile bounds the (batch, tile, px) intermediate; it is performance, not
  # semantics, so every tiling must give the same stack
  f <- .absorbFixture()
  whole <- .absorbBatch(f$W, f$pen, f$nested, f$wt)
  for (tile in c(1L, 7L, 59L, 60L, 1000L)) {
    expect_equal(.absorbBatch(f$W, f$pen, f$nested, f$wt,
                              cell.tile = tile),
                 whole, tolerance = 1e-12,
                 info = sprintf("cell.tile = %d", tile))
  }
})

test_that(".absorbBatch agrees between the base-R and torch branches", {
  skip_if_no_torch()
  f <- .absorbFixture()
  base_S <- .absorbBatch(f$W, f$pen, f$nested, f$wt)
  Wt <- torch::torch_tensor(f$W, dtype = torch::torch_float64())
  wtt <- torch::torch_tensor(f$wt, dtype = torch::torch_float64())
  tor_S <- .absorbBatch(Wt, f$pen, f$nested, wtt)
  expect_true(is_torch_tensor(tor_S))
  expect_equal(as.array(tor_S), base_S, tolerance = gpu_tol())
})

test_that(".absorbBatch refuses columns that are not a partition", {
  f <- .absorbFixture()
  # cell 1 is in group 1 by construction, so give it group 2 as well
  f$W[1, f$px + 2] <- 1
  expect_error(.absorbBatch(f$W, f$pen, f$nested, f$wt),
               "partition")
})
