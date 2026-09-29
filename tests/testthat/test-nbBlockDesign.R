# nbBlockDesign(): a design whose per-group block is held in compact form (one
# n x q matrix and a group per cell) and absorbed group by group, never built
# as the dense n x (p_x + G q) matrix.
#
# The oracle throughout is that dense matrix, as.matrix(D), with the block
# columns absorbed by group through the existing grouped path
# (.newtonSolver()/.newtonSolverBlocked(), polishNB(absorb = <grouping>)), and
# where it applies the full dense inverse with no absorption at all. The
# compact path is an exact re-arrangement of the same arithmetic, so the
# standard is agreement to rounding (1e-8 or tighter), not "close".

# A per-patient intercept and per-patient spatial library-size spline: the
# shape the compact form exists for. `sizes` are the patients' cell counts
# (0 gives a level with no cells); `ztiny` makes patient `tiny` hold fewer cells
# than there are block columns.
.blockFixture <- function(sizes = c(150, 60, 220, 35, 90), px = 3,
                          df = c(2, 2), ngenes = 5, seed = 1) {
  set.seed(seed)
  G <- length(sizes)
  lev <- sprintf("P%02d", seq_len(G))
  block <- factor(rep(lev, sizes), levels = lev)
  n <- length(block)
  perm <- sample.int(n)                       # cells not sorted by patient
  block <- block[perm]
  xy <- cbind(runif(n), runif(n))
  X <- matrix(rnorm(n * px), n, px, dimnames = list(NULL, paste0("n", seq_len(px))))
  l <- rnorm(n, sd = 0.5)                     # centred log library size
  Bs <- tpsBasis(xy[, 1], xy[, 2], df = df)
  Z <- cbind(int = 1, l * Bs[, ])
  colnames(Z) <- c("int", paste0("ls", seq_len(ncol(Bs))))
  D <- nbBlockDesign(X, Z, block)
  q <- ncol(Z)
  beta <- seq(0.3, -0.2, length.out = px)
  a_pat <- rnorm(G, 0, 0.4)
  mu <- exp(0.8 + as.numeric(X %*% beta) + a_pat[as.integer(block)] +
              l * (0.9 + 0.3 * Bs[, 1]))
  Y <- t(vapply(seq_len(ngenes), function(g)
    as.numeric(rnbinom(n, mu = mu * exp(0.2 * g), size = 3 + g)), numeric(n)))
  rownames(Y) <- paste0("g", seq_len(ngenes))
  list(D = D, W = as.matrix(D), group = .blockDenseAbsorb(D), X = X, Z = Z,
       block = block, Y = Y, n = n, px = px, q = q, G = G,
       A0 = matrix(0, ngenes, ncol(D), dimnames = list(rownames(Y), NULL)),
       psi0 = rep(0.4, ngenes),
       # the block intercept is the start column; a small ridge on the spline
       start = c(rep(FALSE, px), TRUE, rep(FALSE, q - 1L)),
       pen = c(rep(0, px), 0, rep(1e-3, q - 1L)))
}

# polishNB()'s outputs that must agree between the two paths
.expect_same_polish <- function(a, b, tol = 1e-8) {
  expect_equal(a$alpha, b$alpha, tolerance = tol)
  expect_equal(a$psi, b$psi, tolerance = tol)
  expect_equal(a$loglik, b$loglik, tolerance = tol)
  for (f in c("restarted", "capped", "singular", "psi_bound", "polished")) {
    expect_identical(a$polish[[f]], b$polish[[f]], info = f)
  }
}

test_that("the compact design is the dense block design it stands for", {
  f <- .blockFixture()
  D <- f$D
  expect_s3_class(D, "nbBlockDesign")
  expect_equal(dim(D), c(f$n, f$px + f$G * f$q))
  expect_equal(ncol(D), ncol(f$W))
  expect_identical(colnames(D)[seq_len(f$px + 2L)], c("n1", "n2", "n3", "P01:int", "P01:ls1"))
  # the dense expansion, written out independently of the class
  W <- matrix(0, f$n, f$px + f$G * f$q)
  W[, seq_len(f$px)] <- f$X
  for (g in seq_len(f$G)) {
    cells <- which(as.integer(f$block) == g)
    W[cells, f$px + (g - 1L) * f$q + seq_len(f$q)] <- f$Z[cells, ]
  }
  expect_equal(unname(f$W), W, tolerance = 0)
  expect_output(print(D), "5 groups x 5 block columns")

  # every product the engines take through the design helpers
  set.seed(9)
  a <- rnorm(ncol(D)); r <- rnorm(f$n)
  A <- matrix(rnorm(3 * ncol(D)), 3); R <- matrix(rnorm(3 * f$n), 3)
  expect_equal(.blockEtaVec(D, a), as.numeric(W %*% a), tolerance = 1e-12)
  expect_equal(.blockScoreVec(D, r), as.numeric(crossprod(W, r)), tolerance = 1e-12)
  expect_equal(.blockEtaBatch(D, A), A %*% t(W), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(.blockScoreBatch(R, D), R %*% W, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(.muBatch(A * 0.1, D), .muBatch(A * 0.1, W), tolerance = 1e-12,
               ignore_attr = TRUE)
  # a batch at or above BLOCK_BATCH_DENSE_MIN genes runs group by group
  A8 <- matrix(rnorm(8 * ncol(D)), 8); R8 <- matrix(rnorm(8 * f$n), 8)
  expect_gte(nrow(A8), BLOCK_BATCH_DENSE_MIN)
  expect_equal(.blockEtaBatch(D, A8), A8 %*% t(W), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(.blockScoreBatch(R8, D), R8 %*% W, tolerance = 1e-12,
               ignore_attr = TRUE)
  for (j in seq_len(ncol(D))) {
    expect_identical(.blockColCells(D, j), which(W[, j] != 0), info = j)
  }
  # a per-column vector for [X | Z] expands to every group's copy
  expect_equal(.blockExpand(D, c(1, 2, 3, 10 + seq_len(f$q)), "v"),
               c(1, 2, 3, rep(10 + seq_len(f$q), f$G)))
  expect_equal(.blockExpand(D, 7, "v"), rep(7, ncol(D)))
  expect_error(.blockExpand(D, 1:4, "lambda.a"), "lambda.a")
})

test_that("the compact solver's step, covariance and Schur pieces match the dense solvers", {
  f <- .blockFixture()
  D <- f$D
  pen <- .blockExpand(D, f$pen, "pen")
  sc <- .newtonSolverCompact(D, pen)
  sd <- .newtonSolver(f$W, pen, f$group)
  set.seed(4)
  for (rep in 1:3) {
    w <- runif(f$n, 0.05, 5)
    s <- rnorm(ncol(D))
    st <- sc$factor(w)
    sd_st <- sd$factor(w)
    expect_equal(st$rank, rep(f$q, f$G))
    # the Schur complement and the cross block are the dense grouped solver's
    expect_equal(st$S, sd_st$S, tolerance = 1e-10)
    expect_equal(st$B, sd_st$B, tolerance = 1e-10, ignore_attr = TRUE)
    # each group's inverse block, H H' = C_g^-1
    for (g in seq_len(f$G)) {
      expect_equal(tcrossprod(st$H[[g]]), chol2inv(sd_st$fac[[g]]),
                   tolerance = 1e-10)
    }
    # the step and the dense-column covariance, against the grouped solver and
    # against the full dense inverse
    info <- crossprod(f$W * sqrt(w))
    diag(info) <- diag(info) + pen
    expect_equal(sc$solve(st, s), sd$solve(sd_st, s), tolerance = 1e-10)
    expect_equal(sc$solve(w, s), as.numeric(solve(info, s)), tolerance = 1e-9)
    expect_equal(sc$xcov(st), sd$xcov(sd_st), tolerance = 1e-10)
    expect_equal(sc$xcov(w), solve(info)[seq_len(f$px), seq_len(f$px)],
                 tolerance = 1e-9)
  }
  # the exported solver dispatches on the class and expands [X | Z] penalties
  ex <- nbNewtonSolver(D, f$pen)
  w <- runif(f$n, 0.05, 5)
  expect_equal(ex$xcov(w), sc$xcov(w), tolerance = 0)
  expect_error(nbNewtonSolver(D, f$pen, absorb = f$group), "must be NULL")
  expect_no_error(nbNewtonSolver(D, f$pen, absorb = rep(FALSE, ncol(D))))
})

test_that("polishNB on a block design matches the dense grouped polish, both engines", {
  f <- .blockFixture()
  for (eng in c("batch", "gene")) {
    for (pm in c("profile", "fixed")) {
      a <- polishNB(f$Y, f$D, f$A0, f$psi0, lambda.a = f$pen,
                    start.cols = f$start, psi.method = pm, engine = eng,
                    tol = 1e-12)
      b <- polishNB(f$Y, f$W, f$A0, f$psi0,
                    lambda.a = .blockExpand(f$D, f$pen, "pen"),
                    absorb = f$group,
                    start.cols = .blockExpand(f$D, f$start, "start"),
                    psi.method = pm, engine = eng, tol = 1e-12)
      expect_true(all(a$polish$polished), info = paste(eng, pm))
      .expect_same_polish(a, b)
    }
  }
  # a converged fit, re-polished warm at a nearby penalty
  a0 <- polishNB(f$Y, f$D, f$A0, f$psi0, lambda.a = f$pen, start.cols = f$start)
  pen2 <- f$pen * 3
  a <- polishNB(f$Y, f$D, a0$alpha, a0$psi, lambda.a = pen2, warm = TRUE,
                tol = 1e-12)
  b <- polishNB(f$Y, f$W, a0$alpha, a0$psi, lambda.a = .blockExpand(f$D, pen2, "pen"),
                absorb = f$group, warm = TRUE, tol = 1e-12)
  .expect_same_polish(a, b)
  # and the profile dispersion at a fixed mean
  expect_equal(nbProfilePsi(f$Y, f$D, a0$alpha, a0$psi),
               nbProfilePsi(f$Y, f$W, a0$alpha, a0$psi), tolerance = 1e-10)
})

test_that("an offset reaches every linear predictor of the block design", {
  f <- .blockFixture(seed = 3)
  set.seed(30)
  off_v <- rnorm(f$n, 0, 0.3)
  off_m <- matrix(rnorm(nrow(f$Y) * f$n, 0, 0.2), nrow(f$Y), f$n)
  for (off in list(off_v, off_m)) {
    for (eng in c("batch", "gene")) {
      a <- polishNB(f$Y, f$D, f$A0, f$psi0, lambda.a = f$pen, offset = off,
                    start.cols = f$start, engine = eng, tol = 1e-12)
      b <- polishNB(f$Y, f$W, f$A0, f$psi0,
                    lambda.a = .blockExpand(f$D, f$pen, "pen"), offset = off,
                    absorb = f$group,
                    start.cols = .blockExpand(f$D, f$start, "start"),
                    engine = eng, tol = 1e-12)
      .expect_same_polish(a, b)
    }
  }
  a <- polishNB(f$Y, f$D, f$A0, f$psi0, lambda.a = f$pen, offset = off_v,
                start.cols = f$start)
  expect_equal(nbProfilePsi(f$Y, f$D, a$alpha, a$psi, offset = off_v),
               nbProfilePsi(f$Y, f$W, a$alpha, a$psi, offset = off_v),
               tolerance = 1e-10)
})

test_that("gene blocking and workers do not move the block design's answer", {
  f <- .blockFixture(ngenes = 6)
  one <- polishNB(f$Y, f$D, f$A0, rep(0.4, 6), lambda.a = f$pen,
                  start.cols = f$start)
  blk <- polishNB(f$Y, f$D, f$A0, rep(0.4, 6), lambda.a = f$pen,
                  start.cols = f$start, block.size = 2, batch.size = 1)
  expect_equal(blk$alpha, one$alpha, tolerance = 1e-12)
  expect_equal(blk$psi, one$psi, tolerance = 1e-12)
})

test_that("a group with no cells: penalised it matches the dense path, unpenalised it keeps its start", {
  f <- .blockFixture(sizes = c(150, 0, 220, 35, 90), seed = 5)
  expect_output(print(f$D), "1 empty")
  # a ridge on every block column: the dense grouped path is well defined too
  pen_all <- c(rep(0, f$px), rep(1e-3, f$q))
  for (eng in c("batch", "gene")) {
    a <- polishNB(f$Y, f$D, f$A0, f$psi0, lambda.a = pen_all,
                  start.cols = f$start, engine = eng, tol = 1e-12)
    b <- polishNB(f$Y, f$W, f$A0, f$psi0,
                  lambda.a = .blockExpand(f$D, pen_all, "pen"), absorb = f$group,
                  start.cols = .blockExpand(f$D, f$start, "start"),
                  engine = eng, tol = 1e-12)
    .expect_same_polish(a, b)
  }

  # the fixture's own penalty leaves the block intercept unpenalised, so the
  # empty group's intercept is not identified: the compact path holds it at
  # its start and fits the rest exactly as the design without that group
  empty <- f$px + f$q + seq_len(f$q)            # group 2's columns
  red_W <- f$W[, -empty]
  red_group <- f$group[-empty]
  red_pen <- .blockExpand(f$D, f$pen, "pen")[-empty]
  red_start <- .blockExpand(f$D, f$start, "start")[-empty]
  A0 <- f$A0
  A0[, empty] <- 0.7                           # a start the fit must not move
  for (eng in c("batch", "gene")) {
    a <- polishNB(f$Y, f$D, A0, f$psi0, lambda.a = f$pen,
                  start.cols = f$start, engine = eng, tol = 1e-12)
    expect_true(all(a$polish$polished))
    expect_equal(unname(a$alpha[, empty[1]]), rep(0.7, nrow(f$Y)))
    b <- polishNB(f$Y, red_W, A0[, -empty], f$psi0, lambda.a = red_pen,
                  absorb = red_group, start.cols = red_start, engine = eng,
                  tol = 1e-12)
    expect_equal(a$alpha[, -empty], b$alpha, tolerance = 1e-8)
    expect_equal(a$psi, b$psi, tolerance = 1e-8)
    expect_equal(a$loglik, b$loglik, tolerance = 1e-8)
  }
})

test_that("an aliased group is solved by a generalised inverse, matching the design without its aliased columns", {
  # patient 4 has 2 cells and 5 unpenalised block columns: rank 2
  f <- .blockFixture(sizes = c(150, 60, 220, 2, 90), seed = 6)
  pen0 <- rep(0, f$px + f$q)
  st <- .newtonSolverCompact(f$D, .blockExpand(f$D, pen0, "pen"))$factor(
    runif(f$n, 0.1, 2))
  expect_equal(st$rank, c(f$q, f$q, f$q, 2L, f$q))

  # the oracle: the dense design with the aliased group reduced to an
  # independent subset of its columns (its first two, on its two cells)
  g4 <- f$px + 3L * f$q + seq_len(f$q)
  cells4 <- which(as.integer(f$block) == 4L)
  expect_equal(qr(f$W[cells4, g4])$rank, 2L)
  drop <- g4[-(1:2)]
  red_W <- f$W[, -drop]
  red_group <- f$group[-drop]
  start_full <- .blockExpand(f$D, f$start, "start")
  for (eng in c("batch", "gene")) {
    a <- polishNB(f$Y, f$D, f$A0, f$psi0, lambda.a = pen0, start.cols = f$start,
                  engine = eng, tol = 1e-12)
    b <- polishNB(f$Y, red_W, f$A0[, -drop], f$psi0, lambda.a = 0,
                  absorb = red_group, start.cols = start_full[-drop],
                  engine = eng, tol = 1e-12)
    expect_true(all(a$polish$polished))
    expect_equal(a$alpha[, seq_len(f$px)], b$alpha[, seq_len(f$px)], tolerance = 1e-8)
    expect_equal(a$psi, b$psi, tolerance = 1e-8)
    expect_equal(a$loglik, b$loglik, tolerance = 1e-8)
    # the fitted means are the same, although patient 4's own coefficients
    # are one solution among many
    expect_equal(tcrossprod(a$alpha, f$W), tcrossprod(b$alpha, red_W),
                 tolerance = 1e-8)
    # the dense-column covariance at the converged fit
    mu <- exp(tcrossprod(a$alpha[1, ], f$W))[1, ]
    w <- mu / (1 + a$psi[1] * mu)
    expect_equal(nbNewtonSolver(f$D, pen0)$xcov(w),
                 nbNewtonSolver(red_W, rep(0, ncol(red_W)), red_group)$xcov(w),
                 tolerance = 1e-8)

    # the dense grouped path on the full design has no generalised inverse:
    # it flags every gene singular and keeps the input fit (it used to stop
    # the whole polish in sqrt(NULL), and the batch engine on a subscript)
    d <- polishNB(f$Y, f$W, f$A0, f$psi0, lambda.a = 0, absorb = f$group,
                  engine = eng)
    expect_true(all(d$polish$singular))
    expect_false(any(d$polish$polished))
    expect_equal(d$alpha, f$A0)
  }
})

test_that("a nearly aliased group matches the dense grouped path", {
  f <- .blockFixture(seed = 7)
  # patient 2's last spline column is its first one plus a 1e-6 perturbation,
  # unpenalised: a condition number of ~1e10 on the unit-diagonal scale, yet
  # every Cholesky pivot stays above rank.tol, so the block is solved exactly
  # on both paths
  Z <- f$Z
  cells <- which(as.integer(f$block) == 2L)
  Z[cells, f$q] <- Z[cells, 2] + 1e-6 * rnorm(length(cells))
  D <- nbBlockDesign(f$X, Z, f$block)
  W <- as.matrix(D)
  pen <- rep(0, f$px + f$q)
  w <- runif(f$n, 0.1, 2)
  st <- .newtonSolverCompact(D, .blockExpand(D, pen, "pen"))$factor(w)
  expect_equal(st$rank, rep(f$q, f$G))
  Cb <- crossprod(Z[cells, ] * sqrt(w[cells]))
  diag(Cb) <- diag(Cb) + pen[f$px + seq_len(f$q)]
  expect_gt(kappa(Cb / tcrossprod(sqrt(diag(Cb))), exact = TRUE), 1e8)
  for (eng in c("batch", "gene")) {
    a <- polishNB(f$Y, D, f$A0, f$psi0, lambda.a = pen, start.cols = f$start,
                  engine = eng, tol = 1e-12)
    b <- polishNB(f$Y, W, f$A0, f$psi0, lambda.a = .blockExpand(D, pen, "pen"),
                  absorb = f$group, start.cols = .blockExpand(D, f$start, "start"),
                  engine = eng, tol = 1e-12)
    expect_true(all(a$polish$polished))
    expect_equal(a$alpha[, seq_len(f$px)], b$alpha[, seq_len(f$px)], tolerance = 1e-8)
    expect_equal(a$psi, b$psi, tolerance = 1e-8)
    expect_equal(a$loglik, b$loglik, tolerance = 1e-8)
    expect_equal(tcrossprod(a$alpha, W), tcrossprod(b$alpha, W), tolerance = 1e-8)
    # the coefficients along the near-null direction are ill determined, so
    # the whole vector is held to a tolerance that allows for the conditioning
    # (measured: agreement to 1e-12 here)
    expect_equal(a$alpha, b$alpha, tolerance = 1e-6)
  }
})

test_that("polishNB refuses what the block design cannot take", {
  f <- .blockFixture(ngenes = 2)
  A0 <- f$A0[1:2, ]
  expect_error(polishNB(f$Y, f$D, A0, 0.4, absorb = f$group), "must be NULL")
  expect_error(polishNB(f$Y, f$D, A0, 0.4,
                        absorb.batch = rep(FALSE, ncol(f$D))), "must be NULL")
  expect_error(polishNB(f$Y, f$D, A0, 0.4, lambda.a = 1:4), "lambda.a")
  expect_error(polishNB(f$Y, f$D, A0, 0.4, start.cols = c(TRUE, FALSE)),
               "start.cols")
  expect_error(polishNB(f$Y, f$D, A0[, -1], 0.4), "one column per column of W")
  local_mocked_bindings(
    checkGPU = function(...) TRUE,
    .requireFloat64 = function(...) invisible(TRUE)
  )
  expect_error(polishNB(f$Y, f$D, A0, 0.4, backend = "gpu"), "CPU only")
})

test_that("the batched and device helpers refuse a block design", {
  f <- .blockFixture(ngenes = 1)
  wt <- matrix(1, 2, f$n)
  expect_error(nbNewtonSolverBatch(f$D, 0), "dense design")
  expect_error(nbGramBatch(f$D, wt), "dense design")
  expect_error(nbAbsorbGramBatch(f$D, 0, NULL, wt), "dense design")
})

test_that("nbBlockDesign validates its inputs", {
  X <- matrix(rnorm(20), 10)
  Z <- cbind(1, rnorm(10))
  b <- rep(c("a", "b"), 5)
  expect_s3_class(nbBlockDesign(X, Z, b), "nbBlockDesign")
  expect_s3_class(nbBlockDesign(X[, 1], Z, b), "nbBlockDesign")   # a vector X
  expect_error(nbBlockDesign(X, Z[-1, ], b), "one row per cell")
  expect_error(nbBlockDesign(X, Z, b[-1]), "one value per cell")
  expect_error(nbBlockDesign(X, Z, replace(b, 3, NA)), "missing")
  expect_error(nbBlockDesign(replace(X, 2, NA), Z, b), "finite")
  expect_error(nbBlockDesign(X, Z[, 0, drop = FALSE], b), "at least one column")
  expect_error(nbBlockDesign(X, Z, b, rank.tol = 2), "rank.tol")
  expect_error(nbBlockDesign(X, matrix("a", 10, 1), b), "numeric")
})
