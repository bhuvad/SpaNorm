# Tests for fitNB()/calculateMu()'s fixed `offset` argument.
#
# The offset enters the linear predictor with a coefficient pinned at 1:
#   log(mu) = gmean + tcrossprod(alpha, W) + offset
# The two properties worth proving are (a) that the coefficient really is 1 --
# it is a known effect, not an estimated one, so it cannot absorb signal that
# happens to correlate with it -- and (b) that the genes x cells offset matrix
# stays aligned with the counts when the fit is carved into gene-blocks. A
# misaligned offset does not error and does not fail to converge; it converges
# to the wrong answer, so "it ran" is not evidence here.
#
# The exact-equivalence construction used below: fitting design
#   W_A = [(1+c)*z, 1, x]                  with no offset
# and design
#   W_B = [z, 1, x]                        with offset c*z
# are the same model written two ways -- the offset lies in the span of the
# design, so the reachable linear predictors are identical. Both fits also
# start from the same alpha init ([1, 0, ...]), hence the same initial log-mean,
# and the coefficient map between them is affine with a positive scale, so
# winsoriseCols() commutes with it. Every step of both fits therefore coincides
# exactly: same log-likelihood, same dispersions, same step-halving, same fitted
# means -- while the coefficients differ by exactly the amount the pinned unit
# coefficient accounts for. Any deviation from a coefficient of exactly 1, in
# any one of the four places the offset enters (log-likelihood, IRLS working
# response, gmean fold, edgeR dispersion offset), breaks the identity.
#
# Two details of that construction are load-bearing, and getting either wrong
# yields a test that passes on broken code:
#   * the scaled column must be COLUMN 1, the only one the fit initialises to a
#     coefficient of 1. Scale any other column and the two fits start from
#     different log-means, the trajectories diverge, and the identity holds only
#     approximately (if at all).
#   * z must VARY ACROSS CELLS. edgeR::estimateDisp() is given an intercept-only
#     design, so it absorbs any offset that is constant within a gene -- under a
#     constant offset the dispersion step is genuinely indifferent to whether it
#     receives the offset at all, and an assertion on psi cannot fail. A
#     cell-varying offset is also the representative case: a library-size linear
#     predictor, the motivating use, varies per cell by construction.

quiet <- function(...) NULL

test_that("offset = NULL and an all-zero offset are strict no-ops", {
  set.seed(101)
  ng <- 12; ns <- 80
  W <- cbind(intercept = 1, cov = as.numeric(scale(rnorm(ns))))
  Y <- matrix(rpois(ng * ns, 8), ng, ns)

  # the seed is reset before each fit: the dispersion subsample is drawn with
  # sample.int(), so an unmatched draw would confound every comparison here
  set.seed(7); f.null <- fitNB(Y, W, maxit.psi = 3, backend = "cpu", verbose = FALSE)
  set.seed(7); f.zero <- fitNB(Y, W, maxit.psi = 3, offset = matrix(0, ng, ns),
                               backend = "cpu", verbose = FALSE)

  expect_identical(f.zero$alpha, f.null$alpha)
  expect_identical(f.zero$psi, f.null$psi)
  expect_identical(f.zero$loglik, f.null$loglik)

  # calculateMu: NULL leaves the pre-existing code path untouched
  expect_identical(calculateMu(f.null$gmean, f.null$alpha, W, offset = NULL),
                   calculateMu(f.null$gmean, f.null$alpha, W))
})

test_that("an offset is exactly a design column whose coefficient is pinned at 1", {
  set.seed(102)
  ng <- 12; ns <- 80
  cc <- 0.5
  # the scaled column MUST be column 1 (the only one the fit initialises to a
  # coefficient of 1), so that both fits start from the same log-mean; and it
  # must vary ACROSS CELLS, or the dispersion-offset check below is vacuous --
  # edgeR fits an intercept per gene, which silently absorbs any offset that is
  # constant within a gene, so a constant offset cannot detect a dispersion
  # step that ignored it.
  z <- as.numeric(scale(rnorm(ns)))
  x <- as.numeric(scale(rnorm(ns)))
  W.B <- cbind(z = z, intercept = 1, x = x)             # + offset cc*z
  W.A <- cbind(z = (1 + cc) * z, intercept = 1, x = x)  # same model, coefficient free
  a.true <- cbind(rnorm(ng, 0.3, 0.15), rnorm(ng, 1.2, 0.3), rnorm(ng, 0, 0.5))
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(exp(tcrossprod(a.true, W.B))), size = 5), ng, ns)
  O <- matrix(rep(cc * z, each = ng), ng, ns)

  set.seed(23); fA <- fitNB(Y, W.A, maxit.psi = 4, backend = "cpu", verbose = FALSE)
  set.seed(23); fB <- fitNB(Y, W.B, maxit.psi = 4, offset = O, backend = "cpu", verbose = FALSE)

  # guard against a vacuous identity: two fits both frozen at the shared
  # alpha initialisation would satisfy every equality below trivially
  expect_gt(max(abs(fA$alpha - matrix(c(1, 0, 0), ng, 3, byrow = TRUE))), 0.5)

  # coefficient identity: alpha_A[,1] * (1+c) == alpha_B[,1] + c*1
  # (i.e. the offset contributed exactly 1 x c*z, not 0, not 2c*z)
  expect_lt(max(abs(fA$alpha[, 1] * (1 + cc) - (fB$alpha[, 1] + cc))), 1e-10)
  expect_lt(max(abs(fA$alpha[, 2:3] - fB$alpha[, 2:3])), 1e-10)

  # dispersions must match too: the offset also enters edgeR::estimateDisp's
  # own offset. Were it dropped there, the coefficients would still look
  # plausible while every gene's psi was wrong -- and the Wald variances that
  # downstream packages compute from psi with them.
  expect_lt(max(abs(fA$psi - fB$psi)), 1e-10)
  expect_equal(fA$loglik, fB$loglik, tolerance = 1e-10)

  # fitted means, through calculateMu's own offset argument
  muA <- calculateMu(fA$gmean, fA$alpha, W.A, winsor = Inf)
  muB <- calculateMu(fB$gmean, fB$alpha, W.B, winsor = Inf, offset = O)
  expect_lt(max(abs(muA - muB) / muA), 1e-10)

  # and the working-weight Wald standard errors built from them (the quantity
  # downstream inference actually consumes)
  seOf <- function(mu, psi, W) {
    t(vapply(seq_len(nrow(mu)), function(g) {
      w <- mu[g, ] / (1 + psi[g] * mu[g, ])
      sqrt(diag(solve(crossprod(W, W * w))))
    }, numeric(ncol(W))))
  }
  seA <- seOf(muA, fA$psi, W.A)
  seB <- seOf(muB, fB$psi, W.B)
  expect_lt(max(abs(seA[, 2:3] - seB[, 2:3]) / seB[, 2:3]), 1e-10)        # shared columns
  expect_lt(max(abs(seA[, 1] * (1 + cc) - seB[, 1]) / seB[, 1]), 1e-10)   # rescaled column
})

test_that("a gene-varying offset shifts the coefficients by exactly its own contribution", {
  # As above but with a per-gene offset (offset = B %*% t(W)), which an offset
  # shared across genes cannot probe: here every gene gets a different row, so
  # this is what catches an offset that is correctly shaped but wrongly indexed
  # down the gene axis. psi is held fixed and the two fits are started at inits
  # differing by exactly B, so the identity alpha_offset == alpha_plain - B must
  # hold at every iteration, not just in the limit.
  set.seed(103)
  ng <- 12; ns <- 80
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  Y <- matrix(rpois(ng * ns, 9), ng, ns)
  psi <- rep(0.2, ng)
  A0 <- cbind(rnorm(ng, 1, 0.2), rnorm(ng, 0, 0.2))
  B <- cbind(rnorm(ng, 0, 0.4), rnorm(ng, 0, 0.4))
  O <- tcrossprod(B, W)

  gA <- fitNBGivenPsi(Y, W, psi, lambda.a = 0, gmean = rep(0, ng), alpha = A0,
                      winsor = Inf, backend = "cpu", msgfun = quiet)
  gB <- fitNBGivenPsi(Y, W, psi, lambda.a = 0, gmean = rep(0, ng), alpha = A0 - B,
                      offset = O, winsor = Inf, backend = "cpu", msgfun = quiet)

  # guard against a vacuous identity (both fits frozen at their inits, which
  # already differ by exactly B): the fit must move, and the counts here are
  # Poisson(9) with no covariate structure, so it must move to log(9)
  expect_gt(max(abs(gA$alpha - A0)), 0.5)
  expect_lt(abs(mean(gA$alpha[, 1]) - log(9)), 0.05)

  expect_lt(max(abs(gB$alpha - (gA$alpha - B))), 1e-10)
  expect_equal(gA$loglik, gB$loglik, tolerance = 1e-10)
})

test_that("a gene-blocked fit with an offset is identical to a single-block fit", {
  # THE alignment test. geneBlockCount() always returns 1L on CPU (that is what
  # keeps the CPU path free of blocking), so the partition is supplied directly
  # here rather than forced through a memory budget; the GPU test below covers
  # the budget-driven route. The offset varies along BOTH axes, with a different
  # magnitude and sign per gene, so a misalignment along either axis moves the
  # answer far more than the tolerance -- as the two controls at the end confirm.
  set.seed(104)
  ng <- 23; ns <- 90 # gene count deliberately not a multiple of the block counts
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  a.true <- cbind(rnorm(ng, 1.2, 0.3), rnorm(ng, 0, 0.5))
  O <- outer(rnorm(ng, 0, 0.6), cos(seq_len(ns) / 7))
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(exp(tcrossprod(a.true, W) + O)), size = 5), ng, ns)
  psi <- rep(0.2, ng)
  A0 <- cbind(rnorm(ng, 1, 0.2), rnorm(ng, 0, 0.2))

  expect_equal(geneBlockCount(ng, ns, ncol(W), backend = "cpu"), 1L)

  fitBlocked <- function(nb, off = O) {
    fitNBGivenPsi(Y, W, psi, lambda.a = 0, gmean = rep(0, ng), alpha = A0, offset = off,
                  maxit.nb = 12, blocks = geneBlockIndices(ng, nb),
                  backend = "cpu", msgfun = quiet)
  }
  ref <- fitBlocked(1)

  # guard against a vacuous comparison: if the fit never moved off its
  # initialisation, every block count would agree trivially and this test would
  # prove nothing. It must move, and toward the generating coefficients.
  expect_gt(max(abs(ref$alpha - A0)), 0.5)
  expect_lt(max(abs(ref$alpha - a.true)), 0.5)

  for (nb in c(2, 5)) {
    blocks <- geneBlockIndices(ng, nb)
    expect_gt(length(blocks), 1L) # sanity: the partition really is split
    expect_lt(max(abs(fitBlocked(nb)$alpha - ref$alpha)), 1e-10)
  }

  # the two blocked reductions that consume the offset, against literal
  # whole-matrix references. These are single-pass, so they also cover the
  # one-gene-per-block edge (where a dropped matrix dimension would bite)
  # without paying for a full IRLS fit at that block count.
  lmu <- tcrossprod(A0, W) + O
  ll.ref <- colSums(dnbinom(Y, mu = exp(lmu), size = 1 / psi, log = TRUE))
  for (nb in c(4, ng)) {
    blocks <- geneBlockIndices(ng, nb)
    expect_lt(max(abs(.blockedNBLoglik(A0, rep(0, ng), Y, W, psi, "cpu", blocks, O) - ll.ref)), 1e-8)
    expect_lt(max(abs(.blockedOffsPsi(A0, W, blocks, O) - lmu)), 1e-12)
  }

  # controls: the agreement above is only meaningful if a misaligned offset
  # would have been caught. Rotating the offset by one gene, or by one cell,
  # must move the fit -- by far more than the tolerances asserted above.
  expect_gt(max(abs(fitBlocked(1, O[c(2:ng, 1), ])$alpha - ref$alpha)), 1e-3)
  expect_gt(max(abs(fitBlocked(1, O[, c(2:ns, 1)])$alpha - ref$alpha)), 1e-3)
})

test_that("the offset is subset by idx in lockstep with the counts", {
  set.seed(105)
  ng <- 12; ns <- 80
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  Y <- matrix(rpois(ng * ns, 8), ng, ns)
  O <- outer(seq_len(ng) / ng, cos(seq_len(ns) / 5))
  idx <- rep(c(TRUE, TRUE, FALSE, TRUE), length.out = ns)

  # fitting the full data with idx, and fitting the pre-subset data with no
  # idx, must agree -- which pins down *which* offset columns were used
  set.seed(31)
  u1 <- fitSpaNormNB(Y, W, idx, offset = O, maxit.psi = 3, lambda.a = 0,
                     is.spanorm = FALSE, backend = "cpu", msgfun = quiet)
  set.seed(31)
  u2 <- fitSpaNormNB(Y[, idx, drop = FALSE], W[idx, , drop = FALSE], rep(TRUE, sum(idx)),
                     offset = O[, idx, drop = FALSE], maxit.psi = 3, lambda.a = 0,
                     is.spanorm = FALSE, backend = "cpu", msgfun = quiet)
  expect_gt(max(abs(u1$alpha - matrix(c(1, 0), ng, 2, byrow = TRUE))), 0.5) # not frozen at init
  expect_lt(max(abs(u1$alpha - u2$alpha)), 1e-10)
  expect_lt(max(abs(u1$psi - u2$psi)), 1e-10)

  # control: taking the first sum(idx) columns instead of the idx ones -- the
  # classic off-by-subset bug -- must not pass
  set.seed(31)
  u3 <- fitSpaNormNB(Y[, idx, drop = FALSE], W[idx, , drop = FALSE], rep(TRUE, sum(idx)),
                     offset = O[, seq_len(sum(idx)), drop = FALSE], maxit.psi = 3,
                     lambda.a = 0, is.spanorm = FALSE, backend = "cpu", msgfun = quiet)
  expect_gt(max(abs(u3$alpha - u1$alpha)), 1e-3)
})

test_that("an out-of-span offset is applied at unit magnitude, not omitted or doubled", {
  # The tests above all place the offset inside the span of the design, where it
  # is exactly checkable. Here it is outside the span, so it genuinely alters
  # the answer: the coefficients are recovered only when the offset is supplied,
  # at exactly its true magnitude. Omitting it, or supplying twice it, biases
  # the estimates -- which is the property that makes an offset the right tool
  # for a known effect a fitted column would otherwise partly absorb.
  set.seed(106)
  ng <- 15; ns <- 120
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  a.true <- cbind(rnorm(ng, 1.2, 0.3), rnorm(ng, 0, 0.5))
  v <- as.numeric(scale(rnorm(ns)))     # not in span(W)
  O <- outer(runif(ng, 0.8, 1.6), v)    # gene-varying magnitude
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(exp(tcrossprod(a.true, W) + O)), size = 8), ng, ns)

  set.seed(41); f.on  <- fitNB(Y, W, maxit.psi = 4, offset = O,     backend = "cpu", verbose = FALSE)
  set.seed(41); f.off <- fitNB(Y, W, maxit.psi = 4,                 backend = "cpu", verbose = FALSE)
  set.seed(41); f.x2  <- fitNB(Y, W, maxit.psi = 4, offset = 2 * O, backend = "cpu", verbose = FALSE)

  rmse <- function(f) sqrt(mean((f$alpha - a.true)^2))
  expect_lt(rmse(f.on), 0.15)
  expect_gt(rmse(f.off), 3 * rmse(f.on))
  expect_gt(rmse(f.x2), 2 * rmse(f.on))
})

test_that("the offset holds up in the small-n, few-iteration regime", {
  # the regime downstream subset fits actually run in: tens of cells,
  # lambda.a = 0, winsor at its default (not disabled), maxit.psi = 2. The
  # offset is cell-varying here for the same reason as in the equivalence test
  # above -- that is both the representative case (a library-size linear
  # predictor varies per cell by construction) and the only one under which the
  # psi assertion has any force.
  set.seed(107)
  ng <- 12; ns <- 30
  cc <- 0.5
  z <- as.numeric(scale(rnorm(ns)))
  x <- as.numeric(scale(rnorm(ns)))
  W.B <- cbind(z = z, intercept = 1, x = x)
  W.A <- cbind(z = (1 + cc) * z, intercept = 1, x = x)
  a.true <- cbind(rnorm(ng, 0.3, 0.15), rnorm(ng, 1.2, 0.3), rnorm(ng, 0, 0.5))
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(exp(tcrossprod(a.true, W.B))), size = 5), ng, ns)
  O <- matrix(rep(cc * z, each = ng), ng, ns)

  set.seed(51); sA <- fitNB(Y, W.A, lambda.a = 0, winsor = 4, maxit.psi = 2,
                            backend = "cpu", verbose = FALSE)
  set.seed(51); sB <- fitNB(Y, W.B, lambda.a = 0, winsor = 4, maxit.psi = 2,
                            offset = O, backend = "cpu", verbose = FALSE)

  expect_gt(max(abs(sA$alpha - matrix(c(1, 0, 0), ng, 3, byrow = TRUE))), 0.5) # not frozen
  expect_lt(max(abs(sA$alpha[, 1] * (1 + cc) - (sB$alpha[, 1] + cc))), 1e-10)
  expect_lt(max(abs(sA$alpha[, 2:3] - sB$alpha[, 2:3])), 1e-10)
  expect_lt(max(abs(sA$psi - sB$psi)), 1e-10)
})

test_that("the is.spanorm model (profiled gmean) carries an offset too", {
  # is.spanorm = TRUE takes a different path through the inner iteration: a
  # profiled per-gene intercept updated by the gmean fold, and a shared
  # unpenalised first column. The offset has to come off the working response
  # in the fold as well, or gmean silently absorbs it.
  set.seed(108)
  ng <- 12; ns <- 80
  cc <- 1.5
  Y <- matrix(rpois(ng * ns, 10), ng, ns)
  # centred and spread out: the gmean fold weights the working response by
  # sig.inv, so an offset the fold failed to subtract would leak into gmean by
  # its *weighted* mean -- zero only if the offset were constant across cells
  logLS <- log(colSums(Y)); logLS <- (logLS - mean(logLS)) * 5
  x <- as.numeric(scale(rnorm(ns)))
  W.B <- cbind(ls = logLS, x = x)
  W.A <- cbind(ls = (1 + cc) * logLS, x = x)
  O <- matrix(rep(cc * logLS, each = ng), ng, ns)

  set.seed(61)
  sA <- fitSpaNormNB(Y, W.A, rep(TRUE, ns), maxit.psi = 2, is.spanorm = TRUE,
                     lambda.a = 0.5, winsor = Inf, backend = "cpu", msgfun = quiet)
  set.seed(61)
  sB <- fitSpaNormNB(Y, W.B, rep(TRUE, ns), maxit.psi = 2, is.spanorm = TRUE,
                     lambda.a = 0.5, winsor = Inf, offset = O, backend = "cpu", msgfun = quiet)

  expect_gt(max(abs(sA$alpha - matrix(c(1, 0), ng, 2, byrow = TRUE))), 0.5) # not frozen at init
  expect_lt(max(abs(sA$gmean - sB$gmean)), 1e-8)
  expect_lt(max(abs(sA$alpha[, 1] * (1 + cc) - (sB$alpha[, 1] + cc))), 1e-8)
  expect_lt(max(abs(sA$alpha[, 2] - sB$alpha[, 2])), 1e-8)
})

test_that("calculateMu adds the offset exactly where gmean goes", {
  set.seed(109)
  ng <- 10; ns <- 40
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  alpha <- cbind(rnorm(ng, 1, 0.3), rnorm(ng, 0, 0.5))
  gmean <- rnorm(ng, 0, 0.2)
  O <- outer(runif(ng), sin(seq_len(ns)))

  lmu <- gmean + tcrossprod(alpha, W) + O
  expect_equal(calculateMu(gmean, alpha, W, winsor = Inf, offset = O), exp(lmu))

  # winsorisation clamps the TOTAL log-mean, i.e. after the offset is added --
  # that is what keeps it bounding the fitted mean rather than a component of it
  clamped <- exp(pmin(lmu, matrixStats::rowMedians(lmu) + 4 * matrixStats::rowMads(lmu)))
  expect_equal(calculateMu(gmean, alpha, W, winsor = 4, offset = O), clamped)
})

test_that("a malformed offset is rejected", {
  set.seed(110)
  ng <- 8; ns <- 30
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  alpha <- matrix(0, ng, 2)
  Y <- matrix(rpois(ng * ns, 5), ng, ns)

  expect_error(fitNB(Y, W, offset = matrix(0, ng, ns - 1), backend = "cpu", verbose = FALSE),
               "genes x cells")
  expect_error(fitNB(Y, W, offset = matrix(0, ng + 1, ns), backend = "cpu", verbose = FALSE),
               "genes x cells")
  expect_error(fitNB(Y, W, offset = rep(0, ns), backend = "cpu", verbose = FALSE),
               "genes x cells")
  expect_error(fitNB(Y, W, offset = matrix(Inf, ng, ns), backend = "cpu", verbose = FALSE),
               "finite")
  expect_error(fitNB(Y, W, offset = matrix(NA_real_, ng, ns), backend = "cpu", verbose = FALSE),
               "finite")
  expect_error(calculateMu(rep(0, ng), alpha, W, offset = matrix(0, 2, 2)), "genes x cells")
})

test_that("a forced tiny gpu.mem.budget keeps an offset aligned across blocks", {
  # the budget-driven route into multi-block fitting, which only engages on an
  # accelerator. UNVERIFIED as of writing: no GPU was available when the offset
  # support was added, so this test has never been executed -- it is here so
  # that the first run on a real device checks the offset alongside everything
  # else. Use an h100-class device: fp64 digamma inside the dispersion step is
  # what breaks on lesser cards, and it does not fail at startup.
  skip_if_no_gpu()
  resetGPUCache()
  on.exit(resetGPUCache())

  set.seed(111)
  ng <- 37; ns <- 60 # not a multiple of the forced block count
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  a.true <- cbind(rnorm(ng, 1.2, 0.3), rnorm(ng, 0, 0.5))
  O <- outer(exp(seq(-2, 2, length.out = ng)), cos(seq_len(ns) / 7))
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(exp(tcrossprod(a.true, W) + O)), size = 5), ng, ns)

  set.seed(71); fit.unblocked <- fitNB(Y, W, offset = O, backend = "gpu",
                                       gpu.mem.budget = Inf, verbose = FALSE)
  tiny.budget <- tinyGpuBudget(ng, ns)
  expect_gt(geneBlockCount(ng, ns, ncol(W), backend = "gpu", budget.bytes = tiny.budget), 1L)
  set.seed(71); fit.blocked <- fitNB(Y, W, offset = O, backend = "gpu",
                                     gpu.mem.budget = tiny.budget, verbose = FALSE)

  expect_equal(fit.blocked$alpha, fit.unblocked$alpha, tolerance = gpu_tol())
  expect_equal(fit.blocked$psi, fit.unblocked$psi, tolerance = gpu_tol())

  set.seed(71); fit.cpu <- fitNB(Y, W, offset = O, backend = "cpu", verbose = FALSE)
  expect_equal(fit.blocked$alpha, fit.cpu$alpha, tolerance = 1e-3)
})
