# Tests for fitNB()'s supplied-dispersion (`psi`) passthrough.
#
# With psi supplied, the outer dispersion loop is bypassed entirely: no
# edgeR::estimateDisp call, no dispersion subsample, one inner IRLS fit at the
# given psi, and the supplied values are returned as-is (no re-estimation, no
# winsorisation -- the caller owns them). The properties worth proving:
#   (a) the supplied psi is the one the fit actually USED, not merely the one
#       echoed back -- a wiring bug that validates and returns psi but fits at
#       something else produces plausible coefficients and a green "psi is
#       returned" check, so the assertions here pin the fit to psi through the
#       log-likelihood identity and through fits at different psi diverging;
#   (b) psi = NULL is a strict no-op; and
#   (c) psi composes with `offset`, the primary downstream use.

quiet <- function(...) NULL

test_that("psi = NULL is a strict no-op", {
  set.seed(201)
  ng <- 12; ns <- 80
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  Y <- matrix(rpois(ng * ns, 8), ng, ns)

  # seed reset before each fit: the psi = NULL path draws a dispersion
  # subsample with sample.int(), so an unmatched draw would confound this
  set.seed(9); f0 <- fitNB(Y, W, maxit.psi = 3, backend = "cpu", verbose = FALSE)
  set.seed(9); fn <- fitNB(Y, W, maxit.psi = 3, psi = NULL, backend = "cpu", verbose = FALSE)

  expect_identical(fn$alpha, f0$alpha)
  expect_identical(fn$psi, f0$psi)
  expect_identical(fn$loglik, f0$loglik)
  expect_identical(fn$sampling, f0$sampling)
})

test_that("a supplied psi is used as-is: returned unchanged, fitted at, never re-estimated", {
  set.seed(202)
  ng <- 15; ns <- 120
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  a.true <- cbind(rnorm(ng, 1.2, 0.3), rnorm(ng, 0, 0.5))
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(exp(tcrossprod(a.true, W))), size = 5), ng, ns)
  p0 <- runif(ng, 0.05, 0.6)

  fp <- fitNB(Y, W, psi = p0, backend = "cpu", verbose = FALSE)
  fp2 <- fitNB(Y, W, psi = 2 * p0, backend = "cpu", verbose = FALSE)

  # echoed back exactly (as-is: no re-estimation, no winsorisation)
  expect_identical(fp$psi, as.numeric(p0))
  expect_identical(fp2$psi, as.numeric(2 * p0))

  # not frozen at the init, and genuinely fitted
  expect_gt(max(abs(fp$alpha - matrix(c(1, 0), ng, 2, byrow = TRUE))), 0.5)

  # bypass really taken: single-element loglik, and no dispersion subsample
  # was drawn, so the sampling factor has no "dispersion" level
  expect_length(fp$loglik, 1L)
  expect_false("dispersion" %in% levels(fp$sampling))

  # the log-likelihood identity, which pins WHICH psi the fit ran at: the
  # returned loglik must equal the literal NB log-likelihood at the returned
  # coefficients and the SUPPLIED psi. A wiring bug that fits at some other
  # psi returns a loglik computed at that other psi, and this breaks.
  ll.lit <- sum(dnbinom(Y, mu = exp(tcrossprod(fp$alpha, W)), size = 1 / p0, log = TRUE))
  expect_lt(abs(fp$loglik - ll.lit), 1e-6)

  # fits at p0 and 2*p0 must DIFFER: the IRLS weights w = mu/(1 + psi*mu)
  # depend on psi, so identical alphas mean the supplied psi never reached the
  # fit (measured difference ~3.5e-3 on this fixture)
  expect_gt(max(abs(fp$alpha - fp2$alpha)), 1e-3)

  # and the difference has the predictable direction: w is decreasing in psi,
  # so the working-weight Wald SE at 2*p0 exceeds the SE at p0 for every gene
  # and coefficient -- not just on average
  seOf <- function(alpha, psi) {
    mu <- exp(tcrossprod(alpha, W))
    t(vapply(seq_len(ng), function(g) {
      w <- mu[g, ] / (1 + psi[g] * mu[g, ])
      sqrt(diag(solve(crossprod(W, W * w))))
    }, numeric(ncol(W))))
  }
  expect_true(all(seOf(fp2$alpha, fp2$psi) > seOf(fp$alpha, fp$psi)))

  # estimateDisp is never called on the supplied-psi path
  called <- FALSE
  suppressMessages(trace(edgeR::estimateDisp,
                         tracer = function() called <<- TRUE, print = FALSE))
  on.exit(suppressMessages(untrace(edgeR::estimateDisp)), add = TRUE)
  invisible(fitNB(Y, W, psi = p0, backend = "cpu", verbose = FALSE))
  expect_false(called)
  # positive control, without which the expect_false above is vacuous (a trace
  # that silently failed to install would also never set `called`): the same
  # trace must fire on the estimating path
  invisible(fitNB(Y, W, maxit.psi = 1, backend = "cpu", verbose = FALSE))
  expect_true(called)
})

test_that("refitting at the estimated psi reproduces the original fit", {
  set.seed(203)
  ng <- 15; ns <- 120
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  a.true <- cbind(rnorm(ng, 1.2, 0.3), rnorm(ng, 0, 0.5))
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(exp(tcrossprod(a.true, W))), size = 5), ng, ns)

  set.seed(9); fa <- fitNB(Y, W, backend = "cpu", verbose = FALSE)
  fb <- fitNB(Y, W, psi = fa$psi, backend = "cpu", verbose = FALSE)

  expect_identical(fb$psi, fa$psi)
  # Same objective, same psi, so both fits approach the same optimum -- but
  # not to machine precision: the first fit's final inner IRLS starts warm
  # (from the previous outer iteration's alpha) and the refit starts cold, and
  # both stop at the inner loop's relative-loglik tolerance of 1e-4, which is
  # not reachable from fitNB()'s arguments (fitNB's `tol` binds the OUTER
  # loop's formal). The measured gap on this fixture is ~4e-3 of pure
  # optimiser stopping slack; 0.02 gives a 4x margin. The psi-identity and
  # loglik-identity assertions above -- not this tolerance -- carry the duty
  # of catching a psi that was ignored or rewired.
  expect_lt(max(abs(fa$alpha - fb$alpha)), 0.02)
  expect_lt(abs(fa$loglik[length(fa$loglik)] - fb$loglik) / abs(fb$loglik), 1e-4)
})

test_that("psi and offset compose: pinned-column equivalence holds at a supplied psi", {
  # the primary downstream use: raw-count fit with the normalisation linear
  # predictor as a fixed offset AND a pooled, pre-estimated dispersion. Same
  # construction as test-fitNB-offset.R (cell-varying offset in the span of
  # the design, scaled column first so both fits share their initialisation);
  # with psi supplied to BOTH fits the trajectories must again coincide
  # exactly, now through the bypass path.
  set.seed(204)
  ng <- 12; ns <- 80
  cc <- 0.5
  z <- as.numeric(scale(rnorm(ns)))
  x <- as.numeric(scale(rnorm(ns)))
  W.B <- cbind(z = z, intercept = 1, x = x)
  W.A <- cbind(z = (1 + cc) * z, intercept = 1, x = x)
  a.true <- cbind(rnorm(ng, 0.3, 0.15), rnorm(ng, 1.2, 0.3), rnorm(ng, 0, 0.5))
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(exp(tcrossprod(a.true, W.B))), size = 5), ng, ns)
  O <- matrix(rep(cc * z, each = ng), ng, ns)
  p0 <- runif(ng, 0.05, 0.6)

  hA <- fitNB(Y, W.A, psi = p0, backend = "cpu", verbose = FALSE)
  hB <- fitNB(Y, W.B, psi = p0, offset = O, backend = "cpu", verbose = FALSE)

  expect_gt(max(abs(hA$alpha - matrix(c(1, 0, 0), ng, 3, byrow = TRUE))), 0.5) # not frozen
  expect_lt(max(abs(hA$alpha[, 1] * (1 + cc) - (hB$alpha[, 1] + cc))), 1e-10)
  expect_lt(max(abs(hA$alpha[, 2:3] - hB$alpha[, 2:3])), 1e-10)
  expect_lt(abs(hA$loglik - hB$loglik), 1e-8)
  expect_identical(hB$psi, as.numeric(p0))
})

test_that("psi is subset-safe: idx plus a supplied psi matches the pre-subset fit", {
  # psi is per-gene, so unlike the offset it needs no idx subsetting -- but the
  # combination has to route through the bypass identically either way
  set.seed(205)
  ng <- 12; ns <- 90
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  Y <- matrix(rpois(ng * ns, 8), ng, ns)
  p0 <- runif(ng, 0.05, 0.6)
  idx <- rep(c(TRUE, TRUE, FALSE), length.out = ns)

  u1 <- fitNB(Y, W, idx = idx, psi = p0, backend = "cpu", verbose = FALSE)
  u2 <- fitNB(Y[, idx, drop = FALSE], W[idx, , drop = FALSE], psi = p0,
              backend = "cpu", verbose = FALSE)

  expect_gt(max(abs(u1$alpha - matrix(c(1, 0), ng, 2, byrow = TRUE))), 0.5) # not frozen
  expect_identical(u1$alpha, u2$alpha)
  # cells not selected by idx are exported as "all", fitted cells as "glm"
  expect_identical(as.vector(table(u1$sampling)[c("all", "glm")]),
                   c(sum(!idx), sum(idx)))
})

test_that("a malformed psi is rejected", {
  set.seed(206)
  ng <- 8; ns <- 30
  W <- cbind(intercept = 1, x = as.numeric(scale(rnorm(ns))))
  Y <- matrix(rpois(ng * ns, 5), ng, ns)
  p0 <- runif(ng, 0.05, 0.6)

  expect_error(fitNB(Y, W, psi = p0[-1], backend = "cpu", verbose = FALSE),
               "length 8")
  expect_error(fitNB(Y, W, psi = letters[1:ng], backend = "cpu", verbose = FALSE),
               "length 8")
  expect_error(fitNB(Y, W, psi = rep(0, ng), backend = "cpu", verbose = FALSE),
               "strictly positive")
  expect_error(fitNB(Y, W, psi = -p0, backend = "cpu", verbose = FALSE),
               "strictly positive")
  expect_error(fitNB(Y, W, psi = c(NA, p0[-1]), backend = "cpu", verbose = FALSE),
               "strictly positive")
  expect_error(fitNB(Y, W, psi = c(Inf, p0[-1]), backend = "cpu", verbose = FALSE),
               "strictly positive")
})
