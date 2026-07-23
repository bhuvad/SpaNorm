test_that("fitSpaNormNB internal function checks works", {
  set.seed(36)
  emat = matrix(rpois(6 * 4e3, 10), 6, 4e3)
  coords = cbind(1:4e3, 1:4e3)

  # maxit
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = 0), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = -1), "greater")
  expect_error(fitSpaNormNB(emat, c(), c(), maxit.psi = 1.1), "integer")
})

test_that("fitNBGivenPsi internal function checks works", {
  Y = matrix(0, 6, 10)
  gmean = psi = rep(0, 6)
  W = matrix(0, 10, 3)
  alpha = matrix(0, 6, 3)

  # maxit
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = 0), "greater")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = -1), "greater")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, maxit.nb = 1.1), "integer")

  # step.factor
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 0), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 1), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = 1.1), "0,1")
  expect_error(fitNBGivenPsi(Y, W, gmean, alpha, psi, lambda.a = 0.1, step.factor = -1.1), "0,1")
})

test_that("fitNB recovers a generic per-gene NB GLM (no intercept term, gmean = 0)", {
  set.seed(10)
  ng <- 40; ns <- 250
  W <- cbind(intercept = 1, a = as.numeric(scale(rnorm(ns))), b = as.numeric(scale(rnorm(ns))))
  true.alpha <- cbind(rnorm(ng, 1, 0.4), rnorm(ng, 0, 0.7), rnorm(ng, 0, 0.7))
  mu <- exp(tcrossprod(true.alpha, W)) # ng x ns, no separate intercept -> encoded in W[,1]
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(mu), size = 5), ng, ns,
              dimnames = list(paste0("g", 1:ng), paste0("c", 1:ns)))

  fit <- fitNB(Y, W, backend = "cpu", verbose = FALSE)

  # generic fit carries no per-gene intercept (log(mu) = W %*% t(alpha))
  expect_true(all(fit$gmean == 0))
  expect_equal(dim(fit$alpha), c(ng, ncol(W)))
  expect_true(all(is.finite(fit$alpha)) && all(is.finite(fit$psi)))
  # coefficients recover the generative truth
  expect_gt(cor(as.vector(fit$alpha), as.vector(true.alpha)), 0.9)
})

test_that("fitNB ridge (lambda.a) shrinks coefficients and winsor clips outliers", {
  set.seed(11)
  ng <- 40; ns <- 120
  W <- cbind(intercept = 1, cov = as.numeric(scale(seq_len(ns))))
  Y <- matrix(rpois(ng * ns, 6), ng, ns, dimnames = list(paste0("g", 1:ng), paste0("c", 1:ns)))
  # two genes with strong, opposite covariate slopes -> outlier alpha[, 2]
  Y[1, ] <- rpois(ns, exp(2 + 3 * W[, "cov"]))
  Y[2, ] <- rpois(ns, exp(2 - 3 * W[, "cov"]))

  # ridge shrinks the covariate coefficient toward 0
  f0 <- fitNB(Y, W, lambda.a = 0,   winsor = Inf, backend = "cpu", verbose = FALSE)
  fr <- fitNB(Y, W, lambda.a = 1e3, winsor = Inf, backend = "cpu", verbose = FALSE)
  expect_lt(mean(abs(fr$alpha[, 2])), mean(abs(f0$alpha[, 2])))

  # winsor clips the outlier coefficients relative to no winsorisation
  fw <- fitNB(Y, W, lambda.a = 0, winsor = 2, backend = "cpu", verbose = FALSE)
  expect_lt(max(abs(fw$alpha[, 2])), max(abs(f0$alpha[, 2])))
})

test_that("fitNB forwards extra fitting parameters (maxn.psi, step.factor) via ...", {
  set.seed(12)
  ng <- 10; ns <- 40
  W <- cbind(intercept = 1, cov = as.numeric(scale(seq_len(ns))))
  Y <- matrix(rpois(ng * ns, 6), ng, ns)

  expect_no_error(fitNB(Y, W, maxn.psi = 20, backend = "cpu", verbose = FALSE))
  expect_no_error(fitNB(Y, W, step.factor = 0.3, backend = "cpu", verbose = FALSE))
})

# ---- Memory-bounded (gene-blocked) GPU fitting -----------------------------
# skip_if_no_gpu()/gpu_tol()/tinyGpuBudget() live in helper-gpu.R (auto-sourced
# by testthat before any test file runs)

test_that("fitNB/fitSpaNormNB forward gpu.mem.budget via ...", {
  set.seed(13)
  ng <- 10; ns <- 40
  W <- cbind(intercept = 1, cov = as.numeric(scale(seq_len(ns))))
  Y <- matrix(rpois(ng * ns, 6), ng, ns)

  expect_no_error(fitNB(Y, W, gpu.mem.budget = 1e9, backend = "cpu", verbose = FALSE))
  # backend = "cpu" never engages blocking, so even a pathologically tiny
  # budget must not error or change anything
  expect_no_error(fitNB(Y, W, gpu.mem.budget = 1, backend = "cpu", verbose = FALSE))
})

test_that("gpu.mem.budget is validated", {
  # backend = "cpu" never evaluates the budget at all (geneBlockCount()
  # short-circuits on the backend check before its lazily-evaluated
  # budget.bytes argument is forced) -- by design, an irrelevant value for a
  # path that never touches the accelerator is silently ignored, not an
  # error. getGPUMemoryBudget()'s own validation is covered directly in
  # test-gpuFunctions.R; here just confirm it actually fires on the path
  # that does consume the value.
  skip_if_no_gpu()
  set.seed(14)
  ng <- 5; ns <- 20
  W <- cbind(intercept = 1, cov = as.numeric(scale(seq_len(ns))))
  Y <- matrix(rpois(ng * ns, 6), ng, ns)

  expect_error(fitNB(Y, W, gpu.mem.budget = -1, backend = "gpu", verbose = FALSE), "positive")
  expect_error(fitNB(Y, W, gpu.mem.budget = 0, backend = "gpu", verbose = FALSE), "positive")
})

test_that("geneBlockCount is 1L (a single whole-gene block) on typical small test fixtures", {
  # sanity check: the existing GPU test suite's small matrices stay within a
  # single block (the fitting machinery always iterates gene-blocks; this
  # confirms the common/small case degenerates to exactly one)
  expect_equal(geneBlockCount(40, 250, 3, backend = "cpu"), 1L)
  skip_if_no_gpu()
  expect_equal(geneBlockCount(40, 250, 3, backend = "gpu", budget.bytes = getGPUMemoryBudget()), 1L)
})

test_that("a forced tiny gpu.mem.budget produces a fit matching the unblocked GPU/CPU fit", {
  skip_if_no_gpu()
  resetGPUCache()
  on.exit(resetGPUCache())

  set.seed(15)
  ng <- 37; ns <- 60 # deliberately not a multiple of the forced block count
  true.alpha <- cbind(rnorm(ng, 1, 0.4), rnorm(ng, 0, 0.7), rnorm(ng, 0, 0.7))
  W <- cbind(intercept = 1, a = as.numeric(scale(rnorm(ns))), b = as.numeric(scale(rnorm(ns))))
  mu <- exp(tcrossprod(true.alpha, W))
  Y <- matrix(rnbinom(ng * ns, mu = as.vector(mu), size = 5), ng, ns)

  # reset the seed before each fit so the (identical) dispersion subsample
  # draw doesn't confound the blocked-vs-unblocked comparison
  set.seed(20); fit.cpu <- fitNB(Y, W, backend = "cpu", verbose = FALSE)
  set.seed(20); fit.gpu.unblocked <- fitNB(Y, W, backend = "gpu", gpu.mem.budget = Inf, verbose = FALSE)

  # force many small blocks: a budget far below what a single block needs
  tiny.budget <- tinyGpuBudget(ng, ns)
  nb <- geneBlockCount(ng, ns, ncol(W), backend = "gpu", budget.bytes = tiny.budget)
  expect_gt(nb, 1L) # sanity: the forced budget actually engages blocking

  set.seed(20); fit.gpu.blocked <- fitNB(Y, W, backend = "gpu", gpu.mem.budget = tiny.budget, verbose = FALSE)

  expect_equal(fit.gpu.blocked$alpha, fit.gpu.unblocked$alpha, tolerance = gpu_tol())
  expect_equal(fit.gpu.blocked$psi, fit.gpu.unblocked$psi, tolerance = gpu_tol())
  expect_equal(fit.gpu.blocked$loglik, fit.gpu.unblocked$loglik, tolerance = gpu_tol())

  # and matches the CPU reference within the same looser cross-device tolerance
  # already used elsewhere in this codebase for GPU-vs-CPU comparisons
  expect_equal(fit.gpu.blocked$alpha, fit.cpu$alpha, tolerance = 1e-3)
})

test_that("a forced tiny gpu.mem.budget matches unblocked for the is.spanorm (shared logLS) model too", {
  skip_if_no_gpu()
  resetGPUCache()
  on.exit(resetGPUCache())

  set.seed(16)
  ng <- 29; ns <- 50
  emat <- matrix(rpois(ng * ns, 10), ng, ns)
  coords <- cbind(runif(ns, 1, 100), runif(ns, 1, 100))
  ls <- colSums(emat) # precomputed library size (bypasses scran::quickCluster,
  # which needs more cells than this small synthetic fixture has)

  quiet <- function(...) NULL
  set.seed(21)
  fit.unblocked <- fitSpaNorm(emat, coords, sample.p = 1, gene.model = "nb", df.tps = 2,
                               batch = NULL, LS = ls, backend = "gpu", gpu.mem.budget = Inf,
                               tol = 1e-2, msgfun = quiet)

  tiny.budget <- tinyGpuBudget(ng, ns)
  nb <- geneBlockCount(ng, ns, 5, backend = "gpu", budget.bytes = tiny.budget)
  expect_gt(nb, 1L)

  set.seed(21)
  fit.blocked <- fitSpaNorm(emat, coords, sample.p = 1, gene.model = "nb", df.tps = 2,
                             batch = NULL, LS = ls, backend = "gpu", gpu.mem.budget = tiny.budget,
                             tol = 1e-2, msgfun = quiet)

  # the shared logLS coefficient (column 1) must be identical across all
  # genes in both fits, and the two fits must agree with each other. This
  # compares a *full* multi-outer-iteration SpaNorm fit (dispersion
  # re-estimation loop wrapping the is.spanorm inner IRLS loop, with its
  # extra Wa1-fold/gmean-fold small-stage arithmetic on top of the blocked
  # passes) rather than a single operation, so float32-on-MPS rounding
  # compounds across iterations well beyond gpu_tol() -- a looser, still
  # tight, tolerance is used here for that reason (see the "reduce-then-
  # broadcast... small floating-point summation-order differences" trade-off
  # noted in the design). all.equal()'s mean-relative-difference check is
  # used rather than expect_equal()'s element-wise waldo comparison, since a
  # single near-zero alpha entry (common for this coefficient) inflates
  # element-wise relative error without reflecting the fit's actual quality.
  fit.tol <- if (getBackendDevice() == "mps") 1e-2 else 1e-6
  expect_equal(length(unique(fit.blocked$alpha[, 1])), 1L)
  expect_true(isTRUE(all.equal(fit.blocked$alpha, fit.unblocked$alpha, tolerance = fit.tol)))
  expect_true(isTRUE(all.equal(fit.blocked$gmean, fit.unblocked$gmean, tolerance = fit.tol)))
  expect_true(isTRUE(all.equal(fit.blocked$psi, fit.unblocked$psi, tolerance = fit.tol)))
})
