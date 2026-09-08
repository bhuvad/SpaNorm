# edgeR v4's quasi-likelihood dispersion, as generic NB-GLM machinery with a
# CPU and a torch backend. edgeR computes the adjusted deviance and effective
# df in C (src/ql_glm.c compute_adjust_vec, src/ql_weights.c compute_weight),
# which a Bioconductor package may not reach into. The definition is not the
# Chebyshev tables -- those are a fast path -- it is the phi >= 4.001 branch:
# per observation, the unit deviance's first two moments under the fitted NB
# are obtained by direct summation over the pmf, then matched to a scaled
# chi-square, giving
#     w0 = 2 E[d] / Var[d]        (rescales the unit deviance)
#     w1 = 2 E[d]^2 / Var[d]      (that observation's effective df)
# edgeR is the oracle for the whole pipeline.

# skip_if_no_torch()/skip_if_no_gpu() live in helper-gpu.R (auto-sourced by
# testthat). These are CPU float64 tensors, so they need a usable torch
# backend, not a GPU.
cpu_tensor <- function(m) {
  skip_if_no_torch()
  torch::torch_tensor(as.matrix(m), dtype = torch::torch_float64())
}

test_that("nbUnitDeviance matches edgeR's unit deviance, zeros and per-gene phi included", {
  set.seed(1)
  ng <- 6L; n <- 50L
  mu <- matrix(exp(rnorm(ng * n, log(0.5), 1.5)), ng, n)
  phi <- c(0.05, 0.2, 1, 5, 20, 5)
  y <- matrix(stats::rnbinom(ng * n, mu = mu, size = 1 / phi), ng, n)
  y[1, 1:10] <- 0
  got <- nbUnitDeviance(y, mu, phi)
  ref <- t(sapply(seq_len(ng), function(g) edgeR::nbinomUnitDeviance(y[g, ], mu[g, ], phi[g])))
  # edgeR's C unit deviance carries its own series approximations: 1e-6 is
  # the oracle's precision, not the formula's
  expect_equal(unname(got), unname(ref), tolerance = 1e-6)
  expect_equal(nbUnitDeviance(y[2, , drop = FALSE], mu[2, , drop = FALSE], 0.2),
               got[2, , drop = FALSE])
})

test_that("nbUnitDeviance on torch tensors equals the matrix result", {
  set.seed(2)
  ng <- 5L; n <- 40L
  mu <- matrix(exp(rnorm(ng * n, log(0.3), 1.2)), ng, n)
  phi <- c(0.1, 0.5, 2, 5, 10)
  y <- matrix(stats::rnbinom(ng * n, mu = mu, size = 1 / phi), ng, n)
  got <- nbUnitDeviance(cpu_tensor(y), cpu_tensor(mu), phi)
  expect_true(is_torch_tensor(got))
  expect_equal(toRMatrix(got), nbUnitDeviance(y, mu, phi), tolerance = 1e-10)
})

test_that("nbDevianceMoments matches a chi-square moment match by construction", {
  # if d ~ c * chisq_nu then E = c*nu and Var = 2c^2 nu, so w1 = nu and
  # w0 = nu / E. Large mu with small phi is the near-Poisson limit where the
  # unit deviance is close to chisq_1, i.e. w1 -> 1.
  m <- nbDevianceMoments(mu = 500, phi = 1e-6)
  expect_equal(m$w1, 1, tolerance = 0.02)
  expect_equal(m$w0, 1, tolerance = 0.02)
})

test_that("nbDevianceMoments halves the df at a sparse operating point", {
  # at mu ~ 0.1 and phi ~ 5 (a CosMx cohort's medians) each observation
  # carries about HALF a degree of freedom, not one, which is why the raw
  # deviance / (n - p) reads far below the Pearson dispersion there
  m <- nbDevianceMoments(mu = 0.1, phi = 5)
  expect_lt(m$w1, 0.6)
  expect_gt(m$w1, 0.4)
  expect_gt(m$w0, 1.5)          # and the unit deviance is scaled up, not down
})

test_that("nbDevianceMoments has converged at the default window", {
  # An sd-based window fails here: at mu = 1, phi = 20 the NB is heavy enough
  # that 12 sd got w1 wrong by 41%. The quantile window holds to ~1e-6 when the
  # tail left outside is tightened by three orders of magnitude.
  g <- expand.grid(mu = c(0.1, 1, 10, 100, 1000), phi = c(0.01, 0.2, 1, 5, 20))
  a <- nbDevianceMoments(g$mu, g$phi, eps = 1e-10)
  b <- nbDevianceMoments(g$mu, g$phi, eps = 1e-13)
  expect_equal(a$w0, b$w0, tolerance = 1e-4)
  expect_equal(a$w1, b$w1, tolerance = 1e-4)
})

test_that("qlDispersion reproduces edgeR::glmQLFit's adjusted deviance and df", {
  set.seed(11)
  ng <- 40L; n <- 8L
  grp <- factor(rep(c("a", "b"), each = n / 2))
  design <- stats::model.matrix(~ grp)
  mu0 <- rep(c(2, 20, 200, 2000), length.out = ng)
  y <- matrix(stats::rnbinom(ng * n, mu = mu0, size = 1 / 0.2), ng, n)
  disp <- 0.2
  fit <- edgeR::glmQLFit(y, design = design, dispersion = disp,
                         legacy = FALSE, abundance.trend = FALSE)
  got <- qlDispersion(y, mu = fit$fitted.values, phi = disp, design = design,
                      prior = fit$average.ql.dispersion,
                      leverage = "exact", moments = "cell")
  # 3e-3, and the slack is edgeR's, not ours: widening our summation window
  # moves our s2 by 5e-7 while the gap to edgeR stays pinned at 1e-3, so we
  # have converged and edgeR's Chebyshev fast path has not
  expect_equal(got$df, as.numeric(fit$df.residual.adj), tolerance = 3e-3)
  expect_equal(got$deviance, as.numeric(fit$deviance.adj), tolerance = 3e-3)
  expect_equal(got$s2, as.numeric(fit$deviance.adj / fit$df.residual.adj),
               tolerance = 3e-3)
})

test_that("qlDispersion's trace approximation tracks the exact leverages", {
  # spiDE's designs have n >> p, where every leverage is ~p/n and the exact
  # per-observation hat values collapse to the trace
  set.seed(12)
  ng <- 20L; n <- 400L
  design <- cbind(1, rnorm(n), rnorm(n))
  y <- matrix(stats::rnbinom(ng * n, mu = 30, size = 1 / 0.3), ng, n)
  mu <- matrix(30, ng, n)
  ex <- qlDispersion(y, mu, phi = 0.3, design = design, leverage = "exact", moments = "cell")
  ap <- qlDispersion(y, mu, phi = 0.3, design = design, leverage = "trace", moments = "cell")
  expect_equal(ap$df, ex$df, tolerance = 1e-3)
  expect_equal(ap$s2, ex$s2, tolerance = 1e-3)
})

test_that("qlDispersion's shared moment table matches per-cell moments across genes with different phi", {
  # The moments depend only on (mu, phi). One table over (log mu, log phi),
  # interpolated, serves every gene: the per-gene grid it replaces cost
  # ngrid pmf summations PER GENE, which was the dominant cost of the
  # pre-pass on small designs.
  set.seed(14)
  ng <- 24L; n <- 500L
  design <- cbind(1, rnorm(n))
  mu <- matrix(exp(rnorm(ng * n, log(0.1), 1.2)), ng, n)
  phi <- exp(seq(log(0.05), log(30), length.out = ng))
  y <- matrix(stats::rnbinom(ng * n, mu = mu, size = 1 / phi), ng, n)
  ex <- qlDispersion(y, mu, phi = phi, design = design, moments = "cell")
  tb <- qlDispersion(y, mu, phi = phi, design = design, moments = "table")
  expect_equal(tb$s2, ex$s2, tolerance = 2e-3)
  expect_equal(tb$df, ex$df, tolerance = 2e-3)
  # and a single-phi call reproduces the old per-gene grid path
  gr <- qlDispersion(y, mu, phi = 5, design = design, moments = "grid")
  tb1 <- qlDispersion(y, mu, phi = 5, design = design, moments = "table")
  expect_equal(tb1$s2, gr$s2, tolerance = 2e-3)
})

test_that("qlDispersion accepts the design's column count in place of the design", {
  set.seed(15)
  ng <- 8L; n <- 300L
  design <- cbind(1, rnorm(n))
  mu <- matrix(exp(rnorm(ng * n, log(1), 1)), ng, n)
  y <- matrix(stats::rnbinom(ng * n, mu = mu, size = 1 / 2), ng, n)
  a <- qlDispersion(y, mu, phi = 2, design = design)
  b <- qlDispersion(y, mu, phi = 2, p = 2L)
  expect_equal(a$s2, b$s2)
  expect_equal(a$df, b$df)
})

test_that("qlDispersion on torch tensors equals the matrix result", {
  set.seed(16)
  ng <- 10L; n <- 300L
  mu <- matrix(exp(rnorm(ng * n, log(0.2), 1.3)), ng, n)
  phi <- exp(seq(log(0.1), log(20), length.out = ng))
  y <- matrix(stats::rnbinom(ng * n, mu = mu, size = 1 / phi), ng, n)
  cpu <- qlDispersion(y, mu, phi = phi, p = 3L)
  dev <- qlDispersion(cpu_tensor(y), cpu_tensor(mu), phi = phi, p = 3L)
  expect_equal(dev$s2, cpu$s2, tolerance = 1e-8)
  expect_equal(dev$df, cpu$df, tolerance = 1e-8)
  expect_equal(dev$deviance, cpu$deviance, tolerance = 1e-8)
})

test_that("a caller-supplied moment table makes qlDispersion invariant to how genes are split", {
  # a blocked caller must get the same per-gene dispersion whether it scores
  # all genes at once or block by block: the moments table is therefore built
  # once over a range the caller fixes, not from each block's own range
  set.seed(17)
  ng <- 20L; n <- 300L
  mu <- matrix(exp(rnorm(ng * n, log(0.3), 1.4)), ng, n)
  phi <- exp(seq(log(0.05), log(20), length.out = ng))
  y <- matrix(stats::rnbinom(ng * n, mu = mu, size = 1 / phi), ng, n)
  tab <- qlMomentTable(lmu_range = log(c(1e-8, max(y) + 1)), lphi_range = log(range(phi)))
  whole <- qlDispersion(y, mu, phi, p = 3L, table = tab)
  parts <- lapply(list(1:7, 8:20), function(i) qlDispersion(y[i, ], mu[i, ], phi[i], p = 3L, table = tab))
  expect_identical(whole$s2, c(parts[[1]]$s2, parts[[2]]$s2))
  expect_identical(whole$df, c(parts[[1]]$df, parts[[2]]$df))
  # and it agrees with the per-cell moments to the interpolation tolerance
  ex <- qlDispersion(y, mu, phi, p = 3L, moments = "cell")
  expect_equal(whole$s2, ex$s2, tolerance = 2e-3)
  # the same table on tensors
  skip_if_no_torch()
  dev <- qlDispersion(cpu_tensor(y), cpu_tensor(mu), phi, p = 3L, table = tab)
  expect_equal(dev$s2, whole$s2, tolerance = 1e-8)
})
