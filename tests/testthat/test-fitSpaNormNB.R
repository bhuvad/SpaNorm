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
