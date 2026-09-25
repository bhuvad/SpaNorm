# ls = "joint": the pooled, profiled Newton step on SpaNorm's shared
# library-size coefficient (R/polishSpaNormJoint.R; spec section 4).

# The total penalised log-likelihood of a SpaNorm-model fit over the genes
# `genes`: unwinsorised mean, each gene's own dispersion and penalty.
.joint_total <- function(fit, Y, genes = seq_len(nrow(Y))) {
  p <- .spaNormPolishProblem(fit)
  A <- cbind(fit$gmean, fit$alpha[, -1])
  sum(vapply(genes, function(g) {
    mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
    sum(dnbinom(Y[g, ], size = 1 / fit$psi[g], mu = mu, log = TRUE)) -
      0.5 * sum(p$pen * A[g, ]^2)
  }, 0))
}

# The pooled score for a1 at a fit, over the genes `genes`.
.joint_a1_score <- function(fit, Y, genes = seq_len(nrow(Y))) {
  p <- .spaNormPolishProblem(fit)
  A <- cbind(fit$gmean, fit$alpha[, -1])
  sum(vapply(genes, function(g) {
    mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
    sum(p$w1 * (Y[g, ] - mu) / (1 + fit$psi[g] * mu))
  }, 0))
}

test_that("joint a1 matches a full optim() of the joint objective", {
  set.seed(30)
  n <- 120; G <- 5
  w1 <- log(runif(n, 0.5, 2)); x2 <- rnorm(n)
  X <- cbind(1, x2); pen <- c(0, 2)
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.3 + g / 5 + 0.8 * w1 + 0.3 * x2), size = 5), numeric(n)))
  psi <- rep(0.2, G)
  negobj <- function(par) {
    a1 <- par[1]; A <- matrix(par[-1], G, 2)
    -sum(vapply(seq_len(G), function(g) {
      mu <- exp(as.numeric(X %*% A[g, ] + a1 * w1))
      sum(dnbinom(Y[g, ], size = 1 / psi[g], mu = mu, log = TRUE)) - 0.5 * sum(pen * A[g, ]^2)
    }, 0))
  }
  opt <- optim(c(1, rep(0, 2 * G)), negobj, method = "BFGS",
               control = list(maxit = 5000, reltol = 1e-14))
  # the oracle itself converged (ruling 6): BFGS reports success and the
  # analytic gradient of the joint objective vanishes at its answer
  expect_identical(opt$convergence, 0L)
  grad <- function(par) {
    a1 <- par[1]; A <- matrix(par[-1], G, 2)
    ga <- 0; gA <- matrix(0, G, 2)
    for (g in seq_len(G)) {
      mu <- exp(as.numeric(X %*% A[g, ] + a1 * w1))
      r <- (Y[g, ] - mu) / (1 + psi[g] * mu)
      ga <- ga + sum(w1 * r)
      gA[g, ] <- as.numeric(crossprod(X, r)) - pen * A[g, ]
    }
    c(ga, gA)
  }
  expect_lt(max(abs(grad(opt$par))), 1e-3)

  prob <- list(X = X, pen = pen, a1 = 1, offset = 1 * w1, w1 = w1,
               A0 = matrix(0, G, 2), psi = psi, cells_idx = rep(TRUE, n))
  pol <- polishNB(Y, X, prob$A0, psi, lambda.a = pen, offset = prob$offset,
                  psi.method = "fixed")
  j <- .polishSharedLS(Y, prob, pol, psi.method = "fixed")
  expect_equal(j$a1, opt$par[1], tolerance = 1e-4)
  expect_equal(j$pol$alpha, matrix(opt$par[-1], G, 2), tolerance = 1e-3)
  expect_lt(abs(j$score), 1e-6 * sum(Y))
  expect_true(j$converged)
  expect_gt(j$iterations, 0L)
  expect_identical(j$singular, 0L)
  # the profiled SE of a1 is the one the Hessian of the profile likelihood gives
  expect_true(is.finite(j$se) && j$se > 0)
})

test_that("joint polish never lowers the total penalised likelihood", {
  spe <- .polish_spe()
  fx <- S4Vectors::metadata(polishSpaNorm(spe, ls = "fixed", verbose = FALSE))$SpaNorm
  jt <- S4Vectors::metadata(polishSpaNorm(spe, ls = "joint", verbose = FALSE))$SpaNorm
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  tot <- function(fit) {
    p <- .spaNormPolishProblem(fit)
    A <- cbind(fit$gmean, fit$alpha[, -1])
    sum(vapply(seq_len(nrow(Y)), function(g) {
      mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
      sum(dnbinom(Y[g, ], size = 1 / fit$psi[g], mu = mu, log = TRUE)) -
        0.5 * sum(p$pen * A[g, ]^2)
    }, 0))
  }
  expect_gte(tot(jt), tot(fx) - 1e-8 * abs(tot(fx)))
  expect_true(length(unique(jt$alpha[, 1])) == 1)        # still shared
})

test_that("the joint fit is stationary in a1 and in every gene's coefficients", {
  # spec section 7, item 4; tol = 1e-12 for the per-gene zero-score check, as
  # in test-polishSpaNorm.R
  spe <- .polish_spe()
  out <- polishSpaNorm(spe, ls = "joint", tol = 1e-12, verbose = FALSE)
  f <- S4Vectors::metadata(out)$SpaNorm
  f0 <- S4Vectors::metadata(out)$SpaNormUnpolished
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  s <- .polishSlot(f)$settings
  expect_identical(s$ls, "joint")
  expect_identical(s$a1.input, f0$alpha[1, 1])
  expect_identical(s$a1, f$alpha[1, 1])
  expect_false(isTRUE(all.equal(s$a1, s$a1.input)))        # a1 moved
  # the pooled score, recomputed here from the stored fit, is at zero
  U <- .joint_a1_score(f, Y)
  expect_lt(abs(U), 1e-6 * sum(Y))
  expect_equal(s$ls.score, U, tolerance = 1e-6 * sum(Y))
  expect_lt(abs(s$ls.score) * s$ls.se, 1e-6)                # the stopping rule
  expect_true(s$ls.converged)
  expect_true(s$ls.iterations >= 1L && s$ls.iterations <= 10L)
  expect_identical(s$ls.singular, 0L)
  expect_true(is.finite(s$ls.se) && s$ls.se > 0)
  # every gene is at its own optimum at the joint a1
  expect_lt(max(.polish_scaled_score(f, Y)), 1e-6)
  expect_equal(f$psi, f0$psi)                                 # psi.method = "fixed"
  # the per-gene log-likelihood is the one at the returned fit
  genes <- .polishSlot(f)$genes
  expect_equal(sum(genes$loglik), .joint_total(f, Y), tolerance = 1e-10)
  expect_true(all(genes$polished))
})

test_that("ls.se is the profiled standard error of a1", {
  # the profiled information is the curvature of the profile log-likelihood
  # l(a1) = max_A l(a1, A): check 1/ls.se^2 against a second difference of
  # the profile, each point a converged polish at that a1
  spe <- .polish_spe()
  f <- S4Vectors::metadata(polishSpaNorm(spe, ls = "joint", tol = 1e-12,
                                         verbose = FALSE))$SpaNorm
  s <- .polishSlot(f)$settings
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  p <- .spaNormPolishProblem(f)
  A <- cbind(f$gmean, f$alpha[, -1])
  prof <- function(a1) {
    pl <- polishNB(Y, p$X, A, f$psi, lambda.a = p$pen, offset = a1 * p$w1,
                   psi.method = "fixed", tol = 1e-14)
    sum(pl$loglik)
  }
  h <- 5 * s$ls.se
  curv <- -(prof(s$a1 + h) - 2 * prof(s$a1) + prof(s$a1 - h)) / h^2
  # Fisher (expected) against observed information: close, not equal
  expect_equal(1 / s$ls.se^2, curv, tolerance = 0.05)
})

test_that("joint: an all-zero gene is held out, keeps its fit, and shares the new a1", {
  # rulings 2 and 8: the zero gene is not polished and does not inform a1; its
  # gmean, alpha[, -1] and psi are its input values bit for bit, and its
  # alpha[, 1] is the new shared a1, so the column stays shared
  spe <- .polish_spe()
  cnt <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  cnt[1, ] <- 0
  SummarizedExperiment::assay(spe, "counts") <- cnt
  out <- polishSpaNorm(spe, ls = "joint", verbose = FALSE)
  f <- S4Vectors::metadata(out)$SpaNorm
  f0 <- S4Vectors::metadata(out)$SpaNormUnpolished
  s <- .polishSlot(f)$settings
  pol <- .polishSlot(f)$genes
  expect_false(pol$polished[1])
  expect_true(all(pol$polished[-1]))
  expect_identical(pol$iterations[1], 0L)
  expect_true(is.na(pol$loglik[1]))
  expect_identical(f$gmean[1], f0$gmean[1])
  expect_identical(f$alpha[1, -1], f0$alpha[1, -1])
  expect_identical(f$psi[1], f0$psi[1])
  expect_true(length(unique(f$alpha[, 1])) == 1)            # still shared
  expect_identical(f$alpha[1, 1], s$a1)
  expect_false(isTRUE(all.equal(s$a1, s$a1.input)))
  # the sums ran over the polished genes only: the score is zero over them
  expect_lt(abs(.joint_a1_score(f, cnt, genes = 2:nrow(cnt))), 1e-6 * sum(cnt))
  expect_equal(s$ls.score, .joint_a1_score(f, cnt, genes = 2:nrow(cnt)),
               tolerance = 1e-6 * sum(cnt))
})

test_that("joint with no polished gene warns and leaves a1 where it was", {
  # ruling 9: an all-zero counts matrix holds every gene out
  f0 <- S4Vectors::metadata(.polish_spe())$SpaNorm
  Y0 <- matrix(0, f0$ngenes, f0$ncells)
  expect_warning(f <- .polishSpaNormFit(f0, Y0, ls = "joint"), "nothing informs")
  s <- .polishSlot(f)$settings
  expect_identical(s$ls, "joint")
  expect_identical(s$ls.iterations, 0L)
  expect_identical(s$a1, s$a1.input)
  expect_identical(s$a1, f0$alpha[1, 1])
  expect_true(all(c("ls.score", "ls.se", "ls.singular", "ls.converged") %in% names(s)))
  expect_identical(s$ls.singular, 0L)
  expect_false(s$ls.converged)
  expect_identical(f$alpha, f0$alpha)
  expect_identical(f$gmean, f0$gmean)
  expect_identical(f$psi, f0$psi)
  expect_false(any(.polishSlot(f)$genes$polished))

  # and the same when the genes are passed but the cold polish converged none
  set.seed(4)
  n <- 50
  X <- cbind(1, rnorm(n)); w1 <- rnorm(n, sd = 0.3)
  Y <- matrix(rpois(3 * n, 5), 3, n)
  prob <- list(X = X, pen = c(0, 1), a1 = 0.7, offset = 0.7 * w1, w1 = w1,
               A0 = matrix(0, 3, 2), psi = rep(0.1, 3), cells_idx = rep(TRUE, n))
  pol <- list(alpha = prob$A0, psi = prob$psi, loglik = rep(NA_real_, 3),
              polish = data.frame(iterations = 0L, psi_fitnb = prob$psi,
                                  restarted = FALSE, capped = FALSE,
                                  singular = TRUE, psi_bound = FALSE,
                                  polished = FALSE))
  expect_warning(j <- .polishSharedLS(Y, prob, pol), "nothing informs")
  expect_identical(j$a1, 0.7)
  expect_identical(j$iterations, 0L)
  expect_identical(j$pol, pol)
})

test_that("a gene with singular information is left out of both U and I", {
  # ruling 4. Gene 1 is Poisson, absent from half the cells and very bright in
  # the other half, with an unpenalised indicator for that half: its
  # information has weights exp(-30) on one half and ~1e6 on the other, so
  # X' D X is numerically rank one and solve() refuses it
  set.seed(8)
  n <- 80
  ind <- rep(0:1, each = n / 2)
  X <- cbind(1, ind); pen <- c(0, 0)
  w1 <- rnorm(n, sd = 0.3)
  A <- rbind(c(-40, 54), c(1, 0.2), c(0.5, -0.1))
  psi <- c(0, 0.1, 0.2)
  Y <- rbind(rpois(n, pmax(exp(A[1, 1] + A[1, 2] * ind + 0.9 * w1), nbMuFloor())),
             rnbinom(n, mu = exp(A[2, 1] + A[2, 2] * ind + 0.9 * w1), size = 10),
             rnbinom(n, mu = exp(A[3, 1] + A[3, 2] * ind + 0.9 * w1), size = 5))
  solver <- nbNewtonSolver(X, pen)
  mu1 <- pmax(exp(as.numeric(X %*% A[1, ] + 0.9 * w1)), nbMuFloor())
  expect_null(solver$solve(solver$factor(mu1), c(1, 1)))      # the premise

  si <- .lsScoreInfo(Y, A, psi, 0.9, X, w1, pen, solver)
  ref <- .lsScoreInfo(Y[2:3, ], A[2:3, ], psi[2:3], 0.9, X, w1, pen, solver)
  expect_identical(si$excluded, 1L)
  expect_identical(si$n, 2L)
  expect_identical(ref$n, 2L)
  expect_equal(si$U, ref$U, tolerance = 1e-12)   # gene 1's score is NOT in U
  expect_equal(si$I, ref$I, tolerance = 1e-12)   # nor its unprofiled w1'Dw1 in I
  # blocking does not change the sums
  sb <- .lsScoreInfo(Y, A, psi, 0.9, X, w1, pen, solver, block = 1L)
  expect_equal(sb[c("U", "I", "n")], si[c("U", "I", "n")], tolerance = 1e-12)
  expect_identical(sb$excluded, si$excluded)

  # through the joint step it is counted in ls.singular
  prob <- list(X = X, pen = pen, a1 = 0.9, offset = 0.9 * w1, w1 = w1,
               A0 = A, psi = psi, cells_idx = rep(TRUE, n))
  pol <- list(alpha = A, psi = psi, loglik = rep(0, 3),
              polish = data.frame(iterations = 1L, psi_fitnb = psi,
                                  restarted = FALSE, capped = FALSE,
                                  singular = FALSE, psi_bound = FALSE,
                                  polished = TRUE))
  # gene 1 stays in the objective (ruling 5), and its warm re-polish is
  # singular too, so it cannot follow a1: here that stalls the line search,
  # which is reported rather than hidden
  expect_warning(j <- .polishSharedLS(Y, prob, pol, psi.method = "fixed"),
                 "no step along the Newton direction")
  expect_identical(j$singular, 1L)
  expect_false(j$converged)
  expect_identical(j$pol$polish$singular, c(TRUE, FALSE, FALSE))
  expect_identical(j$pol$alpha[1, ], A[1, ])                  # held where it was
})

test_that("the joint step reports a cap it hits", {
  set.seed(30)
  n <- 120; G <- 5
  w1 <- log(runif(n, 0.5, 2)); x2 <- rnorm(n)
  X <- cbind(1, x2); pen <- c(0, 2)
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(0.3 + g / 5 + 0.8 * w1 + 0.3 * x2), size = 5), numeric(n)))
  prob <- list(X = X, pen = pen, a1 = 1, offset = w1, w1 = w1,
               A0 = matrix(0, G, 2), psi = rep(0.2, G), cells_idx = rep(TRUE, n))
  pol <- polishNB(Y, X, prob$A0, prob$psi, lambda.a = pen, offset = prob$offset,
                  psi.method = "fixed")
  expect_warning(j <- .polishSharedLS(Y, prob, pol, maxit.ls = 1L), "maxit.ls = 1 reached")
  expect_identical(j$iterations, 1L)
  expect_false(j$converged)
  # one step is still an ascent from the fixed-a1 polish
  expect_gt(sum(j$pol$loglik), sum(pol$loglik))
  # the per-gene diagnostics carry the warm pass's Newton iterations
  expect_true(all(j$pol$polish$iterations > pol$polish$iterations))
})

test_that("the profiled information matches a direct Schur complement", {
  set.seed(12)
  n <- 60; G <- 7
  X <- cbind(1, rnorm(n), rnorm(n)); pen <- c(0, 0.5, 3)
  w1 <- rnorm(n, sd = 0.4)
  A <- cbind(rnorm(G, 1), rnorm(G, 0, 0.2), rnorm(G, 0, 0.2))
  psi <- runif(G, 0.05, 0.5)
  Y <- t(vapply(seq_len(G), function(g)
    rnbinom(n, mu = exp(as.numeric(X %*% A[g, ]) + 0.8 * w1), size = 1 / psi[g]), numeric(n)))
  si <- .lsScoreInfo(Y, A, psi, 0.8, X, w1, pen, nbNewtonSolver(X, pen), block = 3L)
  U <- 0; I <- 0
  for (g in seq_len(G)) {
    mu <- exp(as.numeric(X %*% A[g, ]) + 0.8 * w1)
    d <- mu / (1 + psi[g] * mu)
    U <- U + sum(w1 * (Y[g, ] - mu) / (1 + psi[g] * mu))
    Xa <- cbind(X, w1)
    M <- crossprod(Xa * sqrt(d)) + diag(c(pen, 0))
    I <- I + M[4, 4] - M[4, 1:3] %*% solve(M[1:3, 1:3], M[1:3, 4])
  }
  expect_equal(si$U, U, tolerance = 1e-10)
  expect_equal(si$I, as.numeric(I), tolerance = 1e-10)
  expect_identical(si$n, 7L)
})

test_that("joint under psi.method = 'profile' is stationary at the reported psi", {
  spe <- .polish_spe()
  out <- polishSpaNorm(spe, ls = "joint", psi.method = "profile", tol = 1e-12,
                       verbose = FALSE)
  f <- S4Vectors::metadata(out)$SpaNorm
  f0 <- S4Vectors::metadata(out)$SpaNormUnpolished
  s <- .polishSlot(f)$settings
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  expect_true(s$ls.converged)
  expect_lt(abs(.joint_a1_score(f, Y)), 1e-6 * sum(Y))
  expect_lt(max(.polish_scaled_score(f, Y)), 1e-6)
  expect_false(isTRUE(all.equal(f$psi, f0$psi)))            # psi re-estimated
  # and it does not lose to the fixed-a1 profile polish
  fx <- S4Vectors::metadata(polishSpaNorm(spe, psi.method = "profile", tol = 1e-12,
                                          verbose = FALSE))$SpaNorm
  expect_gte(.joint_total(f, Y), .joint_total(fx, Y) - 1e-8 * abs(.joint_total(fx, Y)))
})

test_that("the null is polished jointly too, at its own a1", {
  # spec section 5: the null is polished with the full fit's settings, so a
  # joint full fit gets a joint null, each at its own optimum
  spe <- suppressWarnings(SpaNormSVG(.polish_spe(), backend = "cpu", verbose = FALSE))
  expect_warning(out <- polishSpaNorm(spe, ls = "joint", verbose = FALSE), "SVG")
  md <- S4Vectors::metadata(out)
  sf <- .polishSlot(md$SpaNorm)$settings
  sn <- .polishSlot(md$SpaNormNull)$settings
  expect_identical(sn$ls, "joint")
  expect_false(isTRUE(all.equal(sn$a1, sn$a1.input)))
  expect_false(isTRUE(all.equal(sn$a1, sf$a1)))              # its own a1
  expect_true(length(unique(md$SpaNormNull$alpha[, 1])) == 1)
  # and SpaNormSVG() accepts the pair as polished alike
  res <- SpaNormSVG(out, backend = "cpu", verbose = FALSE)
  expect_identical(S4Vectors::metadata(res)$SpaNormNull, md$SpaNormNull)
  expect_true(all(c("svg.F", "svg.p", "svg.fdr") %in%
                    colnames(SummarizedExperiment::rowData(res))))
})
