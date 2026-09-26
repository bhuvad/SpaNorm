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
  # tol.ls explicit and tight: this test compares the joint answer against an
  # independent optim() optimum, which needs a tighter stop than the package
  # default (1e-3) to keep its meaning
  j <- .polishSharedLS(Y, prob, pol, psi.method = "fixed", tol.ls = 1e-8)
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
  # the pooled score, recomputed here from the stored fit, is at zero in
  # units of a1's own SE (the stopping rule, checked independently of the
  # recorded score), and the recorded score is that recomputed one. The
  # former bounds, |U| < 1e-3 * sum(Y) and a tolerance of 1e-3 * sum(Y),
  # allowed a1 2.85 SEs from its optimum and any score (final review M-3).
  U <- .joint_a1_score(f, Y)
  expect_lt(abs(U) * s$ls.se, 1e-3)
  expect_lt(abs(s$ls.score - U), 1e-6 * abs(U) + 1e-8)
  expect_lt(abs(s$ls.score) * s$ls.se, 1e-3)                # the stopping rule
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
  # the sums ran over the polished genes only: the score is zero over them,
  # and it is the recorded one (final review M-3)
  U <- .joint_a1_score(f, cnt, genes = 2:nrow(cnt))
  expect_lt(abs(U) * s$ls.se, 1e-3)
  expect_lt(abs(s$ls.score - U), 1e-6 * abs(U) + 1e-8)
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
  expect_true(all(c("ls.score", "ls.se", "ls.singular", "ls.converged",
                    "ls.maxit", "ls.tol") %in% names(s)))
  expect_identical(s$ls.maxit, 10L)
  expect_identical(s$ls.tol, 1e-3)
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
  # fix round 1, amended ruling: the line search compares the total over the
  # SAME genes that formed U and I, so gene 1 (re-polished, but singular and
  # unable to follow a1) no longer stalls the step: the loop converges
  expect_no_warning(j <- .polishSharedLS(Y, prob, pol, psi.method = "fixed"))
  expect_true(j$converged)
  expect_identical(j$singular, 1L)
  expect_identical(j$pol$polish$singular, c(TRUE, FALSE, FALSE))
  expect_identical(j$pol$alpha[1, ], A[1, ])                  # held where it was
  # and the step is the one genes 2-3 alone take: gene 1 informs nothing
  prob2 <- prob
  prob2$A0 <- A[2:3, ]
  prob2$psi <- psi[2:3]
  pol2 <- list(alpha = A[2:3, ], psi = psi[2:3], loglik = rep(0, 2),
               polish = pol$polish[2:3, ])
  j2 <- .polishSharedLS(Y[2:3, ], prob2, pol2, psi.method = "fixed")
  expect_equal(j$a1, j2$a1, tolerance = 1e-10)
  expect_equal(j$pol$alpha[2:3, ], j2$pol$alpha, tolerance = 1e-10)
  expect_identical(j$iterations, j2$iterations)
  expect_equal(j$score, j2$score, tolerance = 1e-8)
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
  U <- .joint_a1_score(f, Y)
  expect_lt(abs(U) * s$ls.se, 1e-3)
  expect_lt(abs(s$ls.score - U), 1e-6 * abs(U) + 1e-8)
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

# ---- fix round 1 --------------------------------------------------------------

test_that("the joint loop runs at its own documented stop, not the per-gene one", {
  # polishSpaNorm()'s maxit/tol used to partial-match .polishSharedLS()'s
  # maxit.ls/tol.ls, so a per-gene tol of 1e-12 became the loop's tol and 50
  # its cap. The loop's stop is |U|/sqrt(I) < 1e-3 within 10 steps.
  out <- polishSpaNorm(.polish_spe(), ls = "joint", psi.method = "profile",
                       tol = 1e-12, verbose = FALSE)
  s <- .polishSlot(S4Vectors::metadata(out)$SpaNorm)$settings
  expect_identical(s$ls.tol, 1e-3)
  expect_identical(s$ls.maxit, 10L)
  expect_lte(s$ls.iterations, 10L)
  expect_identical(s$tol, 1e-12)                              # the per-gene tol
  expect_true(s$ls.converged)
  expect_lt(abs(s$ls.score) * s$ls.se, 1e-3)
})

test_that("every warm re-polish in the joint loop gets the caller's maxit and tol", {
  real <- polishNB
  seen <- list()
  local_mocked_bindings(polishNB = function(...) {
    args <- list(...)
    seen[[length(seen) + 1L]] <<- list(warm = isTRUE(args$warm),
                                       maxit = args$maxit, tol = args$tol)
    real(...)
  })
  out <- polishSpaNorm(.polish_spe(), ls = "joint", psi.method = "profile",
                       maxit = 30L, tol = 1e-5, verbose = FALSE)
  warm <- Filter(function(x) x$warm, seen)
  expect_gt(length(warm), 0L)
  expect_true(all(vapply(warm, function(x) identical(x$maxit, 30L), logical(1))))
  expect_true(all(vapply(warm, function(x) identical(x$tol, 1e-5), logical(1))))
  # the cold pass too, and the loop kept its own defaults
  cold <- Filter(function(x) !x$warm, seen)
  expect_true(all(vapply(cold, function(x) identical(x$maxit, 30L) &&
                           identical(x$tol, 1e-5), logical(1))))
  s <- .polishSlot(S4Vectors::metadata(out)$SpaNorm)$settings
  expect_identical(s[c("maxit", "tol", "ls.maxit", "ls.tol")],
                   list(maxit = 30L, tol = 1e-5, ls.maxit = 10L, ls.tol = 1e-3))
})

test_that("a non-positive profiled information warns and keeps the fixed-a1 polish", {
  # fix round 1, amended ruling: no stop(), the cold pass and a1 are kept
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

  local({
    local_mocked_bindings(.lsScoreInfo = function(Y, A, ...) {
      list(U = 5, I = -1, excluded = integer(0), n = nrow(A))
    })
    expect_warning(j <- .polishSharedLS(Y, prob, pol), "not positive")
    expect_identical(j$a1, 1)
    expect_identical(j$pol, pol)
    expect_identical(j$iterations, 0L)
    expect_false(j$converged)
    expect_true(is.na(j$se))
  })

  # and when it turns non-positive after a step: the cold pass is still what
  # is returned, with the steps taken recorded
  real <- .lsScoreInfo
  calls <- 0L
  local({
    local_mocked_bindings(.lsScoreInfo = function(Y, A, ...) {
      calls <<- calls + 1L
      if (calls == 1L) real(Y, A, ...) else
        list(U = 5, I = -1, excluded = integer(0), n = nrow(A))
    })
    expect_warning(j <- .polishSharedLS(Y, prob, pol), "not positive")
    expect_identical(j$a1, 1)
    expect_identical(j$pol, pol)
    expect_identical(j$iterations, 1L)
    expect_false(j$converged)
  })

  # through polishSpaNorm(): the fit is the fixed-a1 polish, with the record
  spe <- .polish_spe()
  fx <- S4Vectors::metadata(polishSpaNorm(spe, verbose = FALSE))$SpaNorm
  local({
    local_mocked_bindings(.lsScoreInfo = function(Y, A, ...) {
      list(U = 5, I = -1, excluded = integer(0), n = nrow(A))
    })
    expect_warning(out <- polishSpaNorm(spe, ls = "joint", verbose = FALSE),
                   "not positive")
    f <- S4Vectors::metadata(out)$SpaNorm
    s <- .polishSlot(f)$settings
    expect_identical(s$ls, "joint")
    expect_false(s$ls.converged)
    expect_identical(s$ls.iterations, 0L)
    expect_identical(s$a1, s$a1.input)
    expect_identical(f$gmean, fx$gmean)
    expect_identical(f$alpha, fx$alpha)
    expect_identical(f$psi, fx$psi)
    expect_identical(.polishSlot(f)$genes, .polishSlot(fx)$genes)
  })
})

# ---- Task 9b: tol.ls default lowered from 1e-6 to 1e-3 (2026-09-26) --------

test_that("the default tol.ls converges within the cap and matches a tight stop", {
  # Measured on four real YTMA cores (Task 10, jobs 28973894-98): every joint
  # fit hit the ls.maxit = 10 cap and warned at the old 1e-6 default, although
  # a1 had settled by step 3 -- 1e-6 asks for a1 within a millionth of its own
  # profiled SE, far below what the inner per-gene polish can resolve. This
  # checks the new default (1e-3, a thousandth of the SE) both converges
  # inside the cap and lands within its own tolerance of a1's value under a
  # far tighter stop.
  #
  # tol.ls has no argument of polishSpaNorm() to reach: it is forwarded
  # through `...` at every level down to .polishSharedLS(), but the SAME `...`
  # reaches the cold-pass polishNB() call first, which does not accept
  # tol.ls and errors ("unused argument") before .polishSharedLS() is ever
  # called (confirmed by trying it). So the tight run below mocks the shared
  # default `.LS_TOL` instead of passing tol.ls through the public entry
  # point; .polishSharedLS() is not called directly because the point of the
  # test is the default *as delivered by polishSpaNorm()*, not the internal
  # step in isolation.
  spe <- .polish_spe()
  out <- polishSpaNorm(spe, ls = "joint", verbose = FALSE)
  s <- .polishSlot(S4Vectors::metadata(out)$SpaNorm)$settings
  expect_identical(s$ls.tol, 1e-3)
  expect_true(s$ls.converged)
  expect_lt(s$ls.iterations, 10L)

  local_mocked_bindings(.LS_TOL = 1e-9, .package = "SpaNorm")
  tight <- polishSpaNorm(spe, ls = "joint", verbose = FALSE)
  st <- .polishSlot(S4Vectors::metadata(tight)$SpaNorm)$settings
  expect_true(st$ls.converged)
  expect_lt(abs(s$a1 - st$a1), 1e-3 * s$ls.se)
})

# ---- final review I-2: the batch-level hold-out under ls = "joint" ---------

test_that("joint: a gene with no counts in a batch level is held out and shares the new a1", {
  # the same rule as an all-zero gene (rulings 2 and 8): not polished, does not
  # inform a1, keeps gmean, alpha[, -1] and psi bit for bit, takes the new a1
  spe <- .polish_batch_spe()
  g <- .polish_batch_gene
  expect_warning(out <- polishSpaNorm(spe, ls = "joint", verbose = FALSE),
                 "^1 gene with no counts in some batch level")
  f <- S4Vectors::metadata(out)$SpaNorm
  f0 <- S4Vectors::metadata(out)$SpaNormUnpolished
  s <- .polishSlot(f)$settings
  pol <- .polishSlot(f)$genes
  expect_true(s$ls.converged)
  expect_false(isTRUE(all.equal(s$a1, s$a1.input)))          # a1 moved
  expect_false(pol$polished[g])
  expect_identical(pol$held_out[g], "zero-batch-level")
  expect_true(all(pol$polished[-g]))
  expect_identical(f$gmean[g], f0$gmean[g])
  expect_identical(f$alpha[g, -1], f0$alpha[g, -1])
  expect_identical(f$psi[g], f0$psi[g])
  expect_true(length(unique(f$alpha[, 1])) == 1)             # still shared
  expect_identical(f$alpha[g, 1], s$a1)
  # the sums ran over the polished genes only: the score is zero over them
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  U <- .joint_a1_score(f, Y, genes = setdiff(seq_len(nrow(Y)), g))
  expect_lt(abs(U) * s$ls.se, 1e-3)
})
