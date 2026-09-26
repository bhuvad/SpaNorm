# polishSpaNorm(ls = "joint"): the pooled, profiled Newton step on the
# library-size coefficient a1 that every gene of a SpaNorm model shares
# (design spec 2026-09-25-polish-to-spanorm, section 4).
#
# The joint objective is the TOTAL penalised log-likelihood over the polished
# genes, sum_g l_g(a1, A_g), with a1 unpenalised. At a fixed a1 the genes are
# independent, which is polishNB(). a1 is one scalar, so its step is a Newton
# step on the profile log-likelihood p(a1) = sum_g max_{A_g} l_g(a1, A_g):
#
#   U =  p'(a1) = sum_g w1' r_g                                (each A_g at its optimum)
#   I = -p''(a1) = sum_g [ (w1^2)' D_g - c_g' H_g^-1 c_g ]     (A_g profiled out)
#
# with D_g = mu/(1 + psi mu), r_g = (y - mu)/(1 + psi mu), c_g = X'(D_g * w1)
# and H_g = X' D_g X + diag(pen), all at the current polished fit. I is the
# Fisher (expected) information, the same weights the per-gene Newton uses,
# and H_g^-1 c_g comes from nbNewtonSolver()'s factorisation of D_g, the one
# the polish uses. Each gene's term is the Schur complement of H_g in the
# augmented information of (A_g, a1), so it is >= 0: I > 0 unless w1 lies in
# the span of X on every gene's support (a1 is then not identified).

#' The pooled score and profiled information for the shared a1
#'
#' Blocked over genes (\code{block} at a time), so the dense genes x cells
#' matrices it forms are bounded by the block, never the whole gene set. A
#' gene whose penalised information \code{H_g} is singular (the solver returns
#' \code{NULL}, or a non-finite solution) is left out of BOTH sums and
#' reported in \code{excluded}: a zero correction instead would add its
#' unprofiled \code{(w1^2)' D_g} to \code{I} and overstate the information.
#'
#' @param Y counts, genes x cells (dense, sparse or DelayedArray).
#' @param A coefficients, genes x \code{ncol(X)}.
#' @param psi per-gene dispersions.
#' @param a1 the shared library-size coefficient.
#' @param X the per-gene design (\code{.spaNormPolishProblem()$X}).
#' @param w1 the log library size over the same cells.
#' @param pen the per-column penalty of \code{X}.
#' @param solver \code{nbNewtonSolver(X, pen)}.
#' @param block genes per block.
#' @return a list with \code{U}, \code{I}, \code{excluded} (row indices of
#'   the genes left out) and \code{n} (the genes summed).
#' @noRd
.lsScoreInfo <- function(Y, A, psi, a1, X, w1, pen, solver, block = 500L) {
  U <- 0
  I <- 0
  excluded <- integer(0)
  # (D_g * w1)' X is D_g' (w1 * X): scale the design's rows once
  Xw <- X * w1
  w1sq <- w1^2
  for (b in .chunkGenes(nrow(A), block)) {
    Ab <- A[b, , drop = FALSE]
    Mu <- .muBatch(Ab, X, .offsetRows(a1 * w1, seq_along(b), Ab))
    den <- 1 + psi[b] * Mu               # psi[b] recycles down the columns: row k gets psi[b[k]]
    R <- (as.matrix(Y[b, , drop = FALSE]) - Mu) / den
    D <- Mu / den
    rm(Mu, den)
    u <- as.numeric(R %*% w1)            # per-gene w1' r_g
    rm(R)
    d11 <- as.numeric(D %*% w1sq)        # per-gene (w1^2)' D_g
    C <- D %*% Xw                        # row k is c_g'
    for (k in seq_along(b)) {
      h <- solver$solve(solver$factor(D[k, ]), C[k, ])
      if (is.null(h) || !all(is.finite(h))) {
        excluded <- c(excluded, b[k])
        next
      }
      U <- U + u[k]
      I <- I + d11[k] - sum(C[k, ] * h)
    }
  }
  list(U = U, I = I, excluded = excluded, n = nrow(A) - length(excluded))
}

#' Per-gene penalised log-likelihood at a shared a1, blocked over genes
#'
#' The polish's own objective (\code{.muBatch()}'s floored, unwinsorised mean
#' and \code{.nbLoglikBatch()}), one value per gene.
#' @noRd
.lsLoglik <- function(Y, A, psi, a1, X, w1, pen, block = 500L) {
  out <- numeric(nrow(A))
  for (b in .chunkGenes(nrow(A), block)) {
    Ab <- A[b, , drop = FALSE]
    Mu <- .muBatch(Ab, X, .offsetRows(a1 * w1, seq_along(b), Ab))
    out[b] <- .nbLoglikBatch(as.matrix(Y[b, , drop = FALSE]), Mu, psi[b], Ab, pen)
  }
  out
}

# The joint loop's own stop: the standardised score |U|/sqrt(I) and the cap on
# steps. One place, so the loop's defaults and what a fit records agree.
#
# The standardised score is a1's distance from its optimum in units of its
# own profiled SE (score/sqrt(information) = score * se). 1e-6 asked for a1
# to within a millionth of its SE, far below what the inner per-gene polish
# can resolve: measured on four real YTMA cores (Task 10, jobs
# 28973894-98), all 8 joint fits hit the ls.maxit = 10 cap and warned,
# although a1 had settled by step 3. 1e-3 asks for a1 within a thousandth of
# its SE, still far inside its sampling error.
.LS_MAXIT <- 10L
.LS_TOL <- 1e-3

#' The joint step's result when it cannot inform a1
#'
#' No polished gene (ruling 9), every gene's information singular, or a
#' non-positive profiled information (fix round 1): warn, and return the cold
#' pass \code{pol} and the fit's \code{a1} unchanged, so the fit is the
#' fixed-a1 polish with the joint attempt recorded.
#' @noRd
.lsUninformed <- function(pol, a1, why = "no gene was polished",
                          what = "nothing informs the shared library-size coefficient",
                          iterations = 0L, singular = 0L,
                          maxit.ls = .LS_MAXIT, tol.ls = .LS_TOL) {
  warning("ls = \"joint\": ", why, ", so ", what, "; it is left at the fit's ",
          "value and the genes at their fixed-a1 polish", call. = FALSE)
  list(pol = pol, a1 = a1, iterations = as.integer(iterations), score = NA_real_,
       se = NA_real_, singular = as.integer(singular), converged = FALSE,
       maxit.ls = as.integer(maxit.ls), tol.ls = tol.ls)
}

#' The joint polish: alternate the per-gene polish with a pooled step on a1
#'
#' Starts from a converged per-gene polish at the fit's \code{a1}
#' (\code{pol}, from \code{polishNB()}) and takes Newton steps on the profile
#' log-likelihood of \code{a1} (see the file header), each with a line search
#' on the total penalised log-likelihood after a warm re-polish at the
#' candidate \code{a1}. It stops when the standardised score
#' \code{|U| / sqrt(I)} is below \code{tol.ls}, at \code{maxit.ls} steps, or
#' when no halving of the step keeps the total from falling. A non-positive
#' \code{I} (\code{a1} not identified) returns the cold pass and the fit's
#' \code{a1} with a warning, as when no gene is polished.
#'
#' Only the genes \code{pol} polished take part: the sums, the objective and
#' the re-polish run over them, and a gene the cold pass could not polish
#' keeps its input fit (\code{pol}'s fallback), like an all-zero gene held
#' out upstream. Within an iteration the step, the objective it is judged on
#' and the stop use one gene set, the genes that formed \code{U} and
#' \code{I}: a gene left out for a singular information is re-polished at
#' each candidate but not compared. Under \code{psi.method = "profile"} a candidate is a warm
#' re-polish at the held dispersion, a profile of the dispersion at that
#' converged mean (\code{nbProfilePsi()}) and a warm re-polish at it,
#' polishNB()'s own profile-then-re-polish rule: every gene ends at a zero
#' score at its reported dispersion, which the profiled information assumes.
#' The objective is then the profile likelihood in the dispersion too, and the
#' line search is an ascent in it.
#'
#' @param Y counts of the genes \code{pol} covers (genes x cells over the
#'   polished cells).
#' @param prob \code{.spaNormPolishProblem()}'s list (\code{X}, \code{pen},
#'   \code{a1}, \code{w1}).
#' @param pol \code{polishNB()}'s result at \code{prob$a1}.
#' @param psi.method \code{"fixed"} or \code{"profile"}.
#' @param maxit.ls,tol.ls the cap on steps and the standardised-score stop.
#' @param block genes per block of the score, information and objective.
#' @param verbose logical; one line per step.
#' @param ... passed to every warm \code{polishNB()} (\code{maxit},
#'   \code{tol}, \code{engine}, \code{batch.size}, \code{block.size},
#'   \code{backend}, \code{BPPARAM}); \code{psi.range}, \code{block.size} and
#'   \code{BPPARAM} also reach \code{nbProfilePsi()}, which takes no
#'   tolerance. It precedes the options below, which therefore match only by
#'   their exact names.
#' @return a list with \code{pol} (as \code{polishNB()} returns it, at the
#'   joint \code{a1}: \code{loglik} is each gene's penalised log-likelihood at
#'   the returned fit, and the diagnostics add the warm passes' Newton
#'   iterations to \code{iterations} and their flags to \code{capped} and
#'   \code{singular}), \code{a1}, \code{iterations} (accepted steps),
#'   \code{score} (the final \code{U}), \code{se} (\code{1 / sqrt(I)}, the
#'   profiled standard error of \code{a1}), \code{singular} (the number of
#'   genes left out of \code{U} and \code{I} at some step for a singular
#'   information), \code{converged}, and \code{maxit.ls}/\code{tol.ls} as
#'   used.
#' @noRd
.polishSharedLS <- function(Y, prob, pol, ..., psi.method = "fixed",
                            maxit.ls = .LS_MAXIT, tol.ls = .LS_TOL, block = 500L,
                            verbose = FALSE) {
  # `...` comes FIRST so that every option after it matches only by its exact
  # name. Before, the caller's per-gene maxit/tol partial-matched maxit.ls/
  # tol.ls (fix round 1): the loop ran at the polish's tol and cap, and the
  # warm re-polishes never saw the caller's maxit/tol.
  X <- prob$X
  w1 <- prob$w1
  pen <- prob$pen
  a1 <- prob$a1
  ok <- which(pol$polish$polished)
  uninformed <- function(...) {
    .lsUninformed(pol, prob$a1, ..., maxit.ls = maxit.ls, tol.ls = tol.ls)
  }
  if (!length(ok)) return(uninformed())

  Yo <- Y[ok, , drop = FALSE]
  A <- pol$alpha[ok, , drop = FALSE]
  psi <- as.numeric(pol$psi)[ok]
  solver <- nbNewtonSolver(X, pen)
  dots <- list(...)
  arg <- function(name, default) if (is.null(dots[[name]])) default else dots[[name]]
  psi.range <- arg("psi.range", c(1e-3, 1e3))
  BPPARAM <- arg("BPPARAM", BiocParallel::SerialParam())

  # a candidate a1: re-polish every gene warm from the current coefficients,
  # and under "profile" put the dispersion at the converged mean and re-polish
  repolish <- function(A, psi, a1) {
    off <- a1 * w1
    warm <- function(A, psi) {
      polishNB(Yo, X, A, psi, lambda.a = pen, offset = off, psi.method = "fixed",
               warm = TRUE, ...)
    }
    p <- warm(A, psi)
    it <- p$polish$iterations
    capped <- p$polish$capped
    # a warm pass that fails keeps the gene's previous coefficients
    singular <- p$polish$singular | !p$polish$polished
    if (psi.method == "profile") {
      # the profile is a fixed-iteration bisection: it takes no maxit or tol
      psi <- nbProfilePsi(Yo, X, p$alpha, p$psi, psi.range = psi.range,
                          block.size = dots[["block.size"]], BPPARAM = BPPARAM,
                          offset = off)
      p2 <- warm(p$alpha, psi)
      it <- it + p2$polish$iterations
      capped <- capped | p2$polish$capped
      singular <- singular | p2$polish$singular | !p2$polish$polished
      p <- p2
    }
    list(alpha = p$alpha, psi = as.numeric(p$psi), iterations = it,
         capped = capped, singular = singular)
  }
  loglik <- function(A, psi, a1) .lsLoglik(Yo, A, psi, a1, X, w1, pen, block)
  score <- function(A, psi, a1) .lsScoreInfo(Yo, A, psi, a1, X, w1, pen, solver, block)

  ll_g <- loglik(A, psi, a1)
  si <- score(A, psi, a1)
  excluded <- si$excluded
  extra_it <- integer(length(ok))
  extra_capped <- logical(length(ok))
  extra_singular <- logical(length(ok))
  it <- 0L
  stop_why <- NULL
  # defined only for a positive I: the loop returns before stepping otherwise
  zscore <- function(si) {
    if (is.finite(si$I) && si$I > 0) abs(si$U) / sqrt(si$I) else NA_real_
  }
  if (verbose) {
    message(sprintf("  joint library size: a1 = %.6g, |U|/sqrt(I) = %.3g",
                    a1, zscore(si)))
  }
  repeat {
    if (si$n == 0L) {
      # every polished gene's information is singular: nothing informs a1
      if (it == 0L) {
        return(uninformed("every polished gene's information is singular",
                          singular = length(excluded)))
      }
      stop_why <- "every polished gene's information became singular"
      break
    }
    if (!is.finite(si$I) || si$I <= 0) {
      # the profiled information is a sum of Schur complements, so >= 0; a
      # non-positive value means w1 is (numerically) in the span of the
      # per-gene design and a1 is not identified. The steps so far were taken
      # on a meaningless curvature, so the fixed-a1 polish is what is returned.
      return(uninformed(
        sprintf("the profiled information for a1 is not positive (%s)", format(si$I)),
        what = paste0("a1 is not identified (the log library size may be ",
                      "collinear with the model's other terms)"),
        iterations = it, singular = length(excluded)))
    }
    if (zscore(si) <= tol.ls) break
    if (it >= maxit.ls) {
      stop_why <- sprintf("maxit.ls = %d reached", maxit.ls)
      break
    }
    # The step, the objective it is judged on and the stop all refer to ONE
    # gene set: the genes that formed U and I this iteration. A gene left out
    # for a singular information is still re-polished at each candidate, but
    # it is not in the comparison (fix round 1: kept in, a gene that cannot
    # follow a1 stalled the search, or stopped it at a point that was not the
    # maximum of the objective being compared).
    in_set <- setdiff(seq_along(ok), si$excluded)
    ll <- sum(ll_g[in_set])
    step <- si$U / si$I
    accepted <- FALSE
    for (h in 0:10) {
      a1n <- a1 + step / 2^h
      cand <- repolish(A, psi, a1n)
      lln_g <- loglik(cand$alpha, cand$psi, a1n)
      lln <- sum(lln_g[in_set])
      if (is.finite(lln) && lln >= ll - 1e-9 * abs(ll)) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      stop_why <- "no step along the Newton direction raised the total likelihood"
      break
    }
    it <- it + 1L
    a1 <- a1n
    A <- cand$alpha
    psi <- cand$psi
    ll_g <- lln_g
    extra_it <- extra_it + cand$iterations
    extra_capped <- extra_capped | cand$capped
    extra_singular <- extra_singular | cand$singular
    si <- score(A, psi, a1)
    excluded <- union(excluded, si$excluded)
    if (verbose) {
      message(sprintf("  joint library size: step %d, a1 = %.6g, |U|/sqrt(I) = %.3g",
                      it, a1, zscore(si)))
    }
  }
  converged <- is.null(stop_why)
  if (!converged) {
    warning(sprintf(paste0("ls = \"joint\": the step on the shared library-size ",
                           "coefficient stopped with |U|/sqrt(I) = %.3g, above ",
                           "tol.ls = %g (%s)"),
                    zscore(si), tol.ls, stop_why), call. = FALSE)
  }

  pol$alpha[ok, ] <- A
  pol$psi[ok] <- psi
  pol$loglik[ok] <- ll_g
  pol$polish$iterations[ok] <- pol$polish$iterations[ok] + extra_it
  pol$polish$capped[ok] <- pol$polish$capped[ok] | extra_capped
  pol$polish$singular[ok] <- pol$polish$singular[ok] | extra_singular
  list(pol = pol, a1 = a1, iterations = it, score = si$U,
       se = if (is.finite(si$I) && si$I > 0) 1 / sqrt(si$I) else NA_real_,
       singular = length(excluded), converged = converged,
       maxit.ls = as.integer(maxit.ls), tol.ls = tol.ls)
}
