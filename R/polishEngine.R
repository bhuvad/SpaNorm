# The per-gene polish engine, moved from spiDE (R/polish.R): .MU_FLOOR,
# .nbPenLoglik and .polishGene -- damped Newton on one gene's own penalised
# negative binomial log-likelihood, with a profile-ML dispersion. The batched
# engine that computes the same thing for a block of genes is in
# R/polishEngineBatch.R; this one stays as its reference implementation and
# test oracle. The header below is spiDE's, kept as written.
#
# Per-gene convergence of the penalised NB fit.
#
# SpaNorm::fitNB fits every gene in one IRLS loop: it shares a single
# gene-averaged cell weight vector across genes, decides step-halving and
# convergence on the AGGREGATE log-likelihood, and clamps coefficient columns
# across genes. That is what makes a 13,000-gene fit affordable, and for the
# great majority of genes it is indistinguishable from the per-gene optimum. For
# a bright, cell-type-restricted gene -- whose own working weights look nothing
# like the average -- the aggregate criterion is met long before that gene's own
# score is zero: measured on the YTMA cohort, every gene in the top 5% by
# expression sat 1-4 production standard errors from its own penalised-NB
# optimum, with log-likelihood gaps of 1e4-1e6 and an edgeR dispersion ~1.6x too
# large.
#
# This stage removes the fit from the inference question. It runs AFTER fitNB
# (so the cross-gene dispersion moderation still happens on the whole gene set,
# which is why genes must not be blocked at fit time), and because a polished
# gene depends only on its own counts it is blockable and parallel -- the same
# split .blockedInference() uses.

# The fitted mean is floored here and at inference, at the same value. A
# linear predictor below about -745 underflows exp() to exactly 0, and the
# negative binomial quantities built from it then divide by zero: the Pearson
# working dispersion (y - mu)^2 / (mu + psi mu^2) becomes 0/0 = NaN, which
# propagates to the standard error, the statistic and every gene-set test
# downstream. exp(-30) is far below any mean the model can meaningfully
# estimate, so flooring there is a numerical guard, not a statistical clamp --
# and applying the SAME floor in both stages keeps the polish and the inference
# on one mean function, which is the point of not winsorising here.
.MU_FLOOR <- exp(-30)

#' The floor on a fitted negative binomial mean
#'
#' The polish engine floors every fitted mean at this value, and a caller that
#' rebuilds the mean from polished coefficients (for inference, say) should
#' floor at the same value, so both stages work on one mean function. It is a
#' numerical guard against \code{exp()} underflowing to exactly zero, not a
#' statistical clamp: \code{exp(-30)} is far below any mean the model can
#' meaningfully estimate.
#'
#' @return \code{exp(-30)}.
#' @examples
#' nbMuFloor()
#' @keywords internal
#' @export
nbMuFloor <- function() .MU_FLOOR

#' Penalised NB log-likelihood of one gene
#'
#' @param y counts (length ncells).
#' @param mu the fitted mean (length ncells).
#' @param psi the NB dispersion (scalar).
#' @param a the coefficients (length ncol(W)).
#' @param pen the per-column ridge penalty.
#' @return a numeric scalar.
#' @importFrom stats dnbinom
#' @noRd
.nbPenLoglik <- function(y, mu, psi, a, pen) {
  sum(stats::dnbinom(y, size = 1 / psi, mu = mu, log = TRUE)) -
    0.5 * sum(pen * a^2)
}

#' Converge one gene to its own penalised NB optimum
#'
#' Damped Newton on the penalised log-likelihood at fixed \code{psi}, then
#' profile ML for \code{psi} at the converged mean, then a re-polish -- twice.
#' The information matrix is reused for up to three consecutive steps (only the
#' score is recomputed), since it changes slowly and each rebuild is essentially
#' the whole cost of an iteration.
#'
#' Starting from fitNB's coefficients is right for almost every gene, but for a
#' degenerate fit (a fitted mean below exp(-10) somewhere) Newton from that
#' point DIVERGES -- measured: fitted log-means reaching +57 to +358 and
#' predicted one-step gains of 1e6-1e8 against actual gains of 1e2-1e5. Those
#' genes restart from a sane point instead: the per-cell-type log mean, every
#' other coefficient zero, which converges in 5-29 iterations.
#'
#' @param y counts for this gene (length ncells).
#' @param W the design.
#' @param a0,psi0 fitNB's coefficients and dispersion for this gene.
#' @param pen the per-column ridge penalty.
#' @param solver a \code{.newtonSolver()} for this \code{W} and \code{pen}.
#' @param maxit,tol iteration cap and relative log-likelihood tolerance.
#' @param start.cols a logical over the columns of \code{W} marking the cell-type
#'   intercepts (from the design's covtype tags).
#' @param psi.range the search interval for the profile-ML dispersion.
#' @param warm logical; \code{a0}/\code{psi0} are a CONVERGED fit at a nearby
#'   penalty (the re-polish after a variance-component step). A warm polish is
#'   a few damped Newton steps at the held dispersion: no profile-psi search
#'   (its optimum moves at second order in the penalty change) and no
#'   log-mean restart check -- a converged fit legitimately has fitted
#'   log-means below -10 where a gene is absent from a cell type, and treating
#'   that like fitNB's degenerate output threw ~100 of 769 genes back to the
#'   sane start on every re-polish pass at bandwidth 10 on the cohort (3-5 in
#'   the cold pass), each redoing a full cold polish.
#' @return a list with \code{alpha}, \code{psi}, \code{loglik},
#'   \code{iterations}, \code{restarted}, \code{capped}, \code{singular},
#'   \code{psi_bound} and \code{polished}.
#' @importFrom stats optimize
#' @noRd
.polishGene <- function(y, W, a0, psi0, pen, solver, maxit = 50L, tol = 1e-8,
                        start.cols = NULL, psi.range = c(1e-3, 1e3),
                        psi.method = c("profile", "fixed"), warm = FALSE) {
  psi.method <- match.arg(psi.method)
  restarted <- FALSE
  singular <- FALSE

  # the sane start: cell-type (or, absent a cell-type block, overall) log means.
  # `start.cols` comes from the design's own covtype tags. It used to be recovered
  # by a regex on colnames(W), which is a second, weaker parser of a convention
  # .tagCovtype() already owns: a user covariate literally named "CellTypeScore"
  # matched it and was assigned a log mean as though it were an indicator, and a
  # cell-type label containing ":" did not match at all.
  sane_start <- function() {
    a <- numeric(ncol(W))
    ct <- if (is.null(start.cols)) integer(0) else which(start.cols)
    if (length(ct)) {
      for (j in ct) {
        cells <- W[, j] != 0
        a[j] <- if (any(cells)) log(mean(y[cells]) + 1e-3) else 0
      }
    } else {
      a[1] <- log(mean(y) + 1e-3)
    }
    a
  }

  newton <- function(a, psi, maxit) {
    mu <- pmax(as.numeric(exp(W %*% a)), .MU_FLOOR)
    ll <- .nbPenLoglik(y, mu, psi, a, pen)
    it <- 0L
    converged <- FALSE
    stale <- 0L
    w <- NULL
    fac <- NULL
    while (it < maxit) {
      it <- it + 1L
      s <- as.numeric(crossprod(W, (y - mu) / (1 + psi * mu))) - pen * a
      if (is.null(w) || stale >= 3L) {
        w <- mu / (1 + psi * mu)
        # Factor here and only here: the weights are what the information
        # depends on, and between refreshes the same factorisation is exact.
        # A solver without $factor() -- the two-function contract this used to
        # have, which callers and test stubs may still implement -- keeps
        # working: $solve() accepts the weights directly and factors them
        # itself, which is what every step used to do.
        fac <- if (is.function(solver$factor)) solver$factor(w) else w
        stale <- 0L
      }
      d <- solver$solve(fac, s)
      # solve() only ERRORS below rcond ~1e-7; between that and well-conditioned
      # it returns a finite but numerically meaningless answer, which the line
      # search can accept because a badly scaled step in roughly the right
      # direction still raises the objective. A NaN/Inf right-hand side does not
      # error either. Treat both as singular rather than letting a wrong number
      # through as a converged coefficient.
      if (is.null(d) || !all(is.finite(d))) {
        singular <<- TRUE
        break
      }
      step <- 1
      ok <- FALSE
      halvings <- 0L
      while (step > 1e-6) {
        a1 <- a + step * d
        mu1 <- pmax(as.numeric(exp(W %*% a1)), .MU_FLOOR)
        ll1 <- .nbPenLoglik(y, mu1, psi, a1, pen)
        if (is.finite(ll1) && ll1 >= ll - 1e-9 * abs(ll)) {
          ok <- TRUE
          break
        }
        step <- step / 2
        halvings <- halvings + 1L
      }
      if (!ok) {
        # a stale information matrix can give a bad direction; rebuild it once
        # before giving up
        if (stale > 0L) {
          w <- mu / (1 + psi * mu)
          fac <- if (is.function(solver$factor)) solver$factor(w) else w
          stale <- 0L
          next
        }
        break
      }
      # a hard line search or a stale matrix both call for a rebuild next step
      stale <- if (halvings > 2L) 3L else stale + 1L
      gain <- ll1 - ll
      a <- a1
      mu <- mu1
      ll <- ll1
      if (gain < tol * abs(ll)) {
        converged <- TRUE
        break
      }
    }
    list(a = a, mu = mu, ll = ll, it = it, converged = converged)
  }

  # Profile ML for psi at a fixed mean. The optimiser needs a bounded interval,
  # so an under-dispersed or near-empty gene lands ON a bound -- measured: a
  # Poisson gene returns 0.0046 and an all-zero gene 976, neither of which is an
  # estimate. Storing a bound as though it were one is worse than not polishing
  # the dispersion at all, because on the fixed-effects path psi scales the
  # standard error directly. Report it instead, and let the caller keep fitNB's
  # moderated value.
  psi_ml <- function(mu) {
    lo <- log(psi.range[1]); hi <- log(psi.range[2])
    o <- stats::optimize(function(lp) {
      -sum(stats::dnbinom(y, size = 1 / exp(lp), mu = mu, log = TRUE))
    }, c(lo, hi))
    edge <- (o$minimum - lo) < 1e-3 * (hi - lo) ||
      (hi - o$minimum) < 1e-3 * (hi - lo)
    list(psi = exp(o$minimum), at_bound = edge)
  }

  fallback <- function(why) {
    # Hand back fitNB's own estimate, NOT the sane start. Returning the sane
    # start would replace a usable fit with cell-type log means and exact zeros
    # on every tested coefficient, which inference then reports as t = 0 with a
    # finite SE -- a confident null for a gene that was never converged.
    list(alpha = a0, psi = psi0, loglik = NA_real_, iterations = 0L,
         restarted = FALSE, capped = FALSE, singular = identical(why, "singular"),
         psi_bound = FALSE, polished = FALSE)
  }
  # a degenerate START is a fitted log-mean below -10 at a cell that has a
  # count; at a zero-count cell it is the likelihood's own direction (a gene
  # absent from a cell type puts that unpenalised intercept toward -Inf), and
  # restarting such a gene redoes a converged fit for nothing
  degenerate <- function(a) {
    if (!all(is.finite(a))) return(TRUE)
    lp <- as.numeric(W %*% a)
    any(y > 0) && min(lp[y > 0]) < -10
  }
  a <- a0
  psi <- psi0
  if (warm) {
    if (!all(is.finite(a))) return(fallback("nonfinite"))
    f <- newton(a, psi, maxit)
    if (singular || !all(is.finite(f$a)) || !is.finite(f$ll)) {
      return(fallback(if (singular) "singular" else "nonfinite"))
    }
    return(list(alpha = f$a, psi = psi, loglik = f$ll, iterations = f$it,
                restarted = FALSE, capped = !f$converged, singular = FALSE,
                psi_bound = FALSE, polished = TRUE))
  }
  if (degenerate(a)) {
    a <- sane_start()
    restarted <- TRUE
  }
  f <- newton(a, psi, maxit)
  if (!restarted && (singular || !all(is.finite(f$a)) || max(f$mu) > 1e10)) {
    singular <- FALSE
    restarted <- TRUE
    f <- newton(sane_start(), psi, maxit)
  }
  # both starts failed: keep fitNB's fit rather than an unconverged guess
  if (singular || !all(is.finite(f$a)) || !is.finite(f$ll)) {
    return(fallback(if (singular) "singular" else "nonfinite"))
  }
  converged <- f$converged
  it_total <- f$it
  psi_bound <- FALSE
  # "fixed" keeps fitNB's cross-gene moderated dispersion and converges
  # only the mean under it; "profile" re-estimates psi by profile ML at the
  # converged mean and re-polishes, twice
  if (psi.method == "profile") for (k in 1:2) {
    pm <- psi_ml(f$mu)
    if (pm$at_bound) {
      # the dispersion is not identified for this gene; keep the moderated one
      psi_bound <- TRUE
      break
    }
    psi <- pm$psi
    f2 <- newton(f$a, psi, 20L)
    it_total <- it_total + f2$it
    converged <- converged && f2$converged
    f <- f2
  }
  list(alpha = f$a, psi = psi, loglik = f$ll, iterations = it_total,
       restarted = restarted, capped = !converged, singular = singular,
       psi_bound = psi_bound, polished = TRUE)
}
