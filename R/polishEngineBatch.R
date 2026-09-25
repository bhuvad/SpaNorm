# The batched polish engine, moved from spiDE (R/polish-batch.R, all of it),
# plus a copy of spiDE's .covBatchSize() as .polishCovBatchSize() with the
# constants it reads (spiDE keeps its own copy for inference). The header below
# is spiDE's, kept as written.
#
# The per-gene Newton, restructured so a block of genes shares every read of
# the design.
#
# Nothing here is a new estimator. .polishGene() computes the per-gene penalised
# NB optimum and .newtonSolver() already absorbs the nested indicator block; this
# file computes the same thing for B genes at once, and .polishGene() stays in
# the tree as the reference implementation and the test oracle.
#
# What batching buys, measured at the cohort shape (n = 77,454, px = 398;
# FINDINGS.md 2026-09-15). Every gene in a block shares W, so the two
# design-sized matrix-vector products per Newton step -- the linear predictor
# and the score -- become one GEMM each:
#
#   B genes        separate matvecs   one GEMM   per gene
#   1                    0.07 s         0.08 s     79 ms
#   64                   7.38 s         0.11 s    1.7 ms
#   256                 29.33 s         0.22 s    0.8 ms
#
# The gram is NOT batched here. Each gene has its own working weights, so it has
# its own information matrix, and on the CPU that is irreducible; it stays a
# per-gene call into the solver. After the factorisation memoisation it is paid
# on roughly one step in three, and it is what the device path exists to attack.
#
# The dispersion search is also left per gene. At 0.16 s per optimize() against
# ~18 s of gram per gene it is a rounding error, and keeping stats::optimize()
# keeps this a pure restructuring: no deliberate numerical divergence to argue
# about while the control flow is being rebuilt.

# Live gene x cell matrices inside one batched Newton iteration: the counts
# slice, mu, the residual/weight matrix, the candidate coefficients' mean, the
# accepted mean, plus headroom for R's copy-on-modify. Sized by inspection, so
# it errs high.
POLISH_GENE_CELL_MATS <- 6

#' Refuse a single-precision device
#'
#' \code{getBackendDtype()} is float64 on CUDA but float32 on MPS. A
#' penalised Newton solve over hundreds of columns in single precision is not
#' defensible -- the Cholesky of a near-singular information matrix is exactly
#' where the last digits matter -- and the failure would be silent, looking
#' like a convergence problem rather than a precision one. So it errors, with
#' the lever in the message, rather than warning and proceeding.
#'
#' @param dtype the backend dtype, as a string or a torch dtype.
#' @return invisibly TRUE, or an error.
#' @noRd
.requireFloat64 <- function(dtype = getBackendDtype()) {
  nm <- tolower(paste(as.character(dtype), collapse = " "))
  if (!grepl("double|float64", nm)) {
    stop("the polish stage refuses a single precision device (dtype ", nm,
         ").\n  A penalised Newton solve over hundreds of columns in float32 ",
         "is not defensible, and the failure would look like non-convergence ",
         "rather than lost precision.\n  Use backend = \"cpu\", or a device ",
         "with float64 (CUDA has it; MPS does not).", call. = FALSE)
  }
  invisible(TRUE)
}

#' The batched mean, on either backend
#'
#' \code{exp(A W')} floored at \code{.MU_FLOOR}, for a block of genes at once.
#' The floor is not cosmetic: a degenerate start can drive a linear predictor to
#' -50 and the log-likelihood to -Inf, and the per-gene engine has always
#' clamped here.
#'
#' @param A the block's coefficients, \code{genes x p}.
#' @param W the design, \code{cells x p}. Matching types: both matrices, or
#'   both torch tensors.
#' @return \code{genes x cells}, the same type as the inputs.
#' @noRd
.muBatch <- function(A, W) {
  if (is_torch_tensor(A)) {
    eta <- torch::torch_matmul(A, W$transpose(1, 2))
    return(torch::torch_clamp(torch::torch_exp(eta), min = .MU_FLOOR))
  }
  pmax(exp(A %*% t(W)), .MU_FLOOR)
}

#' The batched penalised NB log-likelihood, on either backend
#'
#' \code{sum_cells dnbinom(y; size = 1/psi, mu) - 0.5 * a' diag(pen) a}, per
#' gene. On the torch branch the NB log-pmf is written out, because
#' \code{dnbinom()} has no tensor equivalent:
#'
#'   lgamma(y + r) - lgamma(r) - lgamma(y + 1) + r log(r/(r+mu)) + y log(mu/(r+mu))
#'
#' with \code{r = 1/psi}. A zero count contributes nothing through the last
#' term (\code{mu} is floored strictly positive, so the log is finite and the
#' factor is zero) and \code{lgamma(1) = 0} through the third, which is what
#' keeps an all-zero gene finite rather than \code{NaN}.
#'
#' \code{psi} and \code{pen} may stay plain R vectors even when the counts are
#' tensors: they are per gene and per column, small, and holding them on the
#' host avoids a transfer per line-search trial.
#'
#' @param Y counts, \code{genes x cells}.
#' @param M the fitted means, \code{genes x cells}.
#' @param psi per-gene dispersion, length \code{genes}.
#' @param A the block's coefficients, \code{genes x p}.
#' @param pen the length-\code{p} ridge penalty.
#' @return a length-\code{genes} vector (or tensor) of penalised
#'   log-likelihoods.
#' @noRd
.nbLoglikBatch <- function(Y, M, psi, A, pen) {
  if (is_torch_tensor(M)) {
    dt <- M$dtype
    dev <- M$device
    asT <- function(x) {
      if (is_torch_tensor(x)) x else
        torch::torch_tensor(as.numeric(x), dtype = dt, device = dev)
    }
    Yt <- asT(Y)
    psi_t <- asT(psi)
    # psi = 0 is size = 1/psi = Inf, which dnbinom() treats as POISSON. The
    # written-out log-pmf does not get that for free: r = Inf gives
    # lgamma(y + Inf) - lgamma(Inf) = Inf - Inf = NaN and r log(r/(r+mu)) =
    # Inf log(1) = NaN, so every gene's likelihood comes back NaN and the
    # engine -- which drops a gene whose loglik is not finite -- polishes
    # NOTHING, silently. fitSpiDE() returns psi = 0 on the toy fixture, so this
    # is the common case on small data, not a corner.
    pois <- psi_t == 0
    # clamped so the NB branch stays finite where it will be discarded
    r <- (1 / torch::torch_clamp(psi_t, min = 1e-300))$unsqueeze(2)
    rm_ <- r + M
    ll_nb <- torch::torch_sum(
      torch::torch_lgamma(Yt + r) - torch::torch_lgamma(r) -
        torch::torch_lgamma(Yt + 1) +
        r * torch::torch_log(r / rm_) + Yt * torch::torch_log(M / rm_), dim = 2)
    ll_pois <- torch::torch_sum(
      Yt * torch::torch_log(M) - M - torch::torch_lgamma(Yt + 1), dim = 2)
    ll <- torch::torch_where(pois, ll_pois, ll_nb)
    pen_t <- asT(pen)
    return(ll - 0.5 * torch::torch_matmul(asT(A)$pow(2), pen_t))
  }
  rowSums(stats::dnbinom(Y, size = 1 / psi, mu = M, log = TRUE)) -
    0.5 * as.numeric((A^2) %*% pen)
}

#' The profile dispersion, by fixed-iteration bisection on the score
#'
#' Replaces a per-gene \code{optimize()} with a batched bisection, and the
#' reason is accuracy before it is speed. \code{optimize()}'s default tolerance
#' is \code{.Machine$double.eps^0.25}, about 1.2e-4 ABSOLUTE on a log-interval
#' of width \code{log(1e3) - log(1e-3) = 13.8}, so the dispersion the per-gene
#' engine reports is only determined to ~1e-4. Fifty halvings take that to
#' ~1e-14. Being fixed-iteration it also has no divergent control flow across a
#' batch, which reproducing Brent would have.
#'
#' It bisects the SCORE rather than searching the objective: one
#' \code{digamma} per cell per evaluation against \code{dnbinom}'s three
#' \code{lgamma}. With \code{r = 1/psi},
#'
#'   dl/dr = sum_i digamma(y_i+r) - digamma(r) + log(r/(r+mu_i)) + 1 - (r+y_i)/(r+mu_i)
#'
#' and \code{dl/d log psi = -r dl/dr}, which is positive at the small-psi end
#' of the range for overdispersed counts and negative at the large-psi end, so
#' the root is bracketed. A gene whose optimum lies OUTSIDE the range does not
#' bracket, and is clamped to the endpoint it ran past -- which is then flagged
#' by the same \code{at_bound} rule the per-gene engine applies, verbatim,
#' rather than by a second rule that would have to agree with it.
#'
#' @param Y counts, \code{genes x cells}.
#' @param Mu the fitted means, \code{genes x cells}.
#' @param psi.range the search interval for the dispersion.
#' @param maxit bisection steps. 50 is ~1e-14 on the log interval.
#' @return \code{list(psi, at_bound)}; \code{psi} matches the input's type,
#'   \code{at_bound} is always a plain logical vector.
#' @noRd
.psiProfileBatch <- function(Y, Mu, psi.range = c(1e-3, 1e3), maxit = 50L) {
  lo <- log(psi.range[1])
  hi <- log(psi.range[2])
  edge <- 1e-3 * (hi - lo)

  if (is_torch_tensor(Mu)) {
    dt <- Mu$dtype
    dev <- Mu$device
    Yt <- if (is_torch_tensor(Y)) Y else
      torch::torch_tensor(as.matrix(Y), dtype = dt, device = dev)
    b <- Mu$size(1)
    score <- function(lp) {
      r <- torch::torch_exp(-lp)$unsqueeze(2)
      rm_ <- r + Mu
      g <- torch::torch_sum(
        torch::torch_digamma(Yt + r) - torch::torch_digamma(r) +
          torch::torch_log(r / rm_) + 1 - (r + Yt) / rm_, dim = 2)
      -torch::torch_exp(-lp) * g
    }
    a <- torch::torch_full(c(b), lo, dtype = dt, device = dev)
    z <- torch::torch_full(c(b), hi, dtype = dt, device = dev)
    Sa <- score(a)
    Sz <- score(z)
    for (i in seq_len(maxit)) {
      m <- (a + z) / 2
      pos <- score(m) > 0
      a <- torch::torch_where(pos, m, a)
      z <- torch::torch_where(pos, z, m)
    }
    lp <- (a + z) / 2
    lp <- torch::torch_where(Sa <= 0, torch::torch_full_like(lp, lo), lp)
    lp <- torch::torch_where(Sz >= 0, torch::torch_full_like(lp, hi), lp)
    lpr <- as.numeric(toRMatrix(lp))
    return(list(psi = torch::torch_exp(lp),
                at_bound = (lpr - lo) < edge | (hi - lpr) < edge))
  }

  b <- nrow(Mu)
  score <- function(lp) {
    r <- exp(-lp)
    rm_ <- r + Mu
    -r * rowSums(digamma(Y + r) - digamma(r) + log(r / rm_) + 1 -
                   (r + Y) / rm_)
  }
  a <- rep(lo, b)
  z <- rep(hi, b)
  Sa <- score(a)
  Sz <- score(z)
  for (i in seq_len(maxit)) {
    m <- (a + z) / 2
    pos <- score(m) > 0
    a <- ifelse(pos, m, a)
    z <- ifelse(pos, z, m)
  }
  lp <- (a + z) / 2
  lp[Sa <= 0] <- lo
  lp[Sz >= 0] <- hi
  list(psi = exp(lp), at_bound = (lp - lo) < edge | (hi - lp) < edge)
}

#' Genes per batched Newton, from a memory budget
#'
#' A gene block is sized to bound densification of the counts (2,000 genes); the
#' batched Newton's working set is gene x cell and must be bounded separately.
#' At the cohort's 77,454 cells a 2,000-gene batch would allocate over a
#' terabyte, so the block size cannot be the batch size.
#'
#' @param ncells cells in the design.
#' @param budget bytes available to ONE worker for the batched working set --
#'   per worker, not a total to be divided among them. \code{.covBatchSize()}'s
#'   budget is documented as a total and then claimed independently by every
#'   forked worker, which at 64 workers is a 128 GB claim in a stage already
#'   OOM-killed once at 503 GB. This one says what it means.
#' @return genes per batch, at least 1.
#' @noRd
.polishBatchSize <- function(ncells,
                             # the spiDE option name is the fallback so a
                             # setting made before the move keeps working
                             budget = getOption("SpaNorm.polish.mem.budget",
                                                getOption("spiDE.polish.mem.budget", 1e9))) {
  per_gene <- 8 * as.numeric(ncells) * POLISH_GENE_CELL_MATS
  max(1L, as.integer(floor(budget / per_gene)))
}

#' Converge a block of genes to their own penalised NB optima
#'
#' @param Yb counts for the block (genes x cells), dense.
#' @param W the design (cells x columns).
#' @param A0 fitNB's coefficients for the block (genes x columns).
#' @param psi0 fitNB's dispersions (length nrow(Yb), or recycled).
#' @param pen the per-column ridge penalty.
#' @param solver a \code{.newtonSolver()}.
#' @param maxit,tol,start.cols,psi.range,psi.method,warm as in \code{.polishGene()}.
#'
#' There is no compaction knob. Every batched quantity is built from the active
#' rows, so the work already scales with the active set and there is no
#' masked-but-computed waste to repack. The copy-versus-mask question arises on
#' the device path, where the allocation is the tensor, and belongs there.
#' @return a list of per-gene results in the shape \code{.polishFit()} expects:
#'   \code{alpha} (genes x columns) and the vectors \code{psi}, \code{loglik},
#'   \code{iterations}, \code{restarted}, \code{capped}, \code{singular},
#'   \code{psi_bound}, \code{polished}.
#' @importFrom stats dnbinom optimize
#' @noRd
.polishBatch <- function(Yb, W, A0, psi0, pen, solver, maxit = 50L, tol = 1e-8,
                         start.cols = NULL, psi.range = c(1e-3, 1e3),
                         psi.method = c("profile", "fixed"), warm = FALSE,
                         shared.factor = FALSE, nested = NULL) {
  psi.method <- match.arg(psi.method)
  B <- nrow(Yb)
  p <- ncol(W)
  psi0 <- rep_len(as.numeric(psi0), B)
  tW <- if (is_torch_tensor(W)) W$transpose(1, 2) else t(W)
  has_factor <- is.function(solver$factor)
  # Factorisation accounting, for Phase 2e. The per-gene policy refreshes one
  # gene's information matrix when THAT gene is stale; a single shared tensor
  # factorisation cannot, and must refresh the whole active stack whenever any
  # gene in it is stale. `sync` is what that would have cost, counted as the
  # engine runs, so the design question is answered by measurement rather than
  # by argument. One integer per refresh point; it changes no result.
  n_fac <- 0L
  n_fac_sync <- 0L
  # One factorisation for the active set instead of a list of per-gene ones:
  # the state a device can hold. It is refreshed when ANY active gene is stale,
  # which is a different path from per-gene staleness -- measured at 11% more
  # factorisations at a 128-gene batch (FINDINGS, 2026-09-16) and gated on the
  # objective, not on equality.
  solverB <- if (shared.factor) .newtonSolverBatch(W, pen, nested) else NULL

  # --- batched kernels -------------------------------------------------------
  # psi is length nrow(M): a matrix is column-major, so a per-gene vector
  # recycles down each column and reaches element (i, j) as psi[i]. That is the
  # whole reason these read as if psi were scalar.
  # one definition of each kernel, shared with the tests and with the device
  # path; tW is kept because the base-R branch of .muBatch() transposes W and
  # this loop calls it thousands of times
  mu_of <- function(A) .muBatch(A, W)
  ll_of <- function(Y, M, ps, A) .nbLoglikBatch(Y, M, ps, A, pen)

  # --- the damped Newton, over a set of genes --------------------------------
  # Mirrors .polishGene()'s newton() exactly, including the staleness policy,
  # the 1e-9 acceptance slack, the 1e-6 step floor and the rebuild-retry that
  # consumes an iteration. Every one of those is per gene and is carried as a
  # vector indexed by position within `rows`.
  newton <- function(rows, A, ps, maxit) {
    m <- length(rows)
    Y <- .rowsOf(Yb, rows)
    Ai <- .rowsOf(A, rows)
    pi_ <- ps[rows]
    Mu <- mu_of(Ai)
    # the log-likelihood is length-genes bookkeeping and lives on the host
    ll <- .asHost(ll_of(Y, Mu, pi_, Ai))
    it <- integer(m); conv <- logical(m); sing <- logical(m)
    stale <- integer(m); fac <- vector("list", m)
    st <- NULL; st_rows <- integer(0)
    act <- seq_len(m)

    while (length(act)) {
      it[act] <- it[act] + 1L
      Ya <- .rowsOf(Y, act); Ma <- .rowsOf(Mu, act)
      Aa <- .rowsOf(Ai, act); pa <- pi_[act]
      R <- (Ya - Ma) / (1 + .mulRows(pa, Ma))
      S <- .matmulB(R, W) - .scaleCols(Aa, pen)

      if (shared.factor) {
        # the stack is aligned to `act`; genes only ever LEAVE the active set,
        # so shrinkage is a subset of the stack rather than a rebuild
        if (!is.null(st) && !identical(st_rows, act)) {
          st <- .subsetState(st, match(act, st_rows))
          st_rows <- act
        }
        if (is.null(st) || any(stale[act] >= 3L)) {
          Wt <- Ma / (1 + .mulRows(pa, Ma))
          st <- solverB$factor(Wt)
          st_rows <- act
          stale[act] <- 0L
          n_fac <<- n_fac + length(act)
          n_fac_sync <<- n_fac_sync + length(act)
        }
        D <- solverB$solve(st, S)
        # same verdict as the per-gene path's `!all(is.finite(d))`: NA from a
        # singular slice, but Inf and NaN too
        bad <- !.rowsFinite(D)
        if (any(bad)) {
          D <- .setRows(D, which(bad), .asLike(matrix(0, sum(bad), p), D))
        }
      } else {
      refresh <- vapply(act, function(k) is.null(fac[[k]]), logical(1)) | stale[act] >= 3L
      if (any(refresh)) {
        kk <- act[refresh]
        Wt <- Ma[refresh, , drop = FALSE] / (1 + pa[refresh] * Ma[refresh, , drop = FALSE])
        for (j in seq_along(kk)) {
          w <- Wt[j, ]
          fac[[kk[j]]] <- if (has_factor) solver$factor(w) else w
        }
        stale[kk] <- 0L
        n_fac <<- n_fac + length(kk)
        n_fac_sync <<- n_fac_sync + length(act)
      }

      # the step is per gene: its own information, its own right-hand side
      D <- matrix(0, length(act), p)
      bad <- logical(length(act))
      for (j in seq_along(act)) {
        d <- solver$solve(fac[[act[j]]], S[j, ])
        if (is.null(d) || !all(is.finite(d))) { bad[j] <- TRUE; next }
        D[j, ] <- d
      }
      }
      if (any(bad)) {
        sing[act[bad]] <- TRUE
        act <- act[!bad]
        if (!length(act)) break
        D <- .rowsOf(D, which(!bad)); S <- .rowsOf(S, which(!bad))
        if (shared.factor) {
          st <- .subsetState(st, which(!bad))
          st_rows <- act
        }
      }

      # --- the line search, one trial round for every pending gene ----------
      na <- length(act)
      step <- rep(1, na); halv <- integer(na); ok <- logical(na)
      A1 <- .rowsOf(Ai, act); M1 <- .rowsOf(Mu, act)
      L1 <- ll[act]
      pend <- seq_len(na)
      while (length(pend)) {
        cand <- .rowsOf(Ai, act[pend]) +
          .mulRows(step[pend], .rowsOf(D, pend))
        mu_c <- mu_of(cand)
        # the one deliberate transfer per trial round: the accept test is
        # host-side control flow over length-genes numbers
        ll_c <- .asHost(ll_of(.rowsOf(Yb, rows[act[pend]]), mu_c,
                              pi_[act[pend]], cand))
        ref <- ll[act[pend]]
        acc <- is.finite(ll_c) & ll_c >= ref - 1e-9 * abs(ref)
        if (any(acc)) {
          take <- pend[acc]
          A1 <- .setRows(A1, take, .rowsOf(cand, which(acc)))
          M1 <- .setRows(M1, take, .rowsOf(mu_c, which(acc)))
          L1[take] <- ll_c[acc]
          ok[take] <- TRUE
        }
        fail <- pend[!acc]
        step[fail] <- step[fail] / 2
        halv[fail] <- halv[fail] + 1L
        pend <- fail[step[fail] > 1e-6]
      }

      # --- what each gene does next ------------------------------------------
      # a stale matrix can give a bad direction: rebuild once and retry, which
      # costs an iteration, exactly as the per-gene loop does via `next`
      retry <- !ok & stale[act] > 0L
      if (any(retry)) {
        kk <- act[retry]
        if (shared.factor) {
          # one stack: a retry for any gene rebuilds it for all of them
          Ma2 <- .rowsOf(Mu, act)
          Wt <- Ma2 / (1 + .mulRows(pi_[act], Ma2))
          st <- solverB$factor(Wt)
          st_rows <- act
          stale[act] <- 0L
        } else {
          Wt <- Mu[kk, , drop = FALSE] / (1 + pi_[kk] * Mu[kk, , drop = FALSE])
          for (j in seq_along(kk)) {
            w <- Wt[j, ]
            fac[[kk[j]]] <- if (has_factor) solver$factor(w) else w
          }
          stale[kk] <- 0L
        }
        n_fac <<- n_fac + length(kk)
        n_fac_sync <<- n_fac_sync + length(act)
      }
      stop_now <- !ok & !retry                      # line search exhausted

      leave_converged <- integer(0)
      if (any(ok)) {
        kk <- act[ok]
        gain <- L1[ok] - ll[kk]
        Ai <- .setRows(Ai, kk, .rowsOf(A1, which(ok)))
        Mu <- .setRows(Mu, kk, .rowsOf(M1, which(ok)))
        ll[kk] <- L1[ok]
        stale[kk] <- ifelse(halv[ok] > 2L, 3L, stale[kk] + 1L)
        # the convergence test reads the NEW log-likelihood, as the per-gene
        # loop does (`ll <- ll1` before `if (gain < tol * abs(ll))`)
        done <- gain < tol * abs(ll[kk])
        conv[kk[done]] <- TRUE
        leave_converged <- kk[done]
      }

      act <- setdiff(act, c(act[stop_now], leave_converged))
      act <- act[it[act] < maxit]
      if (shared.factor && length(act) && !identical(st_rows, act)) {
        st <- .subsetState(st, match(act, st_rows))
        st_rows <- act
      }
    }
    list(A = Ai, Mu = Mu, ll = ll, it = it, converged = conv, singular = sing)
  }

  # --- per-gene helpers, vectorised where they are shared --------------------
  sane_start <- function(rows) {
    # built on the host -- it is genes x columns, the one small array here --
    # and moved to the counts' backend on the way out
    A <- matrix(0, length(rows), p)
    ct <- if (is.null(start.cols)) integer(0) else which(start.cols)
    if (length(ct)) {
      for (j in ct) {
        cells <- which(.asHost(.colsOf(W, j)) != 0)
        A[, j] <- if (length(cells)) {
          log(.rowMeansB(.colsOf(.rowsOf(Yb, rows), cells)) + 1e-3)
        } else 0
      }
    } else {
      A[, 1] <- log(.rowMeansB(.rowsOf(Yb, rows)) + 1e-3)
    }
    .asLike(A, Yb)
  }
  degenerate <- function(rows, A) {
    out <- logical(length(rows))
    fin <- .rowsFinite(A)
    out[!fin] <- TRUE
    if (any(fin)) {
      kk <- which(fin)
      Eta <- .matmulB(.rowsOf(A, kk), tW)
      pos <- .rowsOf(Yb, rows[kk]) > 0
      # a gene with no positive count has no such cell: Inf, hence not
      # degenerate, which is the rule the per-gene engine applies
      out[kk] <- .maskedRowMin(Eta, pos) < -10
    }
    out
  }
  # The one deliberate numerical divergence from .polishGene(): a batched
  # fixed-iteration bisection on the score instead of a per-gene optimize().
  # More accurate, not less -- optimize()'s tolerance is ~1.2e-4 on a
  # log-interval of width 13.8 -- and with no divergent control flow across the
  # batch. .polishGene() keeps optimize(), so the two engines differ on psi by
  # about optimize()'s own tolerance; the parity tests take psi out of the
  # comparison with psi.method = "fixed" and measure the psi change alone.
  psi_ml <- function(rows, Mu) {
    .psiProfileBatch(Yb[rows, , drop = FALSE], Mu, psi.range)
  }

  # --- the flow, mirroring .polishGene() ------------------------------------
  # inputs may be tensors; the RETURN is always host -- .polishFit() reads
  # r$alpha[j, ] per gene and @polish is a data frame
  alpha <- .asHostMat(A0); psi <- psi0
  loglik <- rep(NA_real_, B); iters <- integer(B)
  restarted <- capped <- singular <- psi_bound <- polished <- logical(B)

  if (warm) {
    fin <- .rowsFinite(A0)
    if (any(fin)) {
      rows <- which(fin)
      f <- newton(rows, A0, psi0, maxit)
      good <- !f$singular & .rowsFinite(f$A) & is.finite(f$ll)
      kk <- rows[good]
      alpha[kk, ] <- .asHostMat(.rowsOf(f$A, which(good)))
      loglik[kk] <- f$ll[good]
      iters[kk] <- f$it[good]
      capped[kk] <- !f$converged[good]
      polished[kk] <- TRUE
      singular[rows[!good]] <- f$singular[!good]
    }
    return(list(alpha = alpha, psi = psi, loglik = loglik, iterations = iters,
                restarted = restarted, capped = capped, singular = singular,
                psi_bound = psi_bound, polished = polished))
  }

  A <- A0
  deg <- degenerate(seq_len(B), A)
  if (any(deg)) {
    A <- .setRows(A, which(deg), sane_start(which(deg)))
    restarted[deg] <- TRUE
  }
  f <- newton(seq_len(B), A, psi0, maxit)
  fin <- .rowsFinite(f$A)
  mx <- .rowMaxB(f$Mu)
  need <- !restarted & (f$singular | !fin | mx > 1e10)
  if (any(need)) {
    rows <- which(need)
    A <- .setRows(A, rows, sane_start(rows))
    restarted[rows] <- TRUE
    f2 <- newton(rows, A, psi0, maxit)
    # newton() indexes its result by POSITION within `rows`, not by gene id
    f$A <- .setRows(f$A, rows, f2$A)
    f$Mu <- .setRows(f$Mu, rows, f2$Mu)
    f$ll[rows] <- f2$ll; f$it[rows] <- f2$it
    f$converged[rows] <- f2$converged; f$singular[rows] <- f2$singular
  }
  fin <- .rowsFinite(f$A)
  fell <- f$singular | !fin | !is.finite(f$ll)
  singular[fell] <- f$singular[fell]
  keep <- which(!fell)
  if (!length(keep)) {
    return(list(alpha = alpha, psi = psi, loglik = loglik, iterations = iters,
                restarted = restarted, capped = capped, singular = singular,
                psi_bound = psi_bound, polished = polished))
  }

  A <- f$A; Mu <- f$Mu; ll <- f$ll; it_total <- f$it; conv <- f$converged
  ps <- psi0
  if (psi.method == "profile") {
    live <- keep
    for (k in 1:2) {
      if (!length(live)) break
      pm <- psi_ml(live, .rowsOf(Mu, live))
      hit <- pm$at_bound
      psi_bound[live[hit]] <- TRUE
      live <- live[!hit]
      if (!length(live)) break
      ps[live] <- .asHost(pm$psi)[!hit]
      f2 <- newton(live, A, ps, 20L)
      A <- .setRows(A, live, f2$A)
      Mu <- .setRows(Mu, live, f2$Mu)
      ll[live] <- f2$ll
      it_total[live] <- it_total[live] + f2$it
      conv[live] <- conv[live] & f2$converged
    }
  }
  alpha[keep, ] <- .asHostMat(.rowsOf(A, keep))
  psi[keep] <- ps[keep]
  loglik[keep] <- ll[keep]
  iters[keep] <- it_total[keep]
  capped[keep] <- !conv[keep]
  polished[keep] <- TRUE
  structure(list(alpha = alpha, psi = psi, loglik = loglik, iterations = iters,
                 restarted = restarted, capped = capped, singular = singular,
                 psi_bound = psi_bound, polished = polished),
            factorisations = c(pergene = n_fac, sync = n_fac_sync))
}

POLISH_TENSOR_MULT_COV <- 6
# fraction of the memory budget each of the two independently-bounded stages
# (NB math per gene block, covariance per sub-batch) may claim
POLISH_BUDGET_FRACTION <- 0.5
# default cap on the covariance stack when there is no GPU budget to consult
# (the CPU path); override with options(SpaNorm.cov.mem.budget = <bytes>)
POLISH_COV_MEM_BUDGET_CPU <- 2e9

#' Number of genes per covariance sub-batch
#'
#' Bounds the \code{(batch, p, p)} Gram/inverse stack (and, on the GPU, the
#' \code{(batch, ncells, p)} weighted design feeding it) inside
#' \code{.waldCauchyBlock()}. This is the knob that makes wide mixed-effects
#' designs tractable: peak covariance memory is linear in \code{p} and scales
#' with this batch, so a design too wide to process all at once is handled by
#' shrinking the sub-batch rather than failing.
#'
#' Applies on \strong{both} backends. The CPU path needs it just as much as
#' the GPU path -- a single block of 13,348 genes at \code{p = 4906} would
#' otherwise try to allocate a \code{(13348, 4906, 4906)} array -- and unlike
#' \code{.inferenceBlockSize()} it therefore does not return NULL for CPU.
#'
#' @param ncells number of cells.
#' @param p the covariance dimension (\code{ncol} of the Gram design).
#' @param backend the resolved backend.
#' @param gpu.mem.budget \code{NULL} (auto-detect) or a budget in bytes; only
#'   consulted on the GPU path.
#' @param nworkers how many workers will evaluate this concurrently. The budget
#'   is a figure for the machine (or the device), but \code{.waldCauchyBlock()}
#'   runs inside \code{bplapply()} and each forked worker claims it
#'   independently: at 64 workers the CPU default of 2e9 is a 128 GB claim, in
#'   a stage that has already been OOM-killed once at 503 GB MaxRSS. Dividing
#'   here makes the documented budget the total it says it is. The GPU path
#'   takes the same division -- one device, several processes.
#' @return a single integer, genes per covariance sub-batch (at least 1).
#' @noRd
.polishCovBatchSize <- function(ncells, p, backend, gpu.mem.budget = NULL,
                          nworkers = 1L) {
  gpu_active <- backend %in% c("gpu", "auto") && checkGPU()
  budget <- if (gpu_active) {
    getGPUMemoryBudget(gpu.mem.budget)
  } else {
    # the spiDE option name is the fallback so a setting made before the move keeps working
    getOption("SpaNorm.cov.mem.budget",
              getOption("spiDE.cov.mem.budget", POLISH_COV_MEM_BUDGET_CPU))
  }
  if (!is.finite(budget)) {
    budget <- POLISH_COV_MEM_BUDGET_CPU
  }
  budget <- budget / max(1L, as.integer(nworkers))
  bytes <- if (gpu_active) gpuDtypeBytes() else 8
  ncells <- as.numeric(ncells)
  p <- as.numeric(p)

  # GPU: the (batch, ncells, p) weighted design plus its transpose view, then
  # the (batch, p, p) Gram/Cholesky/inverse stack. CPU: the (batch, p, p)
  # stack only -- construction there is one p x p crossprod at a time.
  per_gene <- if (gpu_active) {
    bytes * (ncells * p * 2 + p^2 * POLISH_TENSOR_MULT_COV)
  } else {
    bytes * p^2 * POLISH_TENSOR_MULT_COV
  }
  max(1L, as.integer(floor(budget * POLISH_BUDGET_FRACTION / per_gene)))
}
