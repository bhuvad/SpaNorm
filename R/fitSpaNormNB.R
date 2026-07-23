#' @importFrom stats mad median quantile model.matrix pnbinom dnbinom qnbinom
NULL

# Default winsorisation strength (number of MADs), shared by winsoriseCols(),
# calculateMu(), winsorisePsi(), the fit functions below, and subsetFitGenes()
# in mainSpaNorm.R, so the block-wise normalisation path stays in sync with the
# whole-matrix default from a single source of truth.
DEFAULT_WINSOR <- 4

#' Fit a per-gene negative binomial GLM
#'
#' Fits a per-gene negative binomial regression over an arbitrary design matrix
#' using SpaNorm's IRLS engine (iteratively reweighted least squares with
#' per-gene dispersion estimated by \code{edgeR::estimateDisp}). This exposes the
#' fitting machinery independently of SpaNorm's spatial model: the fitted model
#' is \code{log(mu) = W \%*\% t(alpha)} with no built-in intercept, so encode one
#' as a column of \code{W} when a per-gene baseline is required. It corresponds to
#' the internal fit with \code{is.spanorm = FALSE} (no shared library-size column,
#' ridge applied to every column of \code{W}).
#'
#' @param Y a genes x cells matrix of counts (dense, sparse, or DelayedArray).
#' @param W a cells x covariates numeric design matrix.
#' @param idx a logical vector (length \code{ncol(Y)}) selecting the cells used to
#'   fit the model (default: all cells).
#' @param lambda.a a numeric ridge penalty on the columns of \code{W}: a single
#'   value or a per-column vector (default 0, i.e. unregularised).
#' @param winsor a numeric, the number of MADs at which per-gene coefficients are
#'   winsorised during fitting (default 4). Must be a single positive number;
#'   \code{Inf} disables winsorisation entirely.
#' @param maxit.psi a numeric, the maximum number of dispersion iterations.
#' @param maxit.nb a numeric, the maximum number of NB mean IRLS iterations.
#' @param tol a numeric, the convergence tolerance.
#' @param ... additional fitting parameters forwarded to the internal fitter,
#'   e.g. \code{maxn.psi} (dispersion-estimation subsample size) or
#'   \code{step.factor} (IRLS step-halving factor).
#' @param backend a character, the compute backend ('auto', 'cpu', or 'gpu').
#' @param verbose a logical, whether to print progress messages (default TRUE).
#'
#' @return a list with per-gene coefficients \code{alpha} (genes x covariates),
#'   dispersions \code{psi}, a \code{gmean} element (always zero -- the generic
#'   fit has no intercept term), the \code{sampling} factor, and per-iteration
#'   \code{loglik}.
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rpois(20 * 50, 5), 20, 50)
#' W <- cbind(1, scale(seq_len(50))) # intercept + one covariate
#' fit <- fitNB(Y, W, verbose = FALSE)
#' str(fit$alpha)
#'
#' @export
fitNB <- function(Y, W, idx = rep(TRUE, ncol(Y)), lambda.a = 0, winsor = DEFAULT_WINSOR,
                  maxit.psi = 25, maxit.nb = 50, tol = 1e-4, ...,
                  backend = c("auto", "cpu", "gpu"), verbose = TRUE) {
  backend <- match.arg(backend)
  fitSpaNormNB(Y, W, idx, lambda.a = lambda.a, is.spanorm = FALSE, winsor = winsor,
               maxit.psi = maxit.psi, maxit.nb = maxit.nb, tol = tol, ...,
               backend = backend, msgfun = if (verbose) message else \(...) {})
}

fitSpaNormNB <- function(Y, W, idx, maxit.psi = 25, tol = 1e-4, maxn.psi = 500, ..., is.spanorm = FALSE, winsor = DEFAULT_WINSOR, backend = c("auto", "cpu", "gpu"), gpu.mem.budget = NULL, msgfun = message) {
  backend = match.arg(backend)
  # parameter checks
  if (maxit.psi <= 0) {
    stop("'maxit.psi' should be greater than 0")
  }
  if (maxit.psi - as.integer(maxit.psi) != 0) {
    stop("'maxit.psi' should be integers")
  }

  # message function for inner loop
  msgfunNB = \(x, ...) {
    \(...) msgfun(x, ...)
  }

  # subset data points used for model fitting. Ysub stays a plain, CPU-resident
  # matrix for the whole fit -- the fitting machinery always operates by
  # iterating contiguous gene-blocks (a single, whole-gene block when the
  # problem already fits the accelerator's memory budget, e.g. backend = "cpu"
  # or a small dataset), so no genes x cells tensor is ever fully materialised
  # on the accelerator at once.
  Ysub = as.matrix(Y[, idx, drop = FALSE])
  Wsub = toGPUMatrix(W[idx, , drop = FALSE], backend = backend) # Wsub is always small; never blocked
  nW = ncol(Wsub)
  nsub = sum(idx)
  ngenes = nrow(Ysub)

  # decide once (per fit) how many gene-blocks are needed to stay within the
  # accelerator's memory budget; geneBlockCount() returns 1L for
  # backend = "cpu", when no accelerator is available, or when the problem
  # already fits -- `blocks` is always a real (>= 1 block) partition.
  nblocks = geneBlockCount(ngenes, nsub, nW, backend = backend, budget.bytes = getGPUMemoryBudget(gpu.mem.budget = gpu.mem.budget))
  blocks = geneBlockIndices(ngenes, nblocks)

  # setting initial regression parameter estimates for all genes. alpha/gmean/
  # psi stay plain R objects for the whole fit (they are always small: genes x
  # ncov or length-genes) -- only each block's big intermediates are pushed to
  # the accelerator, inside the block-iteration helpers below. The generic
  # fit has no per-gene intercept, so gmean stays 0 (see fitNBGivenPsi).
  gmean = if (is.spanorm) rowMeans_gpu(log(Ysub + 1)) else rep(0, nrow(Ysub)) # rowMeans(Z)
  alpha = matrix(0, nrow(Ysub), ncol(W))
  # alpha for logLS
  alpha[, 1] = 1

  # subset to a maximum of 50 cells/spots for dispersion parameter estimation (for speed-up)
  nsub.psi = min(maxn.psi, nsub)
  psi.idx = rep(FALSE, nsub)
  psi.idx[sample.int(nsub, size = nsub.psi)] = TRUE
  psi.cols = which(psi.idx)
  # realise only the dispersion subset (genes x nsub.psi), not all sampled cells
  Ysub.psi = as.matrix(Ysub[, psi.cols, drop = FALSE])
  # small (nsub.psi x ncov) slice reused every outer iteration by the offs.psi
  # computation, so the full genes x nsub offset matrix is never built
  Wsub.psi = Wsub[psi.cols, , drop = FALSE]
  # a separate, usually coarser block partition sized for nsub.psi (capped at
  # maxn.psi) rather than reusing the full-fit `blocks` (sized for the much
  # larger nsub) -- offs.psi only ever touches nsub.psi cells, so it would
  # otherwise run far more block iterations than its own memory footprint needs
  nblocks.psi = geneBlockCount(ngenes, nsub.psi, nW, backend = backend, budget.bytes = getGPUMemoryBudget(gpu.mem.budget = gpu.mem.budget))
  blocks.psi = geneBlockIndices(ngenes, nblocks.psi)
  psi = rep(0, nrow(Ysub))

  # convergence trackers
  conv = FALSE
  logl.psi = c()
  iter = 1

  # map sample ids used for estimation for export
  idx = as.numeric(idx)
  idx[idx == 1][psi.idx] = 2
  idx = as.factor(c("all", "glm", "dispersion")[idx + 1])

  # outer loop of IRLS for estimating disp parameter
  while (!conv) {
    msgfun(sprintf("iter:%3d, estimating gene-wise dispersion", iter))
    # calculate/update dispersion parameter (psi) estimates for all genes,
    # computed directly against the small Wsub.psi slice (using the coarser
    # nsub.psi-sized block partition) so the full genes x nsub offset matrix
    # is never built
    offs.psi = .blockedOffsPsi(alpha, Wsub.psi, blocks.psi)
    psi.tmp = tryCatch(
      {
        edgeR::estimateDisp(Ysub.psi, as.matrix(rep(1, nsub.psi)), offset = offs.psi, tagwise = TRUE, robust = TRUE)$tagwise.dispersion
      },
      error = function(e) {
        rep(NA, length(psi))
      }
    )
    valid.psi = !is.na(psi.tmp) & !is.infinite(psi.tmp)
    psi[valid.psi] = psi.tmp[valid.psi]

    # calculate initial loglik for this iteration
    loglik = .blockedNBLoglik(alpha, gmean, Ysub, Wsub, psi, backend, blocks)
    msgfun(sprintf("iter:%3d, log-likelihood: %f", iter, sum(loglik)))

    # fit NB given dispersion estimates (psi) and extract required components
    msgfun(sprintf("iter:%3d, fitting NB model", iter))
    fit.nb = fitNBGivenPsi(Ysub, Wsub, psi, ..., gmean = gmean, alpha = alpha, loglik = loglik, is.spanorm = is.spanorm, winsor = winsor, backend = backend, blocks = blocks, msgfun = msgfunNB(sprintf("iter:%3d, ", iter)))
    gmean = fit.nb$gmean
    alpha = fit.nb$alpha
    loglik = fit.nb$loglik

    # check convergence
    logl.psi = c(sum(loglik), logl.psi) # push to stack
    conv.logl = FALSE
    iter = iter + 1 # similar post-increment to fitNBGivenPsi
    if (iter > 2) { # should be 2 because iter is post-incremented
      conv.logl = ifelse((logl.psi[1] - logl.psi[2]) / abs(logl.psi[2]) < tol, TRUE, FALSE)
    }
    conv = conv.logl | iter > maxit.psi
  }

  # convergence message
  if (iter > maxit.psi) {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (did not converged)", iter, logl.psi[1]))
  } else {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (converged)", iter, logl.psi[1]))
  }

  # final values
  fit.spanorm = list(
    gmean = as.vector(gmean),
    alpha = as.matrix(alpha),
    psi = as.vector(psi),
    sampling = idx,
    loglik = rev(logl.psi)
  )

  return(fit.spanorm)
}

fitNBGivenPsi <- function(Ysub, Wsub, psi, lambda.a, gmean = NULL, alpha = NULL, step.factor = 0.5, maxit.nb = 50, tol = 1e-4, loglik = NULL, backend = c("auto", "cpu", "gpu"), is.spanorm = FALSE, winsor = DEFAULT_WINSOR, blocks = NULL, msgfun = message) {
  backend = match.arg(backend)
  # parameter checks
  if (any(lambda.a < 0)) {
    stop("'lambda.a' should be positive")
  }
  # SpaNorm holds the shared logLS column (column 1) unpenalised, so it ridges
  # only the remaining columns; the generic fit penalises every column of W.
  npen = if (is.spanorm) ncol(Wsub) - 1 else ncol(Wsub)
  if (!length(lambda.a) %in% c(1, npen)) {
    stop("'lambda.a' should either be a single value or a vector of length equal to the number of penalised columns in W")
  }
  if (length(lambda.a) == 1) {
    lambda.a = rep(lambda.a, npen)
  }
  if (step.factor <= 0 | step.factor >= 1) {
    stop("'step.factor' should be in the interval (0,1)")
  }
  if (maxit.nb <= 0) {
    stop("'maxit.nb' should be greater than 0")
  }
  if (maxit.nb - as.integer(maxit.nb) != 0) {
    stop("'maxit.nb' should be integers")
  }

  # Ysub/alpha/gmean/psi stay plain R objects for the whole fit -- they are
  # always small (or, for Ysub, deliberately kept off the accelerator so no
  # genes x cells tensor is ever fully materialised) -- and only per-block
  # slices of the big genes x cells intermediates are pushed to the
  # accelerator, inside the block-iteration helpers below. Wsub is always
  # small and is converted once, up front.
  Ysub = as.matrix(Ysub)
  Wsub = toGPUMatrix(Wsub, backend = backend)
  lambda.a = diag_mat(lambda.a, backend = backend)

  if (is.null(gmean)) {
    # SpaNorm profiles a per-gene intercept; the generic fit has none (gmean = 0,
    # so an intercept must be encoded in W if wanted) -- see the gmean-fold update below
    gmean = if (is.spanorm) rowMeans_gpu(log(Ysub + 1)) else rep(0, nrow(Ysub)) # rowMeans(Z)
  }
  if (is.null(alpha)) {
    alpha = matrix(0, nrow(Ysub), ncol(Wsub))
    # alpha for logLS
    alpha[, 1] = 1
  }

  # decide (once) how many gene-blocks are needed, if the caller didn't
  # already supply a partition -- e.g. a direct call to fitNB()/
  # fitNBGivenPsi() bypassing fitSpaNormNB()'s own (reused-across-iterations)
  # decision; "helpful pattern for external use of the fitting function"
  if (is.null(blocks)) {
    nblocks = geneBlockCount(nrow(Ysub), ncol(Ysub), ncol(Wsub), backend = backend)
    blocks = geneBlockIndices(nrow(Ysub), nblocks)
  }

  # set the IRLS step size to 1
  nsub = ncol(Ysub)
  nW = ncol(Wsub)
  step = rep(1, nsub)

  if (is.null(loglik))  {
    # if not pre-computed (e.g., by outer loop), compute
    loglik = .blockedNBLoglik(alpha, gmean, Ysub, Wsub, psi, backend, blocks)
  }

  # convergence trackers
  conv = FALSE
  logl.beta = c(sum(loglik)) # stack of loglik values - push values
  iter = 1

  # initialise halving counters
  halving = 0

  while (!conv) {
    # message
    msgfun(sprintf("iter:%3d, log-likelihood: %f", iter, logl.beta[1]))

    # saving new best estimate (in case we need to go back in our search)
    best.gmean = gmean
    best.a = alpha

    iter.res = .fitNBGivenPsiIter(alpha, gmean, Ysub, Wsub, psi, step, lambda.a, is.spanorm, winsor, backend, blocks)
    alpha = iter.res$alpha
    gmean = iter.res$gmean
    loglik.tmp = iter.res$loglik

    # check degenerate case
    degener = sum(loglik) > sum(loglik.tmp)
    if (is.na(degener) | degener) {
      # alter stepsize of genes
      check.gene = loglik > loglik.tmp
      step[check.gene] = step[check.gene] * step.factor
      # revert to previous best
      gmean = best.gmean
      alpha = best.a
      # increase halving counter
      halving = halving + 1
      if (halving >= 3) {
        loglik.tmp = loglik
        degener = !degener
      }
    }

    # check degener (again in case halving exceeds maximum, i.e., 3)
    if (!degener) {
      loglik = loglik.tmp
      logl.beta = c(sum(loglik), logl.beta) # push to stack
      iter = iter + 1
      halving = 0 # reset
    }

    # check convergence
    conv.logl = FALSE
    if (iter > 2) { # should be 2 because iter is post-incremented
      conv.logl = ifelse((logl.beta[1] - logl.beta[2]) / abs(logl.beta[2]) < tol, TRUE, FALSE)
    }
    conv = conv.logl | iter > maxit.nb
  }

  # convergence message
  if (iter > maxit.nb) {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (did not converged)", iter, logl.beta[1]))
  } else {
    msgfun(sprintf("iter:%3d, log-likelihood: %f (converged)", iter, logl.beta[1]))
  }

  # final values
  fit.nb = list(
    gmean = as.numeric(gmean),
    alpha = as.matrix(alpha),
    loglik = loglik,
    loglik.iter = logl.beta
  )

  return(fit.nb)
}

# ---- Gene-blocked GPU fitting helpers -------------------------------------
# The IRLS fit always proceeds by iterating contiguous gene-blocks -- a
# single, whole-gene block when the problem already fits the accelerator's
# memory budget (including always, for backend = "cpu"), or several when
# geneBlockCount() (R/gpuFunctions.R) determines more are needed -- so no
# genes x cells tensor larger than one block is ever resident on the
# accelerator at once.
#
# alpha/gmean/psi/b stay plain R objects throughout (they are always small:
# genes x ncov or length-genes) -- only each block's big intermediates are
# pushed to the accelerator (via toGPUMatrix()/the *_gpu helpers' own
# auto-conversion of their other operand), then discarded (toRMatrix() +
# gc(FALSE)) before the next block starts. Wsub/Wsub.wt may be torch tensors
# (Wsub is small and never blocked) -- note torch's `[` interprets negative
# indices Python-style ("from the end"), NOT base R's "exclude this index",
# so column exclusion below always uses an explicit positive range
# (`2:ncov`), never `-1`.

# fitted log-mean for one gene-block: gmean[blk] + tcrossprod(alpha[blk,], Wsub)
.blockLmuHat <- function(alpha, gmean, blk, Wsub, backend) {
  add_vec_mat_gpu(gmean[blk], tcrossprod_gpu(alpha[blk, , drop = FALSE], Wsub), backend = backend)
}

# working response Z for one gene-block: lmu.hat.blk + step-scaled working
# residual, matching the non-blocked path's Zr/Z (fitNBGivenPsi() above)
.blockWorkingResponse <- function(Ysub.blk, lmu.hat.blk, step, backend) {
  Zr.blk <- (Ysub.blk + 0.01) / (exp(lmu.hat.blk) + 0.01) - 1
  lmu.hat.blk + mult_vec_mat_gpu(step, Zr.blk, backend = backend)
}

# blocked replacement for
# `as.numeric(colSums_gpu(dnbinom_gpu(Ysub, mu=exp(lmu.hat), size=1/psi, log=TRUE)))`.
# Reused for: the outer dispersion loop's initial loglik, the inner loop's
# initial loglik when not pre-supplied, and Pass C below.
.blockedNBLoglik <- function(alpha, gmean, Ysub, Wsub, psi, backend, blocks) {
  total <- NULL
  for (blk in blocks) {
    lmu.hat.blk <- .blockLmuHat(alpha, gmean, blk, Wsub, backend)
    Ysub.blk <- toGPUMatrix(Ysub[blk, , drop = FALSE], backend = backend)
    ll.blk <- as.numeric(toRMatrix(colSums_gpu(dnbinom_gpu(Ysub.blk, mu = exp(lmu.hat.blk), size = 1 / psi[blk], log = TRUE))))
    total <- if (is.null(total)) ll.blk else total + ll.blk
    gc(FALSE)
  }
  total
}

# blocked replacement for `tcrossprod_gpu(alpha, Wsub)[, psi.cols]` (the
# outer loop's dispersion-estimation offset): computed directly against the
# small Wsub.psi slice so the full genes x nsub offset matrix is never built.
.blockedOffsPsi <- function(alpha, Wsub.psi, blocks) {
  out <- matrix(0, nrow(alpha), nrow(Wsub.psi))
  for (blk in blocks) {
    out[blk, ] <- toRMatrix(tcrossprod_gpu(alpha[blk, , drop = FALSE], Wsub.psi))
    gc(FALSE)
  }
  out
}

# Pass B: blocked replacement for `matmul_gpu(Zc, Wsub.wt)`, the alpha-update
# numerator (`b`, genes x ncov). Uses this iteration's *starting* alpha/gmean
# (the same ones used in Pass A) -- recomputed rather than cached, since
# Pass B cannot start until Wsub.wt (which needs every block's Pass-A
# contribution) is known, so none of Pass A's big tensors survive to reuse.
.irlsBlockedAlphaNumerator <- function(alpha, gmean, Ysub, Wsub, step, Wsub.wt, backend, blocks) {
  b <- matrix(0, nrow(alpha), ncol(Wsub))
  for (blk in blocks) {
    lmu.hat.blk <- .blockLmuHat(alpha, gmean, blk, Wsub, backend)
    Ysub.blk <- toGPUMatrix(Ysub[blk, , drop = FALSE], backend = backend)
    Z.blk <- .blockWorkingResponse(Ysub.blk, lmu.hat.blk, step, backend)
    Zc.blk <- add_vec_mat_gpu(-gmean[blk], Z.blk, backend = backend)
    b[blk, ] <- toRMatrix(matmul_gpu(Zc.blk, Wsub.wt))
    gc(FALSE)
  }
  b
}

# One-shot (small-matrix, no blocking needed) solve of the new alpha from the
# fully-assembled `b`. Implements the Wa1-fold identity for is.spanorm: since
# Wsub.wt[, 2:ncov] == W_rest.wt already, `b`'s column 1 gives the
# unregularised a1_mean directly and `b`'s remaining columns give the
# regularised-solve numerator after subtracting a single gene-independent
# constant row -- so no separate Wa1/W_rest.wt block computation is needed,
# and (the hard constraint from the SpaNorm model, see fitNBGivenPsi() above)
# no gene's remaining-column alpha is finalised until a1_mean has been
# computed from every gene's contribution to `b`.
.irlsSolveAlphaFromB <- function(b, Wsub, Wsub.wt, lambda.a, is.spanorm) {
  if (is.spanorm) {
    Ainv <- toRMatrix(invert_mat(crossprod_gpu(Wsub.wt, Wsub)))
    alpha.unreg <- b %*% Ainv
    a1.mean <- mean(alpha.unreg[, 1])

    ncov <- ncol(Wsub)
    W.first <- Wsub[, 1, drop = FALSE]
    W.rest <- Wsub[, 2:ncov, drop = FALSE]
    W.rest.wt <- Wsub.wt[, 2:ncov, drop = FALSE]
    const.row <- as.numeric(toRMatrix(crossprod_gpu(W.first, W.rest.wt)))
    b.other <- sweep(b[, 2:ncov, drop = FALSE], 2, a1.mean * const.row, "-")
    Arest.inv <- toRMatrix(invert_mat(crossprod_gpu(W.rest.wt, W.rest) + lambda.a))
    alpha.other <- b.other %*% Arest.inv
    cbind(rep(a1.mean, nrow(b)), alpha.other)
  } else {
    Areg.inv <- toRMatrix(invert_mat(crossprod_gpu(Wsub.wt, Wsub) + lambda.a))
    b %*% Areg.inv
  }
}

# One-shot small-matrix gmean update (is.spanorm only) -- the gmean-fold
# identity: gmean = rowSums(Z.res*sig.inv)/rowSums(sig.inv) where
# Z.res = Z - tcrossprod(alpha.new, Wsub), decomposed into pieces already
# cached during Pass A (rowsum.Zsig, M = sig.inv %*% Wsub, rowsum.sig) so no
# genes x cells tensor needs to be rebuilt after alpha.new is known.
.irlsBlockedGmeanUpdate <- function(alpha.new, rowsum.Zsig, M, rowsum.sig) {
  (rowsum.Zsig - rowSums(alpha.new * M)) / rowsum.sig
}

# Pass A: reduce -- wt.cell (98th-percentile-clipped, matching the
# non-blocked path), plus (is.spanorm) the small per-gene quantities the
# gmean-fold needs later (rowsum.Zsig, M = sig.inv %*% Wsub, rowsum.sig).
# Ysub.blk is only realised on-device for the is.spanorm branch, since the
# generic (non-spanorm) fit never uses it here.
.irlsBlockedWtCell <- function(alpha, gmean, Ysub, Wsub, psi, step, is.spanorm, backend, blocks) {
  ngenes <- nrow(alpha)
  ncov <- ncol(Wsub)

  wtcell.num <- NULL
  rowsum.Zsig <- if (is.spanorm) numeric(ngenes) else NULL
  M <- if (is.spanorm) matrix(0, ngenes, ncov) else NULL
  rowsum.sig <- if (is.spanorm) numeric(ngenes) else NULL

  for (blk in blocks) {
    lmu.hat.blk <- .blockLmuHat(alpha, gmean, blk, Wsub, backend)
    sig.inv.blk <- 1 / add_vec_mat_gpu(psi[blk], exp(-lmu.hat.blk), backend = backend)

    colsum.sig.blk <- as.numeric(toRMatrix(colSums_gpu(sig.inv.blk)))
    wtcell.num <- if (is.null(wtcell.num)) colsum.sig.blk else wtcell.num + colsum.sig.blk

    if (is.spanorm) {
      Ysub.blk <- toGPUMatrix(Ysub[blk, , drop = FALSE], backend = backend)
      Z.blk <- .blockWorkingResponse(Ysub.blk, lmu.hat.blk, step, backend)
      rowsum.Zsig[blk] <- as.numeric(toRMatrix(rowSums_gpu(Z.blk * sig.inv.blk)))
      M[blk, ] <- toRMatrix(matmul_gpu(sig.inv.blk, Wsub))
      rowsum.sig[blk] <- as.numeric(toRMatrix(rowSums_gpu(sig.inv.blk)))
    }
    gc(FALSE)
  }

  wt.cell <- wtcell.num / ngenes
  wt.cell <- pmin(wt.cell, quantile(wt.cell, probs = 0.98))

  list(wt.cell = wt.cell, rowsum.Zsig = rowsum.Zsig, M = M, rowsum.sig = rowsum.sig)
}

# Orchestrates one full inner-IRLS iteration: Pass A (.irlsBlockedWtCell),
# the small Wsub.wt stage, Pass B (.irlsBlockedAlphaNumerator, solve
# numerator), the one-shot alpha/gmean finalisation, and Pass C
# (.blockedNBLoglik, loglik with the new alpha/gmean). Always used by
# fitNBGivenPsi()'s while-loop, regardless of block count.
.fitNBGivenPsiIter <- function(alpha, gmean, Ysub, Wsub, psi, step, lambda.a, is.spanorm, winsor, backend, blocks) {
  pass.a <- .irlsBlockedWtCell(alpha, gmean, Ysub, Wsub, psi, step, is.spanorm, backend, blocks)
  Wsub.wt <- mult_vec_mat_gpu(pass.a$wt.cell, Wsub, backend = backend)

  b <- .irlsBlockedAlphaNumerator(alpha, gmean, Ysub, Wsub, step, Wsub.wt, backend, blocks)

  # small stage: assemble and finalise new alpha (revert-on-bad-values +
  # winsorise, exactly as the non-blocked path does)
  alpha.new <- .irlsSolveAlphaFromB(b, Wsub, Wsub.wt, lambda.a, is.spanorm)
  if (hasBadValues(alpha.new)) alpha.new <- alpha
  alpha.new <- winsoriseCols(alpha.new, k = winsor)

  # small stage: gmean update (is.spanorm only; generic fit keeps gmean = 0)
  gmean.new <- if (is.spanorm) {
    .irlsBlockedGmeanUpdate(alpha.new, pass.a$rowsum.Zsig, pass.a$M, pass.a$rowsum.sig)
  } else {
    gmean
  }

  loglik.new <- .blockedNBLoglik(alpha.new, gmean.new, Ysub, Wsub, psi, backend, blocks)

  list(alpha = alpha.new, gmean = gmean.new, loglik = loglik.new)
}

normaliseLogPAC <- function(Y, scale.factor, fit.spanorm) {
  # calculate mu with all effects
  mu <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  # calculate mu with biology only
  isbio <- fit.spanorm$wtype %in% "biology"
  mu.2 <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], fit.spanorm$W[, isbio, drop = FALSE])
  # winsorize dispersion parameters (large disp slow the process)
  psi <- winsorisePsi(fit.spanorm$psi)

  # logPAC (densify counts once; a sparse/DelayedArray Y is realised transiently
  # here, or per gene-block when called via normaliseBlocked())
  Yd <- as.matrix(Y)
  lb <- pnbinom(Yd - 1, mu = mu, size = 1 / psi)
  ub <- dnbinom(Yd, mu = mu, size = 1 / psi) + lb
  p <- (lb + ub) / 2
  p <- pmax(pmin(p, 0.999), 0.001)

  # return logPAC
  normmat <- log2(qnbinom(p, mu = scale.factor * mu.2, size = 1 / psi) + 1)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
}

normaliseMeanBio <- function(Y, scale.factor, fit.spanorm) {
  isbio <- fit.spanorm$wtype %in% "biology"
  # calculate mu without library size effect
  normmat <- log2(calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], fit.spanorm$W[, isbio, drop = FALSE])) # log2(mu.2)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
}

normaliseMeanBatch <- function(Y, scale.factor, fit.spanorm) {
  isbio = fit.spanorm$wtype %in% "batch"
  # calculate mu without library size effect
  normmat = log2(calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], fit.spanorm$W[, isbio, drop = FALSE])) # log2(mu.2)
  colnames(normmat) = colnames(Y)
  rownames(normmat) = rownames(Y)

  return(normmat)
}

normaliseMeanLS <- function(Y, scale.factor, fit.spanorm) {
  isbio = fit.spanorm$wtype %in% "ls"

  # modify W to compute effect with the median library size
  W = fit.spanorm$W[, isbio, drop = FALSE]
  W = W / W[, 1] * median(W[, 1])

  # calculate mu without library size effect
  normmat = log2(calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], W)) # log2(mu.2)
  colnames(normmat) = colnames(Y)
  rownames(normmat) = rownames(Y)

  return(normmat)
}

normaliseMedianBio <- function(Y, scale.factor, fit.spanorm) {
  # calculate mu without library size effect
  isbio <- fit.spanorm$wtype %in% "biology"
  mu.2 <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, isbio, drop = FALSE], fit.spanorm$W[, isbio, drop = FALSE])
  # winsorize dispersion parameters (large disp slow the process)
  psi <- winsorisePsi(fit.spanorm$psi)

  # median Bio
  normmat <- log2(qnbinom_gpu(0.5, mu = scale.factor * mu.2, size = 1 / psi) + 1)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
}

normalisePearson <- function(Y, scale.factor, fit.spanorm) {
  isbio <- fit.spanorm$wtype %in% "biology"
  # calculate mu
  mu <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  # winsorize dispersion parameters (large disp slow the process)
  psi <- winsorisePsi(fit.spanorm$psi)

  # Pearson
  normmat <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha[, !isbio, drop = FALSE], fit.spanorm$W[, !isbio, drop = FALSE])
  normmat <- (Y - normmat) / sqrt(mu + mu^2 * psi)
  colnames(normmat) <- colnames(Y)
  rownames(normmat) <- rownames(Y)

  return(normmat)
}

# Validate a winsorisation strength (number of MADs): must be a single
# positive number. `Inf` is allowed and means "no winsorisation" -- callers
# must short-circuit on it themselves before computing `winsor * mad(...)`,
# which is NaN whenever the MAD is exactly 0 (e.g. several genes/columns
# sharing an identical fitted value).
.checkWinsor <- function(winsor) {
  if (!is.numeric(winsor) || length(winsor) != 1 || is.na(winsor) || winsor <= 0) {
    stop("'winsor' should be a single positive number")
  }
}

# Winsorise per-gene dispersions to exp(median(log psi) + winsor*MAD) (default
# 4), which bounds runtime for genes with very large dispersion. The threshold
# is a global statistic across all genes: when normalising in gene-blocks it is
# computed once on the full fit and carried as a 'psi.max' attribute on `psi`
# (set by subsetFitGenes) so block results are identical to the whole-matrix
# path; direct callers pass a plain vector and it is computed locally as before.
winsorisePsi <- function(psi, winsor = DEFAULT_WINSOR) {
  .checkWinsor(winsor)
  if (is.infinite(winsor)) return(as.numeric(psi))
  psi.max <- attr(psi, "psi.max")
  if (is.null(psi.max)) {
    psi.max <- exp(median(log(psi)) + winsor * mad(log(psi)))
  }
  pmin(as.numeric(psi), psi.max)
}

#' Compute fitted means from a negative binomial GLM fit
#'
#' Computes the fitted mean matrix \code{mu = exp(gmean + tcrossprod(alpha, W))}
#' from the per-gene coefficients of a negative binomial GLM, with optional
#' per-gene winsorisation of the log-means to \code{median +/- winsor * MAD} to
#' bound the influence of extreme fitted values. Exposed so downstream packages
#' (e.g. spiDE) can reconstruct fitted means from a \code{\link{fitNB}} result;
#' for the generic (intercept-free) fit pass \code{gmean = rep(0, nrow(alpha))}.
#'
#' @param gmean a numeric vector of per-gene intercepts (length \code{nrow(alpha)}).
#' @param alpha a genes x covariates matrix of coefficients.
#' @param W a cells x covariates design matrix.
#' @param winsor a numeric, the number of MADs at which per-gene log-means are
#'   winsorised (default 4); \code{Inf} disables winsorisation.
#' @return a genes x cells matrix of fitted means.
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rpois(20 * 50, 5), 20, 50)
#' W <- cbind(1, scale(seq_len(50)))
#' fit <- fitNB(Y, W, verbose = FALSE)
#' mu <- calculateMu(fit$gmean, fit$alpha, W)
#' dim(mu)
#'
#' @export
calculateMu <- function(gmean, alpha, W, winsor = DEFAULT_WINSOR) {
  .checkWinsor(winsor)
  # calculate mu (rather log of mu)
  # mu <- add_vec_mat_gpu(gmean, tcrossprod_gpu(alpha, W)) # log(mu)
  mu <- gmean + tcrossprod(alpha, W) # log(mu)
  if (is.infinite(winsor)) return(exp(mu))
  # winsorise to median +/- winsor*MAD (per gene)
  lmu.max <- matrixStats::rowMedians(mu) + winsor * matrixStats::rowMads(mu)
  mu <- exp(pmin(mu, lmu.max))

  return(mu)
}
