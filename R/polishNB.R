# polishNB(): the gene-blocked driver of the polish engine, moved from spiDE.
# Its body is spiDE's .polishFit() (R/polish.R); .reprofilePsi() comes from the
# same file. The gene-blocking and BLAS-thread helpers it dispatches through,
# .chunkGenes(), .singleBLAS(), .workerBLAS() and .bplapplySingleBLAS(), are
# copies of spiDE's (R/inference.R), which spiDE keeps for its inference stage.
# The per-gene and batched engines this drives are in R/polishEngine.R and
# R/polishEngineBatch.R; the Newton solvers in R/nbSolver.R.

#' Partition gene indices into blocks
#' @noRd
.chunkGenes <- function(n, block.size = NULL) {
  if (is.null(block.size) || block.size >= n) {
    return(list(seq_len(n)))
  }
  block.size <- max(1L, as.integer(block.size))
  split(seq_len(n), ceiling(seq_len(n) / block.size))
}

#' Should multi-worker dispatch run each worker's BLAS single-threaded?
#' @noRd
.singleBLAS <- function(BPPARAM) BiocParallel::bpnworkers(BPPARAM) > 1L

#' Run this worker's BLAS and OpenMP single-threaded
#'
#' Called inside a \code{bplapply()} worker of a multi-worker BPPARAM, never
#' in the parent.
#' @return the BLAS thread count before the call, so a persistent (snow)
#'   worker can be restored after the block.
#' @noRd
.workerBLAS <- function() {
  prev <- RhpcBLASctl::blas_get_num_procs()
  if (isTRUE(prev > 1L)) RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  invisible(prev)
}

#' bplapply() with single-threaded BLAS inside multi-worker dispatch
#'
#' A forked worker inherits the parent's OpenBLAS thread count, so N workers
#' on N cores run N x threads BLAS threads. Measured at the cohort's design
#' shape (77,454 cells, 345 dense + 660 nested columns; 2026-09-10): 4
#' workers x 4 threads on 4 cores take 2.3-2.6 s per Newton step against
#' 0.24-0.28 s at one thread each, while one worker gains only 1.5x from 4
#' threads -- the per-gene gram is memory-bound. The 0.99.19 cohort runs (4
#' workers x 8 threads on 8 cores) spent 245-395 min on the cold polish pass
#' this way. Both blocked stages -- the polish and the inference, whose
#' per-gene gram is the same BLAS work -- dispatch through this wrapper; the
#' niche construction and the gene-set stage are not gram-bound and do not.
#' A persistent worker gets its thread count back when the element is done;
#' a forked one dies with the call. Serial dispatch is plain
#' \code{bplapply()}.
#' @noRd
.bplapplySingleBLAS <- function(X, FUN, ..., BPPARAM = BiocParallel::SerialParam()) {
  if (!.singleBLAS(BPPARAM)) {
    return(BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM))
  }
  BiocParallel::bplapply(X, function(x, ...) {
    prev <- .workerBLAS()
    on.exit(if (isTRUE(prev > 1L)) RhpcBLASctl::blas_set_num_threads(prev), add = TRUE)
    FUN(x, ...)
  }, ..., BPPARAM = BPPARAM)
}

#' Apply a t tail to a statistic matrix, with per-column df
#'
#' \code{df} may be a scalar (one df for all columns),
#' a length-\code{ncol(tmat)} vector (a df per tested column, broadcast down the
#' rows), or a matrix matching \code{tmat}. \code{tmat} may be a genes x k matrix
#' or a length-k vector (a single gene). Centralises the reference-distribution

#' Converge every gene of a shared negative binomial fit to its own optimum
#'
#' \code{\link{fitNB}()} fits every gene in one IRLS loop: it shares a single
#' gene-averaged cell weight vector across genes, decides step-halving and
#' convergence on the aggregate log-likelihood, and clamps coefficients across
#' genes. For most genes the result is indistinguishable from the gene's own
#' optimum, but a bright, cell-type-restricted gene can be left well short of
#' it, with an inflated dispersion. \code{polishNB()} takes such a fit and
#' converges each gene separately, by damped Newton on that gene's own
#' penalised negative binomial log-likelihood
#' \deqn{\sum_i \log f_{NB}(y_i; \mu_i, \psi) - \frac{1}{2} \sum_j \lambda_j
#'   \alpha_j^2, \qquad \log \mu = W \alpha.}{sum_i log dnbinom(y_i; mu_i, psi)
#'   - 0.5 sum_j lambda_j alpha_j^2, with log(mu) = W alpha.}
#'
#' With \code{psi.method = "profile"} the dispersion is then re-estimated by
#' profile maximum likelihood at the converged mean and the mean re-polished,
#' twice; a gene whose dispersion optimum sits on the search bound keeps its
#' input dispersion (\code{psi_bound}). A gene whose starting point is
#' degenerate (a fitted log-mean below -10 at a cell with a positive count) or
#' whose Newton diverges restarts from a sane point: the log mean over the cells
#' of each \code{start.cols} column, or the overall log mean on the first
#' column. A gene that cannot be polished from either start keeps its input
#' coefficients and dispersion, with \code{polished = FALSE}.
#'
#' Once the shared fit is made the genes are independent, so this stage is
#' exact when blocked: the counts are densified one gene block at a time
#' (\code{block.size}) and the blocks are dispatched with \code{BPPARAM}. The
#' dispersion moderation in \code{fitNB()} works across genes, which is why
#' that fit must see the whole gene set and this one need not.
#'
#' @param Y a genes x cells matrix of integer counts (dense, sparse or
#'   DelayedArray; densified one gene block at a time). The negative binomial
#'   likelihood is undefined on non-integer values, so a non-integer assay is
#'   refused.
#' @param W a cells x p numeric design matrix.
#' @param alpha a genes x p matrix of starting coefficients, typically
#'   \code{fitNB()$alpha}.
#' @param psi the starting per-gene dispersions (length \code{nrow(alpha)}, or
#'   one value for every gene), typically \code{fitNB()$psi}.
#' @param lambda.a the ridge penalty, a single value or one per column of
#'   \code{W}. It is applied as given: each gene's objective subtracts
#'   \code{0.5 * sum(lambda.a * alpha^2)}, with no scaling by the number of
#'   cells or genes. The caller owns the scaling, so a fit made with a scaled
#'   penalty must pass the scaled values here.
#' @param absorb which columns of \code{W} the per-gene Newton solver absorbs
#'   by a Schur complement (exact; it makes a wide indicator block cost one
#'   dense-column gram per iteration). \code{NULL} (the default) absorbs
#'   nothing. A logical over the columns of \code{W} makes each marked column
#'   its own 1x1 block; the marked columns must be 0/1 indicators that
#'   partition the cells. A per-column grouping (integer, factor or character
#'   ids, \code{NA} for a dense column) makes the columns sharing an id one
#'   block; no cell may load on two blocks. See \code{\link{nbNewtonSolver}()}.
#' @param start.cols a logical over the columns of \code{W} marking the
#'   indicator columns (such as cell-type intercepts) that the sane start
#'   fills with the gene's log mean over that column's cells, or \code{NULL}
#'   to put the overall log mean on the first column.
#' @param psi.method how the dispersion is set at the converged mean:
#'   \code{"profile"} (profile maximum likelihood per gene) or \code{"fixed"}
#'   (the input \code{psi} is kept and only the mean is converged).
#' @param warm logical; \code{alpha} and \code{psi} are an already converged
#'   fit at a nearby penalty. A warm polish is a few damped Newton steps at the
#'   held dispersion, with no dispersion search and no restart check.
#' @param maxit,tol the Newton iteration cap and the relative log-likelihood
#'   tolerance, per gene.
#' @param engine \code{"batch"} (the default) runs a batch of genes through
#'   each Newton step so they share every read of the design; \code{"gene"} is
#'   the per-gene reference implementation. Both reach the same optimum; the
#'   batched profile dispersion is a bisection and the per-gene one
#'   \code{optimize()}, so they agree on \code{psi} to \code{optimize()}'s
#'   tolerance.
#' @param batch.size genes per batched Newton (\code{engine = "batch"}), or
#'   \code{NULL} to size it from a per-worker memory budget,
#'   \code{options(SpaNorm.polish.mem.budget = <bytes>)} (default 1e9; the
#'   older \code{spiDE.polish.mem.budget} is read as a fallback).
#' @param block.size genes per block, the unit of densification and dispatch,
#'   or \code{NULL} for at least one block per worker and at most 2,000 genes.
#' @param backend \code{"cpu"} (the default), \code{"auto"} or \code{"gpu"}. A
#'   device needs \code{engine = "batch"} and float64, and forces serial
#'   dispatch.
#' @param gpu.mem.budget accepted for symmetry with SpaNorm's other device
#'   entry points; not currently consulted (the batched working set is sized by
#'   \code{batch.size}).
#' @param BPPARAM a \code{BiocParallelParam} over gene blocks. With more than
#'   one worker, each worker runs its BLAS single-threaded (via RhpcBLASctl),
#'   because forked workers inherit the parent's thread count and
#'   oversubscribing the cores costs an order of magnitude per gene.
#' @param verbose logical; report progress.
#'
#' @return a list with \code{alpha} (genes x p), \code{psi} and \code{loglik}
#'   (the penalised log-likelihood at the returned fit, \code{NA} for a gene
#'   that was not polished), and \code{polish}, a data frame with one row per
#'   gene (row names from \code{alpha}): \code{iterations} (Newton steps),
#'   \code{psi_fitnb} (the input dispersion), \code{restarted} (the sane start
#'   was used), \code{capped} (a Newton pass hit \code{maxit}),
#'   \code{singular} (a singular information matrix), \code{psi_bound} (the
#'   dispersion optimum was on its search bound, so the input value was kept)
#'   and \code{polished} (\code{FALSE} when the gene kept its input fit).
#'
#' @seealso \code{\link{fitNB}()}, \code{\link{nbNewtonSolver}()}.
#' @examples
#' set.seed(1)
#' W <- cbind(1, rnorm(200))
#' Y <- t(replicate(5, rnbinom(200, mu = exp(1 + 0.3 * W[, 2]), size = 3)))
#' fit <- fitNB(Y, W, verbose = FALSE, backend = "cpu")
#' pol <- polishNB(Y, W, fit$alpha, fit$psi)
#' pol$polish
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @export
polishNB <- function(Y, W, alpha, psi, lambda.a = 0, absorb = NULL,
                     start.cols = NULL, psi.method = c("profile", "fixed"),
                     warm = FALSE, maxit = 50L, tol = 1e-8,
                     engine = c("batch", "gene"), batch.size = NULL,
                     block.size = NULL, backend = c("cpu", "auto", "gpu"),
                     gpu.mem.budget = NULL,
                     BPPARAM = BiocParallel::SerialParam(), verbose = FALSE) {
  pen <- lambda.a
  psi.method <- match.arg(psi.method)
  engine <- match.arg(engine)
  backend <- match.arg(backend)
  # short-circuits before checkGPU() probes, as .blockedInference() does
  gpu_active <- backend %in% c("gpu", "auto") && checkGPU()
  if (gpu_active) {
    .requireFloat64()
    if (engine != "batch") {
      stop("backend = \"", backend, "\" needs engine = \"batch\": the per-gene ",
           "engine keeps one factorisation object per gene, which is what ",
           "cannot go to a device.", call. = FALSE)
    }
    if (BiocParallel::bpnworkers(BPPARAM) > 1L) {
      warning("GPU backend active (", getBackendDevice(), "); forcing ",
              "serial dispatch for the polish, since forked workers contend ",
              "for one device", call. = FALSE)
      BPPARAM <- BiocParallel::SerialParam()
    }
  }
  ng <- nrow(alpha)
  if (!length(pen) %in% c(1L, ncol(W))) {
    stop("'lambda.a' must be a single value or one per column of W (",
         ncol(W), " here, ", length(pen), " supplied).", call. = FALSE)
  }
  if (length(pen) == 1L) pen <- rep(pen, ncol(W))
  if (length(psi) == 1L) psi <- rep(psi, ng)
  # The negative binomial likelihood is defined on counts. On a non-integer
  # assay dnbinom() returns -Inf for every cell, so the line search rejects
  # every step, alpha is left exactly at fitNB's value and the dispersion
  # optimiser -- maximising a constant -Inf -- returns its upper bound for EVERY
  # gene. Measured: psi 999.96 across the board, with `capped` and `singular`
  # both reporting success. That is silent, total corruption, and the documented
  # real-cohort object (counts <- 2^logcounts - 1) is exactly such an assay, so
  # refuse it here rather than let it through.
  chk <- as.numeric(Y[seq_len(min(nrow(Y), 20L)), , drop = FALSE])
  chk <- chk[is.finite(chk)]
  if (length(chk) && max(abs(chk - round(chk))) > 1e-8) {
    stop("the polish stage needs integer counts: the negative binomial ",
         "likelihood is undefined otherwise, and every gene's dispersion would ",
         "silently collapse to its upper bound.\n  The assay passed is not ",
         "integer-valued (e.g. a back-transform such as 2^logcounts - 1).\n  ",
         "Use the raw counts.", call. = FALSE)
  }
  # `absorb` is the per-gene solver's absorption: the nested indicators as a
  # logical, or the per-sample grouping of a random-slope fit's whole random
  # block. NULL absorbs nothing.
  nested <- if (!is.null(absorb)) absorb else rep(FALSE, ncol(W))
  solver <- .newtonSolver(W, pen, nested)
  # The shared-factor batched solver (.newtonSolverBatch(), built inside
  # .polishBatch() when a GPU is active) absorbs 1x1 blocks only. A logical
  # `absorb` is exactly that and is passed through unchanged. A grouping cannot
  # be reduced to its 1x1 blocks generically, so on that path it gets the dense
  # batched solver: exact, only slower.
  nested_batch <- if (is.logical(nested)) nested else rep(FALSE, ncol(W))
  # the design goes to the device once, outside the gene blocks
  W_use <- if (gpu_active) toGPUMatrix(W, backend = backend) else W

  # .chunkGenes(ng, NULL) is ONE block, which would hand every gene to a single
  # worker however many BPPARAM has -- a silent loss of the parallelism the
  # caller asked for. Absent an explicit block.size, split at least one block
  # per worker (this stage is exact per gene, so blocking never changes the
  # answer -- test-polish.R asserts that).
  # .chunkGenes(ng, NULL) is ONE block. That would (a) hand every gene to a
  # single worker however many BPPARAM has, and (b) -- worse -- densify the
  # WHOLE counts matrix at line `as.matrix(Y[gi, ])`: 13,348 x 77,454 doubles is
  # 8.3 GB, inside fitSpiDE(), under its own default SerialParam(). The
  # architecture's invariant is that the counts matrix is never eagerly
  # densified, so cap the block regardless of worker count, the way
  # .blockLoglik() does with its 2000-gene default.
  nw <- max(1L, BiocParallel::bpnworkers(BPPARAM))
  if (is.null(block.size)) {
    block.size <- max(1L, min(2000L, ceiling(ng / nw)))
  }
  blocks <- .chunkGenes(ng, block.size)
  single_blas <- .singleBLAS(BPPARAM)
  bsize <- if (engine == "batch") {
    if (is.null(batch.size)) .polishBatchSize(nrow(W)) else batch.size
  } else NA_integer_
  if (verbose) {
    message(sprintf("  %s %d genes %s (%d block%s%s)",
                    if (warm) "re-polishing" else "converging", ng,
                    if (engine == "batch") sprintf("in batches of %d", bsize) else "per gene",
                    length(blocks), if (length(blocks) == 1L) "" else "s",
                    if (single_blas) ", one BLAS thread per worker" else ""))
  }
  # A whole-transcriptome polish is hours of work, so report progress rather
  # than going silent after the opening message. Blocks are timed as they
  # finish; under a parallel BPPARAM they complete out of order, so the count
  # is of blocks retired, not a position in the gene list.
  t0 <- Sys.time()
  nb <- length(blocks)
  step <- max(1L, nb %/% 20L)
  res <- .bplapplySingleBLAS(seq_along(blocks), function(b) {
    gi <- blocks[[b]]
    # densify the whole block once: Y may be sparse or a DelayedArray, where a
    # per-gene read costs a round trip each time (the invariant is that the
    # WHOLE matrix is never densified, not that a block is never densified --
    # .blockedInference() does exactly the same).
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    out <- if (engine == "batch") {
      # a block bounds densification of the counts; the batched Newton's
      # working set is gene x cell and is bounded separately, so a block is
      # walked in sub-batches rather than handed over whole
      unlist(lapply(.chunkGenes(length(gi), bsize), function(ii) {
        Yblk <- Yb[ii, , drop = FALSE]
        Ablk <- alpha[gi[ii], , drop = FALSE]
        if (gpu_active) {
          Yblk <- toGPUMatrix(Yblk, backend = backend)
          Ablk <- toGPUMatrix(Ablk, backend = backend)
        }
        r <- .polishBatch(Yblk, W_use, Ablk,
                          psi[gi[ii]], pen, solver, maxit = maxit, tol = tol,
                          start.cols = start.cols, psi.method = psi.method, warm = warm,
                          shared.factor = gpu_active, nested = nested_batch)
        # back to the per-gene shape the merge below and @polish expect
        lapply(seq_along(ii), function(j) {
          list(alpha = r$alpha[j, ], psi = r$psi[[j]], loglik = r$loglik[[j]],
               iterations = r$iterations[[j]], restarted = r$restarted[[j]],
               capped = r$capped[[j]], singular = r$singular[[j]],
               psi_bound = r$psi_bound[[j]], polished = r$polished[[j]])
        })
      }), recursive = FALSE)
    } else {
      lapply(seq_along(gi), function(i) {
        g <- gi[[i]]
        .polishGene(as.numeric(Yb[i, ]), W, alpha[g, ], psi[[g]], pen, solver,
                    maxit = maxit, tol = tol, start.cols = start.cols,
                    psi.method = psi.method, warm = warm)
      })
    }
    if (verbose && (b %% step == 0L || b == nb)) {
      message(sprintf("    block %d/%d (%.1f min elapsed)", b, nb,
                      as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
    out
  }, BPPARAM = BPPARAM)
  res <- unlist(res, recursive = FALSE)
  if (verbose) {
    message(sprintf(
      "  polished %d/%d genes (%d restarted, %d not converged, %d singular, %d dispersion at a bound)",
      sum(vapply(res, `[[`, logical(1), "polished")), ng,
      sum(vapply(res, `[[`, logical(1), "restarted")),
      sum(vapply(res, `[[`, logical(1), "capped")),
      sum(vapply(res, `[[`, logical(1), "singular")),
      sum(vapply(res, `[[`, logical(1), "psi_bound"))))
  }

  out_alpha <- alpha
  out_psi <- psi
  for (g in seq_len(ng)) {
    out_alpha[g, ] <- res[[g]]$alpha
    out_psi[[g]] <- res[[g]]$psi
  }
  polish <- data.frame(
    iterations = vapply(res, `[[`, integer(1), "iterations"),
    psi_fitnb = psi,
    restarted = vapply(res, `[[`, logical(1), "restarted"),
    capped = vapply(res, `[[`, logical(1), "capped"),
    singular = vapply(res, `[[`, logical(1), "singular"),
    # the dispersion optimum sat on its search bound, so fitNB's moderated psi
    # was kept for this gene rather than a boundary value stored as an estimate
    psi_bound = vapply(res, `[[`, logical(1), "psi_bound"),
    # FALSE when the gene fell back to fitNB's fit entirely
    polished = vapply(res, `[[`, logical(1), "polished"),
    row.names = rownames(alpha)
  )
  list(alpha = out_alpha, psi = out_psi,
       loglik = vapply(res, `[[`, numeric(1), "loglik"), polish = polish)
}

#' Profile-ML dispersion at fixed coefficients, blocked over genes
#'
#' The warm passes of the variance-component loop hold each gene's dispersion,
#' so after the loop the stored psi is the profile optimum at a mean the stage
#' has moved on from. One profile pass at the final coefficients (an
#' optimiser call per gene, no Newton) puts it at the reported mean; a gene
#' whose optimum sits on a search bound keeps its value, as in the cold pass.
#' @noRd
.reprofilePsi <- function(Y, W, alpha, psi, psi.range = c(1e-3, 1e3),
                          block.size = NULL, BPPARAM = BiocParallel::SerialParam()) {
  ng <- nrow(alpha)
  nw <- max(1L, BiocParallel::bpnworkers(BPPARAM))
  if (is.null(block.size)) block.size <- max(1L, min(2000L, ceiling(ng / nw)))
  blocks <- .chunkGenes(ng, block.size)
  # One dispersion optimiser in the package, not two. This used to run its own
  # per-gene optimize(), which would have OVERWRITTEN the batched engine's more
  # accurate bisection at the last step of the default path -- discarding the
  # change for production while keeping it in the engine's internals. Same
  # kernel, same at_bound rule, and the rule's action is unchanged: a gene whose
  # dispersion runs to a bound keeps the one it came in with.
  res <- .bplapplySingleBLAS(blocks, function(gi) {
    Yb <- as.matrix(Y[gi, , drop = FALSE])
    Mu <- .muBatch(alpha[gi, , drop = FALSE], W)
    pm <- .psiProfileBatch(Yb, Mu, psi.range)
    out <- pm$psi
    keep <- pm$at_bound | !is.finite(out)
    out[keep] <- psi[gi][keep]
    out
  }, BPPARAM = BPPARAM)
  as.numeric(unlist(res))
}
