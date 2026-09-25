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

#' Is RhpcBLASctl installed?
#'
#' A package function rather than an inline \code{requireNamespace()}, so a
#' test can mock the probe without mocking base for the whole process.
#' @noRd
.hasBLASctl <- function() requireNamespace("RhpcBLASctl", quietly = TRUE)

#' Should multi-worker dispatch run each worker's BLAS single-threaded?
#'
#' Only with more than one worker, and only when RhpcBLASctl (in Suggests) is
#' installed to do it; without it the dispatch is plain \code{bplapply()} and
#' the workers keep the thread count they inherit.
#' @noRd
.singleBLAS <- function(BPPARAM) {
  BiocParallel::bpnworkers(BPPARAM) > 1L && .hasBLASctl()
}

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
#' a forked one dies with the call. Serial dispatch, and any dispatch without
#' RhpcBLASctl installed, is plain \code{bplapply()}.
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
#'   \alpha_j^2, \qquad \log \mu = W \alpha + \mathrm{offset}.}{sum_i log
#'   dnbinom(y_i; mu_i, psi) - 0.5 sum_j lambda_j alpha_j^2, with log(mu) =
#'   W alpha + offset.}
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
#' The integer-count check reads only the first 20 genes. It is a cheap guard
#' against an assay that is non-integer throughout (a back-transform such as
#' \code{2^logcounts - 1}, on which every gene's dispersion would silently
#' run to its upper bound), not a scan of every value: a non-integer count in
#' a later gene is not detected.
#'
#' @param Y a genes x cells matrix of integer counts (dense, sparse or
#'   DelayedArray; densified one gene block at a time). The negative binomial
#'   likelihood is undefined on non-integer values, so a non-integer assay is
#'   refused.
#' @param W a cells x p numeric design matrix.
#' @param alpha a genes x p matrix of starting coefficients, typically
#'   \code{fitNB()$alpha}.
#' @param psi the starting per-gene dispersions (length \code{nrow(Y)}, or
#'   one value for every gene), typically \code{fitNB()$psi}.
#' @param lambda.a the ridge penalty, a single value or one per column of
#'   \code{W}. It is applied as given: each gene's objective subtracts
#'   \code{0.5 * sum(lambda.a * alpha^2)}, with no scaling by the number of
#'   cells or genes. The caller owns the scaling, so a fit made with a scaled
#'   penalty must pass the scaled values here.
#' @param offset \code{NULL} (the default), a log-scale offset with one value
#'   per cell (the same for every gene), or a genes x cells matrix of them,
#'   added to every linear predictor with its coefficient fixed at 1:
#'   \code{log(mu) = W alpha + offset}. It must be finite. A \code{Matrix}
#'   offset is made a base matrix; a torch tensor is refused (pass it on the
#'   host; on a device the engine moves it). A matrix is held dense, so prefer
#'   a vector when the offset is the same for every gene.
#' @param absorb which columns of \code{W} the per-gene Newton solver absorbs
#'   by a Schur complement (exact; it makes a wide indicator block cost one
#'   dense-column gram per iteration). \code{NULL} (the default) absorbs
#'   nothing. A logical over the columns of \code{W} makes each marked column
#'   its own 1x1 block; the marked columns must be 0/1 indicators that
#'   partition the cells. A per-column grouping (integer, factor or character
#'   ids, \code{NA} for a dense column) makes the columns sharing an id one
#'   block; no cell may load on two blocks. See \code{\link{nbNewtonSolver}()}.
#' @param absorb.batch the absorption for the shared-factor batched solver,
#'   which is what runs on a device (\code{backend = "gpu"}, or \code{"auto"}
#'   with a GPU found): \code{NULL} or a logical over the columns of \code{W}.
#'   That solver absorbs 1x1 blocks only, so it cannot take a grouping with
#'   multi-column blocks. \code{NULL} (the default) passes a logical
#'   \code{absorb} through unchanged and sends a grouping to the dense
#'   batched solver (exact, only slower); a caller that knows the 1x1 subset
#'   of its grouping (e.g. the nested indicators inside a random-slope fit's
#'   per-sample blocks) passes it here. Not used on the CPU.
#' @param start.cols a logical over the columns of \code{W} marking the
#'   indicator columns (such as cell-type intercepts) that the sane start
#'   fills with the gene's log mean over that column's cells, or \code{NULL}
#'   to put the overall log mean on the first column.
#' @param psi.method how the dispersion is set at the converged mean:
#'   \code{"profile"} (profile maximum likelihood per gene) or \code{"fixed"}
#'   (the input \code{psi} is kept and only the mean is converged).
#' @param psi.range the search interval for the profile dispersion, two
#'   increasing positive numbers. A gene whose optimum falls on either end
#'   keeps its input dispersion.
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
#'   one worker and the RhpcBLASctl package installed, each worker runs its
#'   BLAS and OpenMP single-threaded, because forked workers inherit the
#'   parent's thread count and oversubscribing the cores costs an order of
#'   magnitude per gene. Without RhpcBLASctl the workers keep the thread count
#'   they inherit; install it, or set the BLAS threads to 1 (e.g.
#'   \code{OPENBLAS_NUM_THREADS=1}) before starting R, when using several
#'   workers.
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
polishNB <- function(Y, W, alpha, psi, lambda.a = 0, offset = NULL,
                     absorb = NULL, absorb.batch = NULL, start.cols = NULL,
                     psi.method = c("profile", "fixed"),
                     psi.range = c(1e-3, 1e3), warm = FALSE, maxit = 50L,
                     tol = 1e-8,
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
  .polishShapes(Y, W, alpha, psi)
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
  if (!is.null(absorb.batch) &&
      (!is.logical(absorb.batch) || length(absorb.batch) != ncol(W) ||
       anyNA(absorb.batch))) {
    stop("'absorb.batch' must be NULL or a logical with one value per column ",
         "of W (", ncol(W), " here)", call. = FALSE)
  }
  if (!is.numeric(psi.range) || length(psi.range) != 2L ||
      !all(is.finite(psi.range)) || psi.range[1] <= 0 ||
      psi.range[1] >= psi.range[2]) {
    stop("'psi.range' must be two increasing positive numbers", call. = FALSE)
  }
  # NULL, a per-cell vector or a genes x cells base matrix, finite, from here on
  offset <- .polishOffset(offset, Y)
  # `absorb` is the per-gene solver's absorption: the nested indicators as a
  # logical, or the per-sample grouping of a random-slope fit's whole random
  # block. NULL absorbs nothing.
  nested <- if (!is.null(absorb)) absorb else rep(FALSE, ncol(W))
  solver <- .newtonSolver(W, pen, nested)
  # The shared-factor batched solver (.newtonSolverBatch(), built inside
  # .polishBatch() when a GPU is active) absorbs 1x1 blocks only, so it gets
  # `absorb.batch` -- for a slope fit, the nested indicators inside the
  # per-sample blocks the per-gene solver above absorbs whole. Both are exact,
  # so this costs the device path a wider dense block and moves no result.
  # Without it a logical `absorb` is already 1x1 blocks and is passed through
  # unchanged; a grouping cannot be reduced to its 1x1 blocks generically, so
  # it gets the dense batched solver: exact, only slower.
  nested_batch <- if (!is.null(absorb.batch)) {
    absorb.batch
  } else if (is.logical(nested)) {
    nested
  } else {
    rep(FALSE, ncol(W))
  }
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
    # An offset adds genes x cells matrices to the batched working set: the
    # expansion .offsetRows() builds at every mean evaluation and, for a
    # per-gene matrix, the batch's own rows held for the whole batch. No offset
    # adds nothing, so the default count is used exactly.
    mats <- POLISH_GENE_CELL_MATS +
      if (is.null(offset)) 0 else 1 + is.matrix(offset)
    if (is.null(batch.size)) .polishBatchSize(nrow(W), mats = mats) else batch.size
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
        # a per-gene offset matrix is sliced to exactly these genes, from the
        # whole matrix (no block-sized copy is held); a per-cell vector is the
        # same for every gene. It stays on the host: the engine moves it.
        Oblk <- if (is.matrix(offset)) .rowsOf(offset, gi[ii]) else offset
        if (gpu_active) {
          Yblk <- toGPUMatrix(Yblk, backend = backend)
          Ablk <- toGPUMatrix(Ablk, backend = backend)
        }
        r <- .polishBatch(Yblk, W_use, Ablk,
                          psi[gi[ii]], pen, solver, maxit = maxit, tol = tol,
                          start.cols = start.cols, psi.method = psi.method, warm = warm,
                          shared.factor = gpu_active, nested = nested_batch,
                          psi.range = psi.range, offset = Oblk)
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
                    psi.method = psi.method, warm = warm, psi.range = psi.range,
                    offset = if (is.matrix(offset)) offset[g, ] else offset)
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

#' Check that the counts, design, coefficients and dispersions agree in shape
#'
#' A mismatch does not always fail by itself: a design with more rows than
#' \code{Y} has columns recycles each gene's counts against the longer mean on
#' the per-gene engine, which then reports the gene polished with wrong
#' coefficients. So \code{polishNB()} and \code{nbProfilePsi()} check up front
#' and name the argument and both sizes.
#' @param Y,W,alpha,psi as in \code{polishNB()}.
#' @return invisibly \code{TRUE}, or an error.
#' @noRd
.polishShapes <- function(Y, W, alpha, psi) {
  if (nrow(W) != ncol(Y)) {
    stop(sprintf("'W' must have one row per cell (column of Y): nrow(W) = %d, ncol(Y) = %d",
                 nrow(W), ncol(Y)), call. = FALSE)
  }
  if (nrow(alpha) != nrow(Y)) {
    stop(sprintf("'alpha' must have one row per gene (row of Y): nrow(alpha) = %d, nrow(Y) = %d",
                 nrow(alpha), nrow(Y)), call. = FALSE)
  }
  if (ncol(alpha) != ncol(W)) {
    stop(sprintf("'alpha' must have one column per column of W: ncol(alpha) = %d, ncol(W) = %d",
                 ncol(alpha), ncol(W)), call. = FALSE)
  }
  if (!length(psi) %in% c(1L, nrow(Y))) {
    stop(sprintf("'psi' must be one value or one per gene: length(psi) = %d, nrow(Y) = %d",
                 length(psi), nrow(Y)), call. = FALSE)
  }
  invisible(TRUE)
}

#' Check and normalise the offset for the polish
#'
#' One place decides what an offset may be, for \code{polishNB()} and
#' \code{nbProfilePsi()} alike: \code{NULL}, a finite numeric vector with one
#' value per cell (the same for every gene), or a finite numeric genes x cells
#' matrix. A \code{Matrix}-class offset is made a base matrix. A torch tensor
#' is refused: the engines take a host offset and move it themselves. A NaN
#' would otherwise stop both engines with "missing value where TRUE/FALSE
#' needed" and an Inf would come back as the input fit flagged singular.
#'
#' @param offset the offset as given.
#' @param Y the counts, for the expected shape.
#' @return \code{NULL}, a plain numeric vector or a base numeric matrix.
#' @noRd
.polishOffset <- function(offset, Y) {
  if (is.null(offset)) return(NULL)
  ng <- nrow(Y)
  nc <- ncol(Y)
  shape <- sprintf(paste0("NULL, a numeric vector with one value per cell ",
                          "(%d), or a numeric genes x cells matrix (%d x %d)"),
                   nc, ng, nc)
  if (is_torch_tensor(offset)) {
    stop("'offset' must be a base R vector or matrix, not a torch tensor: the ",
         "polish moves it to the device itself. Expected ", shape, ".",
         call. = FALSE)
  }
  if (methods::is(offset, "Matrix")) offset <- as.matrix(offset)
  ok <- if (is.null(dim(offset))) {
    is.numeric(offset) && length(offset) == nc
  } else {
    is.matrix(offset) && is.numeric(offset) && nrow(offset) == ng &&
      ncol(offset) == nc
  }
  if (!ok) stop("'offset' must be ", shape, ".", call. = FALSE)
  if (!all(is.finite(offset))) {
    stop("'offset' must be finite (no NA, NaN or Inf): expected ", shape, ".",
         call. = FALSE)
  }
  if (is.null(dim(offset))) as.numeric(offset) else offset
}

#' Profile-ML dispersion at fixed coefficients, blocked over genes
#'
#' For each gene, the negative binomial dispersion that maximises the
#' likelihood at the mean \code{exp(W alpha + offset)}, the coefficients held
#' fixed: one profile search per gene and no Newton, by the same bisection on
#' the score that \code{\link{polishNB}()}'s batched engine uses. A gene whose
#' optimum sits on a bound of \code{psi.range} (an under-dispersed or
#' near-empty gene), or whose search is not finite, keeps its input
#' dispersion rather than a boundary value stored as an estimate, the rule
#' \code{polishNB()} applies. A caller that re-polishes the mean at a held
#' dispersion (\code{polishNB(warm = TRUE)}) uses this to put the dispersion
#' back at the reported mean.
#'
#' @param Y a genes x cells matrix of counts (dense, sparse or DelayedArray;
#'   densified one gene block at a time).
#' @param W a cells x p numeric design matrix.
#' @param alpha a genes x p matrix of coefficients.
#' @param psi the per-gene dispersions (length \code{nrow(Y)}, or one value
#'   for every gene), kept for a gene whose optimum is on a bound.
#' @param psi.range the search interval for the dispersion.
#' @param block.size,BPPARAM gene blocking and dispatch, as in
#'   \code{polishNB()}.
#' @param offset \code{NULL}, a per-cell log-scale offset (length
#'   \code{ncol(Y)}) or a genes x cells matrix, added to every linear
#'   predictor. The dispersion is profiled at the mean it gives, so the offset
#'   must be the one the coefficients were fitted with.
#' @return a numeric vector of dispersions, one per gene.
#' @seealso \code{\link{polishNB}()}.
#' @examples
#' set.seed(1)
#' W <- cbind(1, rnorm(200))
#' Y <- t(replicate(3, rnbinom(200, mu = exp(1 + 0.3 * W[, 2]), size = 4)))
#' nbProfilePsi(Y, W, matrix(c(1, 0.3), 3, 2, byrow = TRUE), rep(0.5, 3))
#' @keywords internal
#' @export
nbProfilePsi <- function(Y, W, alpha, psi, psi.range = c(1e-3, 1e3),
                         block.size = NULL, BPPARAM = BiocParallel::SerialParam(),
                         offset = NULL) {
  .polishShapes(Y, W, alpha, psi)
  # a single value is every gene's, as in polishNB(); left scalar, `psi[gi]`
  # below would give NA to every at-bound gene after the first
  if (length(psi) == 1L) psi <- rep(psi, nrow(Y))
  offset <- .polishOffset(offset, Y)
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
    Ab <- alpha[gi, , drop = FALSE]
    # the offset rows of exactly these genes (NULL stays NULL, so the
    # no-offset mean is the pre-offset one bit for bit)
    Mu <- .muBatch(Ab, W, .offsetRows(offset, gi, Ab))
    pm <- .psiProfileBatch(Yb, Mu, psi.range)
    out <- pm$psi
    keep <- pm$at_bound | !is.finite(out)
    out[keep] <- psi[gi][keep]
    out
  }, BPPARAM = BPPARAM)
  as.numeric(unlist(res))
}
