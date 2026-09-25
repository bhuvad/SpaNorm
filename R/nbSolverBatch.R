# Batched / accelerator-aware helpers for the penalised NB Newton solve,
# moved from spiDE (R/inference-batch.R): per-gene Gram matrices, their
# batched Schur absorption and Cholesky solve, and the backend-agnostic array
# operations the batched Newton loop needs.
#
# Everything here is shape-and-memory plumbing that must behave identically on
# a base R matrix and a torch tensor -- the two backends differ in ways that
# are easy to get wrong, which is why they are collected in one place rather
# than inlined at their call sites.

#' Select rows of a matrix or a torch tensor
#'
#' Row-subsetting helper for the gene sub-batching in
#' \code{.waldCauchyBlock()}: base R \code{[} does not carry over to torch
#' tensors (which need \code{torch_index_select()}), and the working weights
#' may be either depending on the backend.
#'
#' @param x a matrix or torch tensor.
#' @param ii integer row positions to keep.
#' @return \code{x} restricted to rows \code{ii}, same type as the input.
#' @noRd
.rowsOf <- function(x, ii) {
  if (is_torch_tensor(x)) {
    idx <- torch::torch_tensor(as.integer(ii), dtype = torch::torch_long(),
                               device = x$device)
    return(torch::torch_index_select(x, 1, idx))
  }
  x[ii, , drop = FALSE]
}

#' Batched per-gene weighted Gram matrix, for a sub-batch of genes at once
#'
#' Computes the weighted information matrix \code{crossprod(W * wt_g, W)} =
#' \code{W' diag(wt_g) W} for every gene \code{g} in a sub-batch, returning a
#' \code{(batch, p, p)} array/tensor and optionally adding a shared ridge
#' penalty (mixed-effects) to every slice.
#'
#' On the accelerator this is a single batched matmul over the
#' \code{(batch, ncells, p)} weighted design; on CPU it is a per-gene
#' \code{crossprod()} loop (base R has no batched matmul, and the
#' construction is the cheap part relative to the inversion anyway).
#'
#' Peak memory here is \strong{linear in \code{p}}
#' (\code{batch x ncells x p}) and directly controllable via the sub-batch
#' size. An earlier implementation instead precomputed a Khatri-Rao /
#' face-splitting cross term of \code{W} (\code{ncells x p^2}, built once per
#' bandwidth and shared across all gene blocks), turning every gene's Gram
#' matrix into one row of a single large matmul. That is elegant, but it
#' scales \emph{quadratically} in \code{p} and -- being built once, outside
#' the gene loop -- could not be bounded by any choice of block size. On a
#' realistic mixed-effects design it is fatal: a 602-column random-intercept
#' design over 21,843 cells needs 63 GB, and a 4,906-column random-slope
#' design needs 4.2 TB, which made \code{combine = "cauchy"} unusable with
#' random effects on \emph{both} backends. Do not reintroduce that form
#' without bounding it by the gene sub-batch.
#'
#' @param W the design the Gram matrix is over (cells x p; a matrix or a
#'   torch tensor) -- \code{Wsub} for a fixed-effects fit, the full \code{W}
#'   for a mixed fit.
#' @param wt_block the sub-batch's working weights, \code{batch x ncells}.
#' @param penalty_diag if supplied, a length-\code{p} ridge penalty added to
#'   every slice's diagonal (mixed-effects only).
#' @param backend the resolved backend (unused on the base-R path; kept for
#'   signature symmetry with the other batched helpers).
#' @param cell.tile cells per accumulation tile, or \code{NULL} for all of them
#'   at once. The torch branch's weighted design is \code{batch x ncells x p} --
#'   ~8 GB at a 64-gene batch on the cohort's design, and NOT bounded by the
#'   gene sub-batch that is supposed to bound this stage, since it grows with
#'   \code{ncells}. Accumulating \code{W' diag(w) W} over cell tiles bounds the
#'   peak by the tile plus the \code{(batch, p, p)} stack instead. A Gram matrix
#'   is a sum over cells, so this is exact: every tiling returns the same stack,
#'   which \code{test-inference.R} pins on both branches.
#' @return a \code{(batch, p, p)} array (or torch tensor).
#' @noRd
.gramBatch <- function(W, wt_block, penalty_diag = NULL, backend = "cpu",
                       cell.tile = NULL) {
  is_t <- is_torch_tensor(W)
  n <- if (is_t) W$size(1) else nrow(W)
  p <- if (is_t) W$size(2) else ncol(W)
  b <- if (is_torch_tensor(wt_block)) wt_block$size(1) else nrow(wt_block)
  tile <- if (is.null(cell.tile)) n else max(1L, min(as.integer(cell.tile), n))
  starts <- seq(1L, n, by = tile)

  if (is_t) {
    info <- NULL
    for (s in starts) {
      ii <- s:min(s + tile - 1L, n)
      idx <- torch::torch_tensor(as.integer(ii), dtype = torch::torch_long(),
                                 device = W$device)
      Wi <- torch::torch_index_select(W, 1, idx)
      wi <- torch::torch_index_select(wt_block, 2, idx)
      # Weighting both factors by sqrt(wt) (rather than one by wt) keeps the
      # product exactly symmetric, which linalg_cholesky() in
      # invert_mat_batched() relies on; the working weights 1/(1/mu + psi) are
      # strictly positive, so the sqrt is safe.
      Wg <- torch::torch_sqrt(wi)$unsqueeze(3) * Wi$unsqueeze(1)
      part <- torch::torch_matmul(Wg$transpose(2, 3), Wg)
      info <- if (is.null(info)) part else info + part
    }
    if (!is.null(penalty_diag)) {
      pen <- torch::torch_tensor(as.numeric(penalty_diag),
                                 dtype = info$dtype, device = info$device)
      info <- info + torch::torch_diag(pen)$unsqueeze(1)
    }
    return(info)
  }

  pen_mat <- if (is.null(penalty_diag)) NULL else diag(penalty_diag, nrow = p)
  info <- array(0, dim = c(b, p, p))
  for (g in seq_len(b)) {
    ig <- matrix(0, p, p)
    for (s in starts) {
      ii <- s:min(s + tile - 1L, n)
      Wi <- W[ii, , drop = FALSE]
      ig <- ig + crossprod(Wi * wt_block[g, ii], Wi)
    }
    if (!is.null(pen_mat)) ig <- ig + pen_mat
    info[g, , ] <- ig
  }
  info
}

#' Sum the rows of a matrix within groups, on either backend
#'
#' \code{rowsum()} with no tensor equivalent is the stated reason the nested
#' indicator block is absorbed on the CPU only (\code{.blockedInference()}), so
#' the GPU path falls back to a dense \code{p x p} gram -- on the cohort's
#' design, 1,107 columns where 398 would do. That reason does not hold: torch
#' has \code{torch_index_add()}, which is exactly this operation, and
#' \code{torch_scatter_add()} and \code{torch_segment_reduce()} besides
#' (checked 2026-09-15).
#'
#' @param M a matrix or tensor, \code{n x k}.
#' @param gidx the group of each row, integers in \code{1:G}.
#' @param G the number of groups; groups absent from \code{gidx} come back as
#'   zero rows, which \code{rowsum()} would drop.
#' @return a \code{G x k} matrix or tensor of within-group column sums.
#' @noRd
.segmentSum <- function(M, gidx, G) {
  if (is_torch_tensor(M)) {
    idx <- torch::torch_tensor(as.integer(gidx), dtype = torch::torch_long(),
                               device = M$device)
    out <- torch::torch_zeros(c(G, M$size(2)), dtype = M$dtype, device = M$device)
    return(out$index_add(1, idx, M))
  }
  # levels = seq_len(G) so an absent group is a zero row rather than a missing
  # one: the caller indexes the result positionally
  gf <- factor(gidx, levels = seq_len(G))
  out <- matrix(0, G, ncol(M))
  r <- rowsum(M, group = gf, reorder = TRUE)
  out[as.integer(rownames(r)), ] <- r
  out
}

#' Batched Schur absorption of a nested indicator block, on either backend
#'
#' The nested (sample x cell type) columns are a 0/1 partition of the cells, so
#' \code{C = Z' diag(w) Z} is diagonal and the block can be absorbed exactly:
#' the covariance restricted to the dense columns is \code{S^-1}, with
#' \code{S = A - B C^-1 B'}. \code{.newtonSolver()}'s \code{parts()} does this
#' one gene at a time in base R and is the oracle this is tested against
#' (\code{test-absorb-batch.R}).
#'
#' Two things this adds. It is \strong{batched}, so the GPU inference path can
#' absorb instead of falling back to a dense \code{p x p} gram -- 1,107 columns
#' where 398 would do on the cohort's design, 7.7x the flops. And it is
#' \strong{tiled over cells}: the weighted design that feeds both \code{A} and
#' \code{B} is \code{batch x ncells x px}, ~8 GB at a 64-gene batch on that
#' design, so it is accumulated a cell-tile at a time and the peak is the
#' \code{(batch, px, px)} stack plus one tile, whatever \code{ncells} is.
#'
#' @param W the full design (cells x p), a matrix or a torch tensor.
#' @param pen the length-\code{p} ridge penalty.
#' @param nested logical, length \code{p}: the nested indicator columns.
#' @param wt_block the batch's working weights, \code{batch x ncells}, matching
#'   \code{W}'s type.
#' @param cell.tile cells per accumulation tile; \code{NULL} is all of them.
#'   Performance only -- every tiling returns the same stack, which
#'   \code{test-absorb-batch.R} pins.
#' @param parts if \code{TRUE}, return the pieces the back-substitution needs
#'   (\code{S}, \code{B}, \code{cvec}, \code{xi}, \code{zi}) rather than
#'   \code{S} alone. The inference path only ever wants \code{S}, since
#'   \code{S^-1} IS the covariance on the dense columns; the Newton step also
#'   needs \code{B} and \code{cvec} to recover the nested coefficients.
#'   \code{B} is carried as \code{(batch, G, px)}, the transpose of
#'   \code{.newtonSolver()}'s \code{px x G}, because that is the orientation
#'   both matmuls want.
#' @return a \code{(batch, px, px)} array or tensor of Schur complements, with
#'   \code{px = sum(!nested)}; or, with \code{parts = TRUE}, a list.
#' @noRd
.absorbBatch <- function(W, pen, nested, wt_block, cell.tile = NULL,
                         parts = FALSE) {
  is_t <- is_torch_tensor(W)
  colsOf <- function(M, ii) {
    if (is_torch_tensor(M)) {
      idx <- torch::torch_tensor(as.integer(ii), dtype = torch::torch_long(),
                                 device = M$device)
      return(torch::torch_index_select(M, 2, idx))
    }
    M[, ii, drop = FALSE]
  }
  tr2 <- function(M) if (is_torch_tensor(M)) M$transpose(1, 2) else t(M)

  xi <- which(!nested)
  zi <- which(nested)
  px <- length(xi)
  G <- length(zi)
  n <- if (is_t) W$size(1) else nrow(W)
  b <- if (is_torch_tensor(wt_block)) wt_block$size(1) else nrow(wt_block)
  pen_x <- pen[xi]
  pen_z <- pen[zi]

  X <- colsOf(W, xi)
  Zblk <- colsOf(W, zi)
  # the absorption is exact only if every cell belongs to exactly one group
  rs <- if (is_t) as.numeric(torch::torch_sum(Zblk, dim = 2)) else rowSums(Zblk)
  if (anyNA(rs) || max(abs(rs - 1)) > 1e-8) {
    stop("the nested random-effect columns are not 0/1 indicators partitioning ",
         "the cells; .absorbBatch() cannot absorb them", call. = FALSE)
  }
  # round(), not as.integer(): a floating-point product of 7 can come back as
  # 6.9999999, which as.integer() truncates to the WRONG group, silently
  gidx <- if (is_t) {
    sq <- torch::torch_tensor(as.numeric(seq_len(G)), dtype = Zblk$dtype,
                              device = Zblk$device)
    round(as.numeric(torch::torch_matmul(Zblk, sq)))
  } else {
    round(as.numeric(Zblk %*% seq_len(G)))
  }

  tile <- if (is.null(cell.tile)) n else max(1L, min(as.integer(cell.tile), n))
  A <- if (is_t) torch::torch_zeros(c(b, px, px), dtype = X$dtype, device = X$device)
       else array(0, c(b, px, px))
  Bacc <- if (is_t) torch::torch_zeros(c(b, G, px), dtype = X$dtype, device = X$device)
          else array(0, c(b, G, px))
  cvec <- if (is_t) torch::torch_zeros(c(b, G), dtype = X$dtype, device = X$device)
          else matrix(0, b, G)

  for (s in seq(1L, n, by = tile)) {
    ii <- s:min(s + tile - 1L, n)
    Xt <- .rowsOf(X, ii)
    wt <- colsOf(wt_block, ii)                       # batch x tile
    gt <- gidx[ii]
    A <- A + .gramBatch(Xt, wt)
    cvec <- cvec + tr2(.segmentSum(tr2(wt), gt, G))  # batch x G
    if (is_t) {
      # (batch, tile, px): bounded by the tile, which is the whole point
      M <- wt$unsqueeze(3) * Xt$unsqueeze(1)
      idx <- torch::torch_tensor(as.integer(gt), dtype = torch::torch_long(),
                                 device = M$device)
      Bacc <- Bacc$index_add(2, idx, M)
    } else {
      for (g in seq_len(b)) {
        Bacc[g, , ] <- Bacc[g, , ] + .segmentSum(Xt * wt[g, ], gt, G)
      }
    }
  }

  if (is_t) {
    A <- A + torch::torch_diag(torch::torch_tensor(as.numeric(pen_x),
                                                   dtype = A$dtype,
                                                   device = A$device))$unsqueeze(1)
    cvec <- cvec + torch::torch_tensor(as.numeric(pen_z), dtype = cvec$dtype,
                                       device = cvec$device)$unsqueeze(1)
    S <- A - torch::torch_matmul(Bacc$transpose(2, 3), Bacc / cvec$unsqueeze(3))
    if (!parts) return(S)
    return(list(S = S, B = Bacc, cvec = cvec, xi = xi, zi = zi))
  }
  pen_mat <- diag(pen_x, nrow = px)
  cvec <- cvec + rep(pen_z, each = b)
  S <- array(0, c(b, px, px))
  for (g in seq_len(b)) {
    Bg <- matrix(Bacc[g, , ], nrow = G, ncol = px)   # G x px
    S[g, , ] <- A[g, , ] + pen_mat - crossprod(Bg, Bg / cvec[g, ])
  }
  if (!parts) return(S)
  list(S = S, B = Bacc, cvec = cvec, xi = xi, zi = zi)
}

#' Per-slice Cholesky of a (batch, p, p) stack, with a per-slice verdict
#'
#' The trap this exists to avoid is the one \code{.waldCauchyBlock()} fell
#' into: a batched Cholesky over a stack containing one singular slice fails
#' the WHOLE stack, and the only recourse offered was telling the user to
#' shrink the batch. A per-slice \code{ok} lets a singular gene drop to its
#' fallback exactly as the per-gene engine does, leaving its neighbours alone.
#'
#' The Schur complement of a penalised information matrix is symmetric positive
#' definite where the design has full rank, so Cholesky is the right
#' factorisation and its failure IS the singularity test. \code{torch} reports
#' it per slice through \code{linalg_cholesky_ex()}'s \code{info} without
#' raising; base R needs a \code{tryCatch} per slice.
#'
#' @param S a \code{(batch, p, p)} array or tensor.
#' @return \code{list(L, ok)}; \code{L} is upper-triangular \code{R} with
#'   \code{S = R'R} on the base-R branch and lower-triangular \code{L} with
#'   \code{S = LL'} on the torch branch, matching each backend's own
#'   convention, and the solvers below respect that.
#' @noRd
.cholBatch <- function(S) {
  if (is_torch_tensor(S)) {
    r <- torch::linalg_cholesky_ex(S)
    info <- as.numeric(toRMatrix(r[[2]]))
    return(list(L = r[[1]], ok = info == 0))
  }
  b <- dim(S)[1]
  p <- dim(S)[2]
  L <- array(0, c(b, p, p))
  ok <- logical(b)
  for (g in seq_len(b)) {
    cg <- tryCatch(chol(matrix(S[g, , ], p, p)), error = function(e) NULL)
    if (!is.null(cg) && all(is.finite(cg))) {
      L[g, , ] <- cg
      ok[g] <- TRUE
    }
  }
  list(L = L, ok = ok)
}

#' Solve a batch of Cholesky-factorised systems, NA where the slice is singular
#'
#' @param ch the \code{list(L, ok)} from \code{.cholBatch()}.
#' @param rhs a \code{(batch, p)} matrix or tensor of right-hand sides.
#' @return a \code{(batch, p)} matrix or tensor; rows for \code{!ok} slices
#'   are \code{NA} (base R) or \code{NaN} (torch), which the caller reads as
#'   "this gene falls back".
#' @noRd
.cholSolveBatch <- function(ch, rhs) {
  if (is_torch_tensor(ch$L)) {
    out <- torch::torch_cholesky_solve(rhs$unsqueeze(3), ch$L)$squeeze(3)
    if (!all(ch$ok)) {
      bad <- which(!ch$ok)
      idx <- torch::torch_tensor(as.integer(bad), dtype = torch::torch_long(),
                                 device = out$device)
      out <- out$index_fill(1, idx, NaN)
    }
    return(out)
  }
  b <- dim(ch$L)[1]
  p <- dim(ch$L)[2]
  out <- matrix(NA_real_, b, p)
  for (g in seq_len(b)) {
    if (!ch$ok[g]) next
    R <- matrix(ch$L[g, , ], p, p)
    out[g, ] <- backsolve(R, backsolve(R, rhs[g, ], transpose = TRUE))
  }
  out
}

#' The block operations newton() needs, on either backend
#'
#' \code{.polishBatch()}'s Newton is written once and runs on a matrix or a
#' tensor, so the handful of array operations it does on \code{genes x cells}
#' and \code{genes x columns} blocks get one signature each. Every base-R
#' branch is plain indexing: converting the loop to these must not move the CPU
#' path at all, which test-polish-batch.R's strict parity tests are what prove.
#'
#' The division of labour they encode: the \strong{arrays} go to the device,
#' the \strong{bookkeeping} stays on the host. Which genes are active, which
#' accepted their step, how many halvings each has taken -- those are
#' length-\code{genes} integers and logicals, and pushing them to a device buys
#' nothing and costs a synchronisation per branch. \code{.asHost()} is the one
#' deliberate transfer per line-search round: the accept test needs the
#' log-likelihoods as R numbers.
#'
#' @name batch-array-ops
#' @noRd
NULL

#' @rdname batch-array-ops
#' @noRd
.setRows <- function(X, ii, V) {
  if (is_torch_tensor(X)) {
    idx <- torch::torch_tensor(as.integer(ii), dtype = torch::torch_long(),
                               device = X$device)
    return(X$index_copy(1, idx, V))
  }
  X[ii, ] <- V
  X
}

#' @rdname batch-array-ops
#' @noRd
.mulRows <- function(v, M) {
  if (is_torch_tensor(M)) {
    vt <- if (is_torch_tensor(v)) v else
      torch::torch_tensor(as.numeric(v), dtype = M$dtype, device = M$device)
    return(vt$unsqueeze(2) * M)
  }
  v * M
}

#' @rdname batch-array-ops
#' @noRd
.scaleCols <- function(A, s) {
  if (is_torch_tensor(A)) {
    st <- if (is_torch_tensor(s)) s else
      torch::torch_tensor(as.numeric(s), dtype = A$dtype, device = A$device)
    return(A * st$unsqueeze(1))
  }
  sweep(A, 2L, s, `*`)
}

#' @rdname batch-array-ops
#' @noRd
.asHost <- function(x) {
  if (is_torch_tensor(x)) return(as.numeric(toRMatrix(x)))
  x
}

#' @rdname batch-array-ops
#' @noRd
.matmulB <- function(X, Y) {
  if (is_torch_tensor(X)) return(torch::torch_matmul(X, Y))
  X %*% Y
}

#' @rdname batch-array-ops
#' @noRd
.rowsFinite <- function(X) {
  if (is_torch_tensor(X)) {
    ok <- torch::torch_isfinite(X)$all(dim = 2)
    return(as.logical(as.numeric(toRMatrix(ok))))
  }
  apply(X, 1L, function(r) all(is.finite(r)))
}

#' @rdname batch-array-ops
#' @noRd
.rowMeansB <- function(X) {
  if (is_torch_tensor(X)) {
    return(as.numeric(toRMatrix(torch::torch_mean(X, dim = 2))))
  }
  rowMeans(X)
}

#' @rdname batch-array-ops
#' @noRd
.colsOf <- function(X, jj) {
  if (is_torch_tensor(X)) {
    idx <- torch::torch_tensor(as.integer(jj), dtype = torch::torch_long(),
                               device = X$device)
    return(torch::torch_index_select(X, 2, idx))
  }
  X[, jj, drop = FALSE]
}

#' The smallest value in each row, over the masked entries only
#'
#' \code{degenerate()} asks for the smallest linear predictor among the cells
#' with a POSITIVE count. A gene with no positive count has no such cell, and
#' the answer there is \code{Inf} -- which is what makes it NOT degenerate,
#' since the test is \code{min < -10}. Returning an empty minimum (base R's
#' \code{min(numeric(0))} warns and returns \code{Inf}) or \code{NA} would
#' both be wrong in ways that only show on an all-zero gene.
#'
#' @rdname batch-array-ops
#' @noRd
.maskedRowMin <- function(X, mask) {
  if (is_torch_tensor(X)) {
    big <- torch::torch_full_like(X, Inf)
    return(as.numeric(toRMatrix(
      torch::torch_where(mask, X, big)$amin(dim = 2))))
  }
  vapply(seq_len(nrow(X)), function(i) {
    p <- mask[i, ]
    if (!any(p)) Inf else min(X[i, p])
  }, numeric(1))
}

#' @rdname batch-array-ops
#' @noRd
.rowMaxB <- function(X) {
  if (is_torch_tensor(X)) {
    return(as.numeric(toRMatrix(X$amax(dim = 2))))
  }
  apply(X, 1L, max)
}

#' The boundary back to the host
#'
#' \code{.polishBatch()} may take its counts and design as tensors, but it
#' always RETURNS host matrices and vectors: \code{.polishFit()} reads
#' \code{r$alpha[j, ]} per gene and \code{@polish} is a data frame. Keeping
#' the returned coefficients on the host also costs nothing -- they are
#' genes x columns, the one small array in the loop.
#'
#' @rdname batch-array-ops
#' @noRd
.asHostMat <- function(X) {
  if (is_torch_tensor(X)) return(as.matrix(toRMatrix(X)))
  X
}

#' @rdname batch-array-ops
#' @noRd
.asLike <- function(M, ref) {
  if (!is_torch_tensor(ref)) return(M)
  if (is_torch_tensor(M)) return(M)
  torch::torch_tensor(as.matrix(M), dtype = ref$dtype, device = ref$device)
}

#' Restrict a batched factorisation state to a subset of its genes
#'
#' The active set only ever shrinks inside \code{newton()}, so a gene leaving
#' it is a subset of the stack rather than a reason to refactorise. \code{xi}
#' and \code{zi} are column indices, shared by every slice, and must not be
#' touched.
#'
#' @param st a state from \code{.newtonSolverBatch()$factor()}.
#' @param ii positions to keep, within the state's current order.
#' @return the state restricted to \code{ii}.
#' @noRd
.subsetState <- function(st, ii) {
  sub1 <- function(x) {
    if (is.null(x)) return(NULL)
    if (is_torch_tensor(x)) {
      idx <- torch::torch_tensor(as.integer(ii), dtype = torch::torch_long(),
                                 device = x$device)
      return(torch::torch_index_select(x, 1, idx))
    }
    if (is.array(x) && length(dim(x)) == 3L) return(x[ii, , , drop = FALSE])
    if (is.matrix(x)) return(x[ii, , drop = FALSE])
    x[ii]
  }
  st$S <- sub1(st$S)
  st$L <- sub1(st$L)
  st$B <- sub1(st$B)
  st$cvec <- sub1(st$cvec)
  st$ok <- st$ok[ii]
  st
}

#' A batched .newtonSolver(): one factorisation object for a block of genes
#'
#' \code{.newtonSolver()} returns closures over a single gene's weights and
#' \code{.polishBatch()}'s \code{newton()} keeps a LIST of them, one per gene,
#' refreshed under a per-gene staleness counter. That list is what cannot go to
#' a device, and a per-gene \code{solve()} is a kernel launch per gene per
#' iteration. This is the same mathematics with the state batched: one
#' \code{(batch, px, px)} Cholesky for the whole block.
#'
#' Measured before it was written (FINDINGS, 2026-09-16): refreshing the whole
#' active stack whenever any gene is stale costs 11% more factorisations at a
#' 128-gene batch and 6% at 64, so a shared factorisation keeps essentially all
#' of Phase 0b's memoisation. That is why this returns one object rather than
#' trying to keep per-gene states on device.
#'
#' @param W the design (cells x p), a matrix or a torch tensor.
#' @param pen the length-\code{p} ridge penalty.
#' @param nested logical, length \code{p}: the nested indicator columns.
#' @return \code{list(factor, solve, xcov)}. \code{factor(wt_block)} takes the
#'   block's \code{batch x ncells} weights and returns a state carrying
#'   \code{ok}; \code{solve(state, Sc)} takes a \code{batch x p} score and
#'   returns the \code{batch x p} steps; \code{xcov(state)} returns the
#'   \code{(batch, px, px)} covariance on the dense columns.
#' @noRd
.newtonSolverBatch <- function(W, pen, nested = NULL) {
  # This path implements the DIAGONAL absorption only: C = Z' diag(w) Z with one
  # column per group, so C^-1 is a reciprocal and the whole thing batches as an
  # elementwise divide. .newtonSolver() also accepts a block grouping, for the
  # random-slope case where a sample's columns are not mutually orthogonal and
  # C^-1 is a per-block Cholesky; that has no batched equivalent written yet.
  # Refuse it here rather than let `which()` fail on a non-logical, because the
  # cost of guessing would be a wrong Newton step rather than an error.
  blk <- .absorbBlocks(nested, ncol(W))
  if (anyDuplicated(blk[!is.na(blk)])) {
    stop("the batched Newton solver absorbs 1x1 blocks only (one column per ",
         "group, C diagonal); this grouping has a block with more than one ",
         "column, which needs a per-block Cholesky. Use .newtonSolver() on the ",
         "CPU, or pass a logical `nested`.", call. = FALSE)
  }
  nested <- !is.na(blk)
  has_nested <- any(nested)
  xi <- which(!nested)
  zi <- which(nested)

  colsOf <- function(M, ii) {
    if (is_torch_tensor(M)) {
      idx <- torch::torch_tensor(as.integer(ii), dtype = torch::torch_long(),
                                 device = M$device)
      return(torch::torch_index_select(M, 2, idx))
    }
    M[, ii, drop = FALSE]
  }

  list(
    factor = function(wt_block, cell.tile = NULL) {
      if (!has_nested) {
        S <- .gramBatch(W, wt_block, penalty_diag = pen, cell.tile = cell.tile)
        ch <- .cholBatch(S)
        return(list(S = S, L = ch$L, ok = ch$ok, xi = xi, zi = zi,
                    B = NULL, cvec = NULL))
      }
      pr <- .absorbBatch(W, pen, nested, wt_block, cell.tile = cell.tile,
                         parts = TRUE)
      ch <- .cholBatch(pr$S)
      c(pr, list(L = ch$L, ok = ch$ok))
    },
    solve = function(state, Sc) {
      if (!length(state$zi)) return(.cholSolveBatch(state, Sc))
      s_x <- colsOf(Sc, state$xi)
      s_z <- colsOf(Sc, state$zi)
      if (is_torch_tensor(Sc)) {
        sc <- s_z / state$cvec
        rhs <- s_x - torch::torch_matmul(state$B$transpose(2, 3),
                                         sc$unsqueeze(3))$squeeze(3)
        dx <- .cholSolveBatch(state, rhs)
        dz <- (s_z - torch::torch_matmul(state$B, dx$unsqueeze(3))$squeeze(3)) /
          state$cvec
        out <- torch::torch_zeros(c(Sc$size(1), Sc$size(2)), dtype = Sc$dtype,
                                  device = Sc$device)
        ix <- torch::torch_tensor(as.integer(state$xi),
                                  dtype = torch::torch_long(), device = out$device)
        iz <- torch::torch_tensor(as.integer(state$zi),
                                  dtype = torch::torch_long(), device = out$device)
        out <- out$index_copy(2, ix, dx)$index_copy(2, iz, dz)
        return(out)
      }
      b <- nrow(Sc)
      rhs <- matrix(0, b, length(state$xi))
      for (g in seq_len(b)) {
        Bg <- matrix(state$B[g, , ], nrow = length(state$zi),
                     ncol = length(state$xi))          # G x px
        rhs[g, ] <- s_x[g, ] - crossprod(Bg, s_z[g, ] / state$cvec[g, ])
      }
      dx <- .cholSolveBatch(state, rhs)
      out <- matrix(NA_real_, b, ncol(Sc))
      for (g in seq_len(b)) {
        if (!state$ok[g]) next
        Bg <- matrix(state$B[g, , ], nrow = length(state$zi),
                     ncol = length(state$xi))
        out[g, state$xi] <- dx[g, ]
        out[g, state$zi] <- (s_z[g, ] - as.numeric(Bg %*% dx[g, ])) /
          state$cvec[g, ]
      }
      out
    },
    xcov = function(state) {
      if (is_torch_tensor(state$L)) {
        return(torch::torch_cholesky_inverse(state$L))
      }
      b <- dim(state$L)[1]
      p <- dim(state$L)[2]
      out <- array(NA_real_, c(b, p, p))
      for (g in seq_len(b)) {
        if (!state$ok[g]) next
        out[g, , ] <- chol2inv(matrix(state$L[g, , ], p, p))
      }
      out
    }
  )
}
