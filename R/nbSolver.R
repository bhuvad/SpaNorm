# Per-gene and per-block penalised NB Newton solvers, moved from spiDE
# (R/polish.R): .absorbBlocks, .newtonSolverBlocked, .newtonSolver. These form
# the penalised information W' diag(w) W + diag(pen) and solve it, absorbing a
# nested indicator (or per-sample block-diagonal random) block by a Schur
# complement so the per-iteration cost is one dense-column gram regardless of
# how many groups exist.

#' Normalise an absorption specification into a per-column block id
#'
#' Returns an integer vector, \code{NA} for a column that stays in the dense
#' block and otherwise the id of the block of \code{C = Z' diag(w) Z} the
#' column belongs to.
#' @noRd
.absorbBlocks <- function(nested, p) {
  if (is.null(nested)) return(rep(NA_integer_, p))
  if (length(nested) != p) {
    stop("`nested` must have one entry per column of W", call. = FALSE)
  }
  if (is.logical(nested)) {
    out <- rep(NA_integer_, p)
    out[which(nested)] <- which(nested)   # each absorbed column its own block
    return(out)
  }
  out <- rep(NA_integer_, p)
  keep <- !is.na(nested)
  if (any(keep)) out[keep] <- as.integer(factor(as.character(nested[keep])))
  out
}

#' Schur absorption of a random block that is block-diagonal, not diagonal
#'
#' Same identity as the scalar path, with \code{C^-1} a per-block solve. The
#' blocks must be mutually orthogonal under any weight, which holds iff no cell
#' loads on two blocks -- checked here rather than assumed, because a violated
#' assumption would not error, it would return a wrong step.
#' @noRd
.newtonSolverBlocked <- function(W, pen, blk, xi, zi, X, pen_x, pen_z) {
  bid <- blk[zi]
  cols <- split(seq_along(zi), bid)              # positions WITHIN zi
  Z <- W[, zi, drop = FALSE]

  # No cell may load on two blocks. Counting hits per cell is O(n x p_z), the
  # same order as the gram that follows, and it is done once per solver.
  hits <- integer(nrow(W))
  for (cc in cols) {
    hits <- hits + (rowSums(abs(Z[, cc, drop = FALSE])) > 0)
  }
  if (max(hits) > 1) {
    stop("the absorbed random-effect columns are not block-orthogonal: ",
         sum(hits > 1), " cells load on more than one block, so ",
         "C = Z'WZ is not block-diagonal and .newtonSolver() cannot absorb it",
         call. = FALSE)
  }

  parts <- function(w) {
    sw <- sqrt(w)
    Xw <- X * sw
    A <- crossprod(Xw)
    diag(A) <- diag(A) + pen_x
    Zw <- Z * sw
    B <- crossprod(Xw, Zw)                       # px x p_z
    # one dense C block per group, Cholesky-factorised once per weight vector
    fac <- lapply(cols, function(cc) {
      Cb <- crossprod(Zw[, cc, drop = FALSE])
      diag(Cb) <- diag(Cb) + pen_z[cc]
      tryCatch(chol(Cb), error = function(e) NULL)
    })
    if (any(vapply(fac, is.null, logical(1)))) return(NULL)
    S <- A
    for (k in seq_along(cols)) {
      Bb <- B[, cols[[k]], drop = FALSE]
      S <- S - Bb %*% chol2inv(fac[[k]]) %*% t(Bb)
    }
    structure(list(S = S, B = B, fac = fac, cols = cols), class = "spiDE_nfac")
  }
  as_fac <- function(x) if (inherits(x, "spiDE_nfac")) x else parts(x)
  # C^-1 v, block by block
  cinv <- function(p, v) {
    out <- numeric(length(v))
    for (k in seq_along(p$cols)) {
      cc <- p$cols[[k]]
      out[cc] <- backsolve(p$fac[[k]],
                           backsolve(p$fac[[k]], v[cc], transpose = TRUE))
    }
    out
  }

  list(
    factor = parts,
    solve = function(w, s) {
      p <- as_fac(w)
      if (is.null(p)) return(NULL)
      rhs <- s[xi] - as.numeric(p$B %*% cinv(p, s[zi]))
      dx <- tryCatch(solve(p$S, rhs), error = function(e) NULL)
      if (is.null(dx)) return(NULL)
      dz <- cinv(p, s[zi] - as.numeric(crossprod(p$B, dx)))
      out <- numeric(ncol(W))
      out[xi] <- dx
      out[zi] <- dz
      out
    },
    xcov = function(w) {
      p <- as_fac(w)
      if (is.null(p)) return(NULL)
      tryCatch(solve(p$S), error = function(e) NULL)
    }
  )
}

#' Newton solver for the penalised NB information, with indicator columns absorbed
#'
#' The nested (sample x cell type) random intercepts are 0/1 indicators that
#' partition the cells, so their blocks of the information matrix are cheap:
#' with per-cell weights \code{w} and the design split into dense \code{X} and
#' indicators \code{Z}, \code{A = X' diag(w) X + diag(pen_x)},
#' \code{B = X' diag(w) Z} is one \code{rowsum} pass, and
#' \code{C = Z' diag(w) Z + diag(pen_z)} is DIAGONAL. The step then comes from
#' the Schur complement \code{S = A - B C^-1 B'}, whose cost is one
#' \code{ncol(X)}-column gram regardless of how many groups there are -- the
#' difference between a 345-column and a 1,000-column gram per iteration on a
#' real cohort. \code{S^-1} is also the X-block of the full penalised
#' covariance, which is what the tested columns need.
#'
#' @param W the full design (cells x columns).
#' @param pen the per-column ridge penalty.
#' @param nested which columns to absorb, and how they block: a logical over
#'   the columns of \code{W} marking the indicator block (each absorbed column
#'   its own 1x1 block, \code{C} diagonal), or a per-column grouping
#'   (integer, factor or character block ids, \code{NA} for a dense column)
#'   whose columns sharing an id form one dense block of \code{C} -- the
#'   per-sample grouping of a random-slope fit's whole random block, from
#'   \code{.absorbSpec()}. All-FALSE, all-NA or \code{NULL} selects the plain
#'   dense path.
#' @return a list with \code{factor(w)} (the penalised information at those
#'   weights, as an opaque state), \code{solve(w, s)} (the Newton step over all
#'   columns) and \code{xcov(w)} (the X-block of the penalised covariance).
#'   \code{solve()} and \code{xcov()} take either a weight vector or a state
#'   from \code{factor()}; the first two return \code{NULL} on a singular
#'   system.
#'
#' Why \code{factor()} is separate: the caller already knows when the weights
#' change. \code{.polishGene()} refreshes them only every third Newton step, but
#' until 2026-09-15 it still paid for a full rebuild on every step, because
#' \code{solve()} recomputed the information from the (unchanged) weights it was
#' handed -- the gram, the group sums and the Schur complement, which the
#' documentation below calls "essentially the whole cost of an iteration".
#' Splitting the two lets a caller reuse a factorisation it has already paid
#' for. It is exactly reproducible: the state is built from the same weights
#' \code{solve()} would have used, and the step is the same \code{solve(S, rhs)}.
#' @noRd
.newtonSolver <- function(W, pen, nested = NULL) {
  # `nested` says which columns are absorbed and, now, how they BLOCK.
  #   NULL / all-FALSE  -> nothing absorbed, dense solve
  #   logical           -> each TRUE column is its own 1x1 block. This is the
  #                        nested (sample x cell type) intercept: 0/1 indicators
  #                        partitioning the cells, so C is diagonal.
  #   integer/factor/character with NA for the dense columns
  #                     -> columns sharing a level form ONE dense block of C.
  #                        Random slopes need this: a sample's slope columns are
  #                        indicator x covariate, so they are not orthogonal to
  #                        each other nor to that sample's intercept, and C is
  #                        block-diagonal by sample rather than diagonal.
  # The absorption identity is the same either way; only C^-1 changes, from a
  # reciprocal to a per-block solve.
  blk <- .absorbBlocks(nested, ncol(W))
  nested <- !is.na(blk)
  if (!any(nested)) {
    factor_dense <- function(w) {
      info <- crossprod(W * sqrt(w))
      diag(info) <- diag(info) + pen
      structure(list(info = info), class = "spiDE_nfac")
    }
    as_fac <- function(x) if (inherits(x, "spiDE_nfac")) x else factor_dense(x)
    return(list(
      factor = factor_dense,
      solve = function(w, s) {
        tryCatch(solve(as_fac(w)$info, s), error = function(e) NULL)
      },
      xcov = function(w) {
        tryCatch(solve(as_fac(w)$info), error = function(e) NULL)
      }
    ))
  }

  xi <- which(!nested)
  zi <- which(nested)
  X <- W[, xi, drop = FALSE]
  pen_x <- pen[xi]
  pen_z <- pen[zi]
  # A block carrying more than one column cannot use the reciprocal path.
  if (anyDuplicated(blk[zi])) {
    return(.newtonSolverBlocked(W, pen, blk, xi, zi, X, pen_x, pen_z))
  }
  # The indicator each cell belongs to. The columns are 0/1 and partition the
  # cells, so this dot product recovers the group index -- via round(), not
  # as.integer(): a floating-point product of 7 can come back as 6.9999999,
  # which as.integer() would truncate to the WRONG group, silently.
  Zblk <- W[, zi, drop = FALSE]
  # the absorption is exact only if C = Z' diag(w) Z is diagonal, i.e. only if
  # every cell belongs to exactly one group. Check that before relying on it.
  rs <- rowSums(Zblk)
  if (anyNA(rs) || max(abs(rs - 1)) > 1e-8) {
    stop("the nested random-effect columns are not 0/1 indicators partitioning ",
         "the cells; .newtonSolver() cannot absorb them", call. = FALSE)
  }
  gidx <- round(as.numeric(Zblk %*% seq_along(zi)))
  gf <- factor(gidx, levels = seq_along(zi))

  parts <- function(w) {
    # The gram is formed SYMMETRICALLY, crossprod(X * sqrt(w)) rather than
    # crossprod(X, X * w). This reverses an earlier decision, and the
    # measurement is the reason (FINDINGS.md, 2026-09-15, at the cohort shape
    # n = 77,454, px = 398):
    #
    #                       1 BLAS thread   8 BLAS threads
    #   crossprod(X, X*w)       0.426 s         2.436 s
    #   crossprod(X*sqrt(w))    0.247 s         0.149 s
    #
    # The old comment was right that the memory traffic, not the flop count, is
    # what this costs -- and wrong about which form pays it. dsyrk halves the
    # flops and, more to the point, SCALES with threads where the dgemm form
    # degrades: 0.43 -> 2.44 s as threads are added, against 0.25 -> 0.15 s.
    # That is the mechanism behind the recorded "one worker gains only 1.5x
    # from four BLAS threads", and it is why the polish has been run 64 workers
    # x 1 thread. It also makes the gram exactly symmetric, which the Cholesky
    # downstream assumes.
    #
    # The sqrt-weighted copy is then re-weighted in place to give X * w for the
    # group sums, so the peak is the same two n x px temporaries the previous
    # form reached through `Xw` plus the rowsum's own copy.
    sw <- sqrt(w)
    Xw <- X * sw
    A <- crossprod(Xw)
    diag(A) <- diag(A) + pen_x
    Xw <- Xw * sw                                    # now X * w
    cvec <- as.numeric(rowsum(w, group = gf, reorder = TRUE)) + pen_z
    B <- t(rowsum(Xw, group = gf, reorder = TRUE))   # ncol(X) x G
    structure(list(S = A - B %*% (t(B) / cvec), B = B, cvec = cvec),
              class = "spiDE_nfac")
  }
  as_fac <- function(x) if (inherits(x, "spiDE_nfac")) x else parts(x)

  list(
    factor = parts,
    solve = function(w, s) {
      p <- as_fac(w)
      rhs <- s[xi] - as.numeric(p$B %*% (s[zi] / p$cvec))
      dx <- tryCatch(solve(p$S, rhs), error = function(e) NULL)
      if (is.null(dx)) return(NULL)
      dz <- (s[zi] - as.numeric(crossprod(p$B, dx))) / p$cvec
      out <- numeric(ncol(W))
      out[xi] <- dx
      out[zi] <- dz
      out
    },
    xcov = function(w) {
      tryCatch(solve(as_fac(w)$S), error = function(e) NULL)
    }
  )
}

#' Penalised NB Newton solvers and batched grams (for package developers)
#'
#' Low-level building blocks of \code{\link{polishNB}}, exported so that a
#' downstream package (e.g. spiDE) can form the penalised covariance of a
#' polished fit from the SAME factorisation the polish used. The information is
#' \code{crossprod(W, w * W) + diag(pen)}. \code{absorb} marks an indicator block whose
#' part of the information is block-diagonal and is eliminated by a Schur
#' complement: \code{NULL} (dense), a logical over \code{W}'s columns (each
#' marked column its own 1x1 block), or an integer block id per column
#' (\code{NA} = dense).
#'
#' @param W a cells x p design (base matrix, or torch tensor for the batch forms).
#' @param pen a length-p ridge penalty.
#' @param absorb see Description.
#' @param wt_block a genes x cells weight matrix (one row per gene).
#' @param penalty_diag \code{NULL}, or a length-p ridge penalty added to the
#'   diagonal of every gene's gram.
#' @param backend the resolved backend; unused on base R matrices, kept for
#'   symmetry with the other batched helpers.
#' @param cell.tile cells per accumulation tile, or \code{NULL} for all at
#'   once. A gram is a sum over cells, so every tiling gives the same result;
#'   a tile bounds the \code{batch x cells x p} weighted design on a device.
#' @param parts if \code{TRUE}, return the pieces a Newton step's
#'   back-substitution needs (\code{S}, \code{B}, \code{cvec}, \code{xi},
#'   \code{zi}) rather than the Schur complements \code{S} alone.
#' @return \code{nbNewtonSolver()}: a list of closures \code{factor(w)},
#'   \code{solve(state, s)} (the Newton step, \code{NULL} if singular) and
#'   \code{xcov(state)} (the dense-block covariance). The batch forms return
#'   the batched grams or Schur complements.
#' @examples
#' W <- cbind(1, seq(-1, 1, length.out = 20))
#' s <- nbNewtonSolver(W, pen = c(0, 1))
#' st <- s$factor(rep(1, 20))
#' s$solve(st, c(1, 0))
#' @keywords internal
#' @export
nbNewtonSolver <- function(W, pen, absorb = NULL) .newtonSolver(W, pen, absorb)

#' @rdname nbNewtonSolver
#' @keywords internal
#' @export
nbNewtonSolverBatch <- function(W, pen, absorb = NULL) .newtonSolverBatch(W, pen, absorb)

#' @rdname nbNewtonSolver
#' @keywords internal
#' @export
nbGramBatch <- function(W, wt_block, penalty_diag = NULL, backend = "cpu", cell.tile = NULL) {
  .gramBatch(W, wt_block, penalty_diag, backend, cell.tile)
}

#' @rdname nbNewtonSolver
#' @keywords internal
#' @export
nbAbsorbGramBatch <- function(W, pen, absorb, wt_block, cell.tile = NULL, parts = FALSE) {
  .absorbBatch(W, pen, absorb, wt_block, cell.tile, parts)
}
