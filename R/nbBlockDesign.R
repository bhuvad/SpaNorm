# A compact representation of a design whose absorbed part is block-diagonal
# by cell group, and the Newton solver that uses it.
#
# The grouped absorption in .newtonSolver() (R/nbSolver.R) already handles a
# random block that is block-diagonal by sample -- columns sharing a group id
# form one dense block of C = Z' diag(w) Z, eliminated by a Schur complement.
# But it takes that block DENSE: an n x (G * q) matrix of which only n * q
# entries are non-zero, so the linear predictor, the score and the cross block
# B = X' diag(w) Z all run over G times more columns than carry information. At
# 300,000 cells, 45 groups and a 10-column block per group that is a 450-column
# dense block to carry 10 columns' worth of numbers. The compact form stores
# the block as it is: one n x q matrix Z (each cell's own block values) and a
# block id per cell. Every product is then O(n (p_x + q)) and the per-gene gram
# is O(n (p_x + q)^2), independent of the number of groups.

#' A design with a per-group block, in compact form
#'
#' Describes the design matrix
#' \deqn{W = [\, X \mid Z_1 \mid Z_2 \mid \dots \mid Z_G \,],}{W = [X | Z_1 | ... | Z_G],}
#' where \code{X} holds the ordinary (dense) columns and the per-group blocks
#' \eqn{Z_g} are the rows of \code{Z} belonging to group \eqn{g} and zero
#' elsewhere: every group has its own copy of the \code{q} columns of \code{Z}.
#' This is, for example, a per-patient intercept and per-patient spatial
#' library-size spline, \code{Z = cbind(1, l * tpsBasis(x, y, df = c(3, 3)))},
#' with \code{block} the patient. The linear predictor of a gene with
#' coefficients \code{alpha} is
#' \code{X \%*\% alpha_x + rowSums(Z * A[block, ])}, with \code{A} the
#' \code{G x q} matrix of that gene's block coefficients.
#'
#' The object is accepted as the design \code{W} by \code{\link{polishNB}()},
#' \code{\link{nbProfilePsi}()} and \code{\link{nbNewtonSolver}()}, which then
#' never build the dense \code{n x (p_x + G q)} matrix. The \code{G} blocks are
#' absorbed exactly by a Schur complement, as \code{polishNB(absorb = )} absorbs
#' a per-column grouping of a dense design, so a Newton step costs one
#' \code{n x (p_x + q)} gram whatever \code{G} is. The results are those of the
#' dense design \code{as.matrix(W)} with the block columns absorbed by group
#' (\code{absorb = c(rep(NA, p_x), rep(seq_len(G), each = q))}).
#'
#' \strong{Coefficient layout.} The implied design has
#' \code{p = p_x + G * q} columns, in the order of \code{as.matrix()}: the
#' columns of \code{X}, then group 1's \code{q} columns, then group 2's, and so
#' on (groups in the order of \code{levels(block)}). \code{alpha},
#' \code{lambda.a} and \code{start.cols} follow that layout; the per-group
#' coefficients of gene \code{g} are
#' \code{matrix(alpha[g, -(1:p_x)], G, q, byrow = TRUE)}. \code{lambda.a} and
#' \code{start.cols} may also be given for \code{[X | Z]} only (length
#' \code{p_x + q}), in which case the \code{Z} part is used for every group.
#'
#' \strong{Degenerate blocks.} Each group's block of the information matrix,
#' \eqn{C_g = Z_g' \mathrm{diag}(w) Z_g + \mathrm{diag}(\lambda_g)}{C_g =
#' Z_g' diag(w) Z_g + diag(lambda_g)}, is scaled to unit diagonal and
#' Cholesky-factorised. A group whose block is rank-deficient or nearly so --
#' fewer cells than columns, a column constant or zero over the group's cells,
#' a level of \code{block} with no cells at all, any of these with a zero
#' penalty -- is recognised by its smallest squared Cholesky pivot falling
#' below \code{rank.tol} (on the unit-diagonal scale that pivot is
#' \eqn{1 - R^2}{1 - R^2} of a column regressed on the columns before it, so
#' the default \code{1e-10} flags a variance inflation above \eqn{10^{10}}).
#' Such a block is inverted by a truncated eigendecomposition instead (a
#' generalised inverse; eigenvalues below \code{rank.tol} times the largest are
#' dropped), and a column with no weight and no penalty in the group is
#' dropped outright. The directions this drops are combinations of the
#' group's coefficients that do not change the linear predictor on its cells
#' and carry no penalty, so the objective is flat along them: the Newton step
#' is the minimum-norm step in the unit-diagonal scaling, which leaves the
#' coefficients along those directions at their starting values (an
#' unpenalised column of a group with no cells keeps its start exactly). The
#' likelihood, the fitted means, the dispersion, the dense-column
#' coefficients and their covariance are those of the same design with the
#' aliased columns removed; the individual block
#' coefficients of an aliased group are not identified and are one solution
#' among many. The dense grouped path has no such rule and reports a gene with
#' an unpenalised aliased block as singular.
#'
#' @param X a numeric cells x \code{p_x} matrix of the non-absorbed columns
#'   (\code{p_x >= 1}).
#' @param Z a numeric cells x \code{q} matrix: each cell's values of its own
#'   group's block columns.
#' @param block the group of each cell, a factor (or a vector coerced with
#'   \code{factor()}) with no missing values. Every level is a group, including
#'   a level with no cells, whose coefficients are then identified by the
#'   penalty alone (or held at their start, see Degenerate blocks).
#' @param rank.tol the tolerance below which a group's block is treated as
#'   rank-deficient (see Degenerate blocks).
#' @param x an \code{nbBlockDesign}.
#' @param ... unused.
#' @return an object of class \code{nbBlockDesign}: a list holding \code{X},
#'   \code{Z}, \code{block} and precomputed indexing, with \code{dim()} and
#'   \code{dimnames()} methods that describe the implied design (so
#'   \code{nrow()}, \code{ncol()} and \code{colnames()} work) and an
#'   \code{as.matrix()} method that builds it densely (for checking, or for a
#'   small problem).
#' @seealso \code{\link{polishNB}()}, \code{\link{nbNewtonSolver}()},
#'   \code{\link{tpsBasis}()}.
#' @examples
#' set.seed(1)
#' n <- 300
#' patient <- factor(sample(c("P1", "P2", "P3"), n, replace = TRUE))
#' X <- cbind(niche = rnorm(n))
#' l <- rnorm(n)                                  # centred log library size
#' Z <- cbind(1, l)                               # per-patient intercept and slope
#' D <- nbBlockDesign(X, Z, patient)
#' dim(D)                                         # 300 x (1 + 3 * 2)
#' mu <- exp(1 + 0.3 * X[, 1] + 0.5 * l)
#' Y <- t(replicate(4, rnbinom(n, mu = mu, size = 5)))
#' a0 <- matrix(0, 4, ncol(D))
#' fit <- polishNB(Y, D, a0, rep(0.2, 4), lambda.a = c(0, 1e-3, 1e-3),
#'                 start.cols = c(FALSE, TRUE, FALSE))
#' fit$alpha[, 1]                                 # the niche coefficient
#' @importFrom Matrix sparseMatrix
#' @export
nbBlockDesign <- function(X, Z, block, rank.tol = 1e-10) {
  X <- .asNumericMatrix(X, "X")
  Z <- .asNumericMatrix(Z, "Z")
  n <- nrow(X)
  if (ncol(X) < 1L) stop("'X' must have at least one column", call. = FALSE)
  if (ncol(Z) < 1L) stop("'Z' must have at least one column", call. = FALSE)
  if (nrow(Z) != n) {
    stop(sprintf("'X' and 'Z' must have one row per cell: nrow(X) = %d, nrow(Z) = %d",
                 n, nrow(Z)), call. = FALSE)
  }
  if (length(block) != n) {
    stop(sprintf("'block' must have one value per cell: length(block) = %d, nrow(X) = %d",
                 length(block), n), call. = FALSE)
  }
  if (anyNA(block)) stop("'block' must not contain missing values", call. = FALSE)
  if (!is.factor(block)) block <- factor(block)
  if (!is.numeric(rank.tol) || length(rank.tol) != 1L || !is.finite(rank.tol) ||
      rank.tol <= 0 || rank.tol >= 1) {
    stop("'rank.tol' must be a single number in (0, 1)", call. = FALSE)
  }
  px <- ncol(X)
  q <- ncol(Z)
  G <- nlevels(block)
  bid <- as.integer(block)
  # each group's cells (ascending), and its rows of X and Z held pre-split: the
  # per-gene factorisation scales each group's rows once and takes that
  # group's pieces of the information from them, with no row subsetting (which
  # measured 1.4-2.5x slower at 60,000-300,000 cells, 45 groups)
  cells <- split(seq_len(n), block)
  names(cells) <- NULL
  Xl <- lapply(cells, function(i) X[i, , drop = FALSE])
  Zl <- lapply(cells, function(i) Z[i, , drop = FALSE])
  # the block part of the implied design as a sparse n x (G q) matrix with
  # exactly n q stored entries: what the linear predictor and the score use
  Zsp <- Matrix::sparseMatrix(
    i = rep(seq_len(n), q),
    j = (rep(bid, q) - 1L) * q + rep(seq_len(q), each = n),
    x = as.numeric(Z), dims = c(n, G * q))
  xn <- colnames(X)
  if (is.null(xn)) xn <- paste0("X", seq_len(px))
  zn <- colnames(Z)
  if (is.null(zn)) zn <- paste0("Z", seq_len(q))
  cn <- c(xn, paste(rep(levels(block), each = q), rep(zn, G), sep = ":"))
  structure(list(
    X = X, Z = Z, block = block, bid = bid,
    cells = cells, Xl = Xl, Zl = Zl, Zsp = Zsp,
    n = n, px = px, q = q, G = G,
    xi = seq_len(px), zi = px + seq_len(G * q),
    colnames = cn, rank.tol = rank.tol
  ), class = "nbBlockDesign")
}

#' A numeric base matrix, or an error naming the argument
#' @noRd
.asNumericMatrix <- function(M, nm) {
  if (methods::is(M, "Matrix")) M <- as.matrix(M)
  if (is.null(dim(M))) M <- matrix(M, ncol = 1L)
  if (!is.matrix(M) || !is.numeric(M)) {
    stop(sprintf("'%s' must be a numeric matrix", nm), call. = FALSE)
  }
  if (!all(is.finite(M))) {
    stop(sprintf("'%s' must be finite (no NA, NaN or Inf)", nm), call. = FALSE)
  }
  storage.mode(M) <- "double"
  M
}

#' @rdname nbBlockDesign
#' @export
dim.nbBlockDesign <- function(x) c(x$n, x$px + x$G * x$q)

#' @rdname nbBlockDesign
#' @export
dimnames.nbBlockDesign <- function(x) list(NULL, x$colnames)

#' @rdname nbBlockDesign
#' @export
as.matrix.nbBlockDesign <- function(x, ...) {
  # filled group by group rather than through as.matrix(<sparse>), which warns
  # on every large coercion
  out <- matrix(0, x$n, x$px + x$G * x$q)
  out[, x$xi] <- x$X
  for (g in seq_len(x$G)) {
    if (length(x$cells[[g]])) {
      out[x$cells[[g]], x$px + (g - 1L) * x$q + seq_len(x$q)] <- x$Zl[[g]]
    }
  }
  colnames(out) <- x$colnames
  out
}

#' @rdname nbBlockDesign
#' @export
print.nbBlockDesign <- function(x, ...) {
  sz <- lengths(x$cells)
  cat(sprintf("nbBlockDesign: %d cells x %d columns (%d dense + %d groups x %d block columns)\n",
              x$n, x$px + x$G * x$q, x$px, x$G, x$q))
  cat(sprintf("  cells per group: min %d, median %g, max %d%s\n",
              min(sz), stats::median(sz), max(sz),
              if (any(sz == 0L)) sprintf(" (%d empty)", sum(sz == 0L)) else ""))
  invisible(x)
}

#' Is this design the compact block form?
#' @noRd
.isBlockDesign <- function(W) inherits(W, "nbBlockDesign")

#' The absorbed grouping of the dense design a compact one stands for
#'
#' \code{as.matrix(D)} absorbed with this grouping is the problem the compact
#' path solves; the tests use it as the oracle.
#' @noRd
.blockDenseAbsorb <- function(D) {
  c(rep(NA_integer_, D$px), rep(seq_len(D$G), each = D$q))
}

#' Expand a per-column vector given for X and Z to the implied design
#'
#' A vector of length \code{p_x + G q} is returned as is; one of length
#' \code{p_x + q} has its \code{Z} part repeated for every group; a single
#' value is recycled. Anything else is refused with \code{what} named.
#' @noRd
.blockExpand <- function(D, v, what) {
  p <- D$px + D$G * D$q
  if (length(v) == 1L) return(rep(v, p))
  if (length(v) == p) return(v)
  if (length(v) == D$px + D$q) {
    return(c(v[seq_len(D$px)], rep(v[D$px + seq_len(D$q)], D$G)))
  }
  stop(sprintf(paste0("'%s' must be a single value, one per column of the ",
                      "block design (%d), or one per column of [X | Z] (%d); ",
                      "%d supplied"), what, p, D$px + D$q, length(v)),
       call. = FALSE)
}

#' The cells on which column j of the implied design is non-zero, ascending
#' @noRd
.blockColCells <- function(D, j) {
  if (j <= D$px) return(which(D$X[, j] != 0))
  jj <- j - D$px - 1L
  g <- jj %/% D$q + 1L
  k <- jj %% D$q + 1L
  cells <- D$cells[[g]]
  cells[D$Zl[[g]][, k] != 0]
}

#' Linear predictor of one gene over the implied design (no offset)
#' @noRd
.blockEtaVec <- function(D, a) {
  as.numeric(D$X %*% a[D$xi]) + as.numeric(D$Zsp %*% a[D$zi])
}

#' t(W) r over the implied design, for one gene
#' @noRd
.blockScoreVec <- function(D, r) {
  c(as.numeric(crossprod(D$X, r)), as.numeric(Matrix::crossprod(D$Zsp, r)))
}

# Above this many genes the batched products run group by group through BLAS
# (tcrossprod(A_g, Z_g) and R[, cells_g] %*% Z_g) rather than through the
# sparse block matrix, whose kernel is not BLAS: measured at 300,000 cells, 45
# groups, q = 10, the per-group form is 2-4x faster at 32-69 genes, level at 3
# and slower for one gene.
BLOCK_BATCH_DENSE_MIN <- 5L

#' Linear predictors of a block of genes, genes x cells (no offset)
#' @noRd
.blockEtaBatch <- function(D, A) {
  if (is_torch_tensor(A)) .blockNoDevice()
  E <- tcrossprod(A[, D$xi, drop = FALSE], D$X)
  Az <- A[, D$zi, drop = FALSE]
  if (nrow(A) < BLOCK_BATCH_DENSE_MIN) {
    return(E + as.matrix(Matrix::tcrossprod(Az, D$Zsp)))
  }
  q <- D$q
  for (g in seq_len(D$G)) {
    cg <- D$cells[[g]]
    if (length(cg)) {
      E[, cg] <- E[, cg, drop = FALSE] +
        tcrossprod(Az[, (g - 1L) * q + seq_len(q), drop = FALSE], D$Zl[[g]])
    }
  }
  E
}

#' R W over the implied design, for a block of genes (R is genes x cells)
#' @noRd
.blockScoreBatch <- function(R, D) {
  if (is_torch_tensor(R)) .blockNoDevice()
  Sx <- R %*% D$X
  if (nrow(R) < BLOCK_BATCH_DENSE_MIN) return(cbind(Sx, as.matrix(R %*% D$Zsp)))
  q <- D$q
  Sz <- matrix(0, nrow(R), D$G * q)
  for (g in seq_len(D$G)) {
    cg <- D$cells[[g]]
    if (length(cg)) {
      Sz[, (g - 1L) * q + seq_len(q)] <- R[, cg, drop = FALSE] %*% D$Zl[[g]]
    }
  }
  cbind(Sx, Sz)
}

#' The compact block design runs on the CPU only
#' @noRd
.blockNoDevice <- function() {
  stop("an nbBlockDesign runs on the CPU only (backend = \"cpu\"); the ",
       "device path takes a dense design matrix", call. = FALSE)
}

# Design-agnostic operations for the per-gene engine. On a base matrix each is
# the expression .polishGene() has always evaluated, so the matrix path is
# unchanged bit for bit.
.designEta <- function(W, a) {
  if (.isBlockDesign(W)) .blockEtaVec(W, a) else W %*% a
}
.designCross <- function(W, r) {
  if (.isBlockDesign(W)) .blockScoreVec(W, r) else crossprod(W, r)
}
.designColNonzero <- function(W, j) {
  if (!.isBlockDesign(W)) return(W[, j] != 0)
  out <- logical(W$n)
  out[.blockColCells(W, j)] <- TRUE
  out
}

#' A generalised inverse half of one group's information block
#'
#' Returns \code{H} (\code{q x r}) with \code{C^- = H H'}: the exact inverse
#' through a Cholesky factor when the unit-diagonal-scaled block is well
#' conditioned, otherwise a truncated eigendecomposition (see
#' \code{nbBlockDesign()}, Degenerate blocks). \code{NULL} when the block is not
#' finite.
#' @noRd
.blockHalfInverse <- function(Cb, rank.tol) {
  q <- nrow(Cb)
  if (!all(is.finite(Cb))) return(NULL)
  d2 <- diag(Cb)
  keep <- d2 > 0
  H <- matrix(0, q, 0L)
  if (!any(keep)) return(H)
  d <- sqrt(d2[keep])
  k <- sum(keep)
  Cn <- Cb[keep, keep, drop = FALSE] / tcrossprod(d)
  R <- tryCatch(chol(Cn), error = function(e) NULL)
  if (!is.null(R) && all(is.finite(R)) && min(diag(R))^2 > rank.tol) {
    # Cb = D Cn D and Cn = R'R, so Cb^-1 = (D^-1 R^-1)(D^-1 R^-1)'
    Hk <- backsolve(R, diag(k)) / d
  } else {
    e <- eigen(Cn, symmetric = TRUE)
    ok <- e$values > rank.tol * max(e$values)
    Hk <- sweep(e$vectors[, ok, drop = FALSE], 2L, sqrt(e$values[ok]), `/`) / d
  }
  H <- matrix(0, q, ncol(Hk))
  H[keep, ] <- Hk
  H
}

#' Newton solver for the penalised NB information of a compact block design
#'
#' The same contract as \code{.newtonSolver()}: \code{factor(w)} returns the
#' penalised information at weights \code{w} as a state, \code{solve(w, s)}
#' the Newton step over all \code{p_x + G q} columns and \code{xcov(w)} the
#' dense-column block of the penalised covariance, the last two taking a
#' weight vector or a state and returning \code{NULL} on a singular Schur
#' complement. With \code{A = X' W X + diag(pen_x)},
#' \code{B_g = X_g' W_g Z_g} and \code{C_g = Z_g' W_g Z_g + diag(pen_g)}
#' (\code{X_g}, \code{Z_g}, \code{W_g} the group's rows), the step is
#' \code{dx = S^-1 (s_x - sum_g B_g C_g^- s_g)} with
#' \code{S = A - sum_g B_g C_g^- B_g'}, and
#' \code{dz_g = C_g^- (s_g - B_g' dx)}.
#'
#' The state carries \code{S} (\code{p_x x p_x}), \code{B} (\code{p_x x G q},
#' the groups' \code{B_g} side by side in the coefficient order -- the same
#' matrix the dense grouped solver's state calls \code{B}), \code{H} (a list of
#' \code{q x r_g} matrices with \code{C_g^- = H_g H_g'}) and \code{rank}
#' (\code{r_g}, per group).
#' @noRd
.newtonSolverCompact <- function(D, pen) {
  px <- D$px
  q <- D$q
  G <- D$G
  p <- px + G * q
  if (length(pen) != p) stop("internal: `pen` must be expanded to ncol(W)", call. = FALSE)
  xi <- D$xi
  pen_x <- pen[xi]
  pen_z <- matrix(pen[D$zi], q, G)
  cells <- D$cells
  Xl <- D$Xl
  Zl <- D$Zl
  rank.tol <- D$rank.tol
  zcols <- function(g) (g - 1L) * q + seq_len(q)

  parts <- function(w) {
    if (!all(is.finite(w))) return(NULL)
    # One pass over each group's rows: scaled once by sqrt(w), they give that
    # group's share of the dense gram A, its cross block B_g and its own block
    # C_g -- X' W Z = (X sqrt(w))' (Z sqrt(w)) -- every gram formed
    # symmetrically, as .newtonSolver() forms its own. A and the Schur
    # correction are accumulated separately and differenced once.
    A <- matrix(0, px, px)
    diag(A) <- pen_x
    corr <- matrix(0, px, px)
    Bm <- matrix(0, px, G * q)
    H <- vector("list", G)
    rk <- integer(G)
    for (g in seq_len(G)) {
      if (length(cells[[g]])) {
        sw <- sqrt(w[cells[[g]]])
        Xg <- Xl[[g]] * sw
        Zg <- Zl[[g]] * sw
        A <- A + crossprod(Xg)
        Bg <- crossprod(Xg, Zg)
        Cb <- crossprod(Zg)
      } else {
        Bg <- matrix(0, px, q)
        Cb <- matrix(0, q, q)
      }
      diag(Cb) <- diag(Cb) + pen_z[, g]
      Hg <- .blockHalfInverse(Cb, rank.tol)
      if (is.null(Hg)) return(NULL)
      if (ncol(Hg)) corr <- corr + tcrossprod(Bg %*% Hg)
      Bm[, zcols(g)] <- Bg
      H[[g]] <- Hg
      rk[g] <- ncol(Hg)
    }
    S <- A - corr
    structure(list(S = S, B = Bm, H = H, rank = rk),
              class = c("nbBlockFactor", "spiDE_nfac"))
  }
  as_fac <- function(x) {
    if (is.null(x) || inherits(x, "spiDE_nfac")) x else parts(x)
  }
  # C_g^- v for every group at once; v is q x G
  cinv <- function(st, V) {
    out <- matrix(0, q, G)
    for (g in seq_len(G)) {
      Hg <- st$H[[g]]
      if (ncol(Hg)) out[, g] <- Hg %*% crossprod(Hg, V[, g])
    }
    out
  }

  list(
    factor = parts,
    solve = function(w, s) {
      st <- as_fac(w)
      if (is.null(st)) return(NULL)
      sz <- matrix(s[D$zi], q, G)
      rhs <- s[xi] - as.numeric(st$B %*% as.numeric(cinv(st, sz)))
      dx <- tryCatch(solve(st$S, rhs), error = function(e) NULL)
      if (is.null(dx)) return(NULL)
      dz <- cinv(st, sz - matrix(crossprod(st$B, dx), q, G))
      c(as.numeric(dx), as.numeric(dz))
    },
    xcov = function(w) {
      st <- as_fac(w)
      if (is.null(st)) return(NULL)
      tryCatch(solve(st$S), error = function(e) NULL)
    }
  )
}
