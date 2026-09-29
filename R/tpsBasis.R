#' Thin-plate-style tensor spline basis over spatial coordinates
#'
#' The spatial basis SpaNorm builds its biology and library-size functions
#' from, exported so that a downstream model can use the same smooth functions
#' of position. It is the tensor product of a natural cubic spline basis
#' (\code{splines::ns()}) along each axis, with its columns centred.
#'
#' The basis is defined by a set of \emph{reference} coordinates (by default
#' the evaluation coordinates themselves) and can be evaluated anywhere:
#' \itemize{
#'   \item the per-axis degrees of freedom come from the reference ranges (see
#'     \code{df});
#'   \item each axis's interior knots sit at the quantiles of the reference
#'     coordinates along that axis and its boundary knots at their range,
#'     exactly as \code{splines::ns(ref, df = )} places them; a point outside
#'     the reference range is extrapolated linearly, as natural splines do;
#'   \item \strong{centring}: when \code{center = TRUE}, every column has the
#'     mean of that column over the \emph{reference} cells subtracted (the
#'     \code{scale(center = TRUE, scale = FALSE)} of the reference basis). The
#'     columns therefore average to zero over the reference cells, not over the
#'     evaluation cells: evaluated at a subset of the reference, each column is
#'     the same function of position that the whole reference sees, so two
#'     subsets of one section share one basis. The subtracted means are
#'     returned in the \code{"scaled:center"} attribute.
#' }
#'
#' With the default reference (\code{ref.x = x}, \code{ref.y = y}) and a single
#' \code{df}, the result is identical to the internal basis SpaNorm's fit
#' uses (\code{bs.tps()}). Note that \code{SpaNorm()} rescales each axis to
#' unit range before building its basis, so there the aspect-ratio rule below
#' gives \code{df} knots-worth of degrees of freedom on both axes; to reproduce
#' that on raw coordinates pass \code{df = c(df, df)}.
#'
#' @param x,y numeric vectors of the coordinates to evaluate the basis at.
#' @param df the degrees of freedom. A single positive integer applies
#'   SpaNorm's aspect-ratio rule to the reference ranges: with
#'   \code{gap = max(range_x, range_y) / df}, the axes get
#'   \code{ceiling(range_x / gap)} and \code{ceiling(range_y / gap)} degrees of
#'   freedom, so the longer axis gets \code{df}. Two positive integers give the
#'   x and y degrees of freedom directly. The basis has
#'   \code{df_x * df_y} columns (for example 9 at \code{df = c(3, 3)}).
#' @param ref.x,ref.y numeric vectors of the reference coordinates that define
#'   the knots, the per-axis degrees of freedom and the centring (for example
#'   the whole tissue section). Default: \code{x} and \code{y}.
#' @param center logical; subtract the reference column means (default
#'   \code{TRUE}).
#'
#' @return a \code{length(x)} by \code{df_x * df_y} numeric matrix. Column
#'   \code{(i - 1) * df_y + j} is the product of the \code{i}-th x-axis and the
#'   \code{j}-th y-axis natural spline. Attributes: \code{"df.tps"}, the
#'   per-axis degrees of freedom \code{c(df_x, df_y)}, and, when centred,
#'   \code{"scaled:center"}, the subtracted reference means.
#'
#' @examples
#' set.seed(1)
#' x <- runif(200, 0, 2)
#' y <- runif(200, 0, 1)
#' B <- tpsBasis(x, y, df = 4)          # 4 x 2 degrees of freedom, 8 columns
#' dim(B)
#' # the same basis, built on the whole section and evaluated on one region
#' sub <- x < 0.5
#' Bs <- tpsBasis(x[sub], y[sub], df = 4, ref.x = x, ref.y = y)
#' all.equal(Bs, B[sub, ], check.attributes = FALSE)
#' @export
tpsBasis <- function(x, y, df = 6, ref.x = x, ref.y = y, center = TRUE) {
  chk <- function(v, nm) {
    if (!is.numeric(v) || !length(v) || anyNA(v) || !all(is.finite(v))) {
      stop(sprintf("'%s' must be a non-empty finite numeric vector", nm),
           call. = FALSE)
    }
  }
  chk(x, "x"); chk(y, "y"); chk(ref.x, "ref.x"); chk(ref.y, "ref.y")
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length", call. = FALSE)
  }
  if (length(ref.x) != length(ref.y)) {
    stop("'ref.x' and 'ref.y' must have the same length", call. = FALSE)
  }
  if (!is.numeric(df) || !length(df) %in% c(1L, 2L) || anyNA(df) ||
      any(df <= 0) || any(df != round(df))) {
    stop("'df' must be one or two positive integers", call. = FALSE)
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("'center' must be TRUE or FALSE", call. = FALSE)
  }

  # per-axis degrees of freedom: bs.tps()'s aspect-ratio rule on the
  # reference ranges, or given directly
  if (length(df) == 1L) {
    xrng <- diff(range(ref.x))
    yrng <- diff(range(ref.y))
    gap <- max(xrng, yrng) / df
    df.x <- ceiling(xrng / gap)
    df.y <- ceiling(yrng / gap)
  } else {
    df.x <- as.integer(df[1])
    df.y <- as.integer(df[2])
  }
  if (!is.finite(df.x) || !is.finite(df.y) || df.x < 1 || df.y < 1) {
    stop("the reference coordinates span no range along an axis, so no ",
         "basis can be built along it; pass df = c(df_x, df_y) or other ",
         "reference coordinates", call. = FALSE)
  }

  # the per-axis bases on the reference, built exactly as bs.tps() builds them
  ns.x <- splines::ns(ref.x, df = df.x)
  ns.y <- splines::ns(ref.y, df = df.y)
  tensor <- function(bx, by) {
    out <- matrix(0, nrow = nrow(bx), ncol = df.x * df.y)
    for (i in seq_len(df.x)) {
      for (j in seq_len(df.y)) {
        out[, (i - 1) * ncol(by) + j] <- bx[, i] * by[, j]
      }
    }
    out
  }
  ref_tensor <- tensor(ns.x, ns.y)

  # at the reference itself the basis is the reference tensor; elsewhere the
  # same splines (same knots, same boundary knots) evaluated at the new points
  at_ref <- identical(x, ref.x) && identical(y, ref.y)
  bs.xy <- if (at_ref) ref_tensor else {
    tensor(unclass(stats::predict(ns.x, x)), unclass(stats::predict(ns.y, y)))
  }

  if (center) {
    # the reference column means, as scale(center = TRUE, scale = FALSE)
    # computes them and subtracts them: at the reference this is bs.tps()'s
    # own scale() call, bit for bit
    cm <- colMeans(ref_tensor, na.rm = TRUE)
    bs.xy <- sweep(bs.xy, 2L, cm, check.margin = FALSE)
    attr(bs.xy, "scaled:center") <- cm
  }
  attr(bs.xy, "df.tps") <- c(df.x, df.y)
  bs.xy
}
