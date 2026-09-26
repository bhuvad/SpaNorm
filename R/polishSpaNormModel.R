# Maps a SpaNorm-model fit onto the generic per-gene polish problem that
# polishNB() (R/polishNB.R) solves. Groundwork for polishSpaNorm().
#
# SpaNorm's mean is log mu_gi = gmean_g + a1 * W_i1 + sum_{j>=2} W_ij alpha_gj:
# W[, 1] is log library size, and its coefficient a1 is SHARED by every gene
# (all rows of alpha[, 1] are equal); gmean is a per-gene intercept held
# outside W; columns 2..p are ridge-penalised by wtype (biology ->
# lambda.a[1], ls -> lambda.a[2], batch -> 0), times ncells.

#' The penalty vector SpaNorm's fitters use
#'
#' Mirrors the penalty construction in \code{fitSpaNorm()}
#' (\code{R/mainSpaNorm.R}, and from 1.7.12 \code{fitSpaNormTechnical()}):
#' \code{biology -> lambda.a[1]}, \code{ls -> lambda.a[2]}, \code{batch ->
#' 0}, times \code{ncells}. Column 1 (the shared log library size) is not
#' penalised and is not in the per-gene problem at all, hence the leading 0.
#' @param fit a SpaNormFit.
#' @return a numeric vector of length \code{ncol(fit$W)}: \code{0} followed
#'   by the penalty for \code{fit$W[, -1]}'s columns.
#' @noRd
.spaNormPenalty <- function(fit) {
  lam <- rep_len(fit$lambda.a, 2L)
  v <- numeric(length(fit$wtype))
  v[fit$wtype == "biology"] <- lam[1]
  v[fit$wtype == "ls"] <- lam[2]
  c(0, v[-1] * fit$ncells)          # leading 0: the unpenalised gmean column
}

#' Map a SpaNorm-model fit onto the generic per-gene polish problem
#'
#' \code{polishNB()} solves a generic penalised-NB problem
#' (\code{X}/\code{pen}/\code{offset}/start); this maps a SpaNorm-model fit
#' (\code{fitSpaNorm()}'s \code{alpha}/\code{gmean}/\code{W}/\code{wtype}
#' parameterisation) onto it. The shared library-size coefficient
#' (\code{alpha[, 1]}, identical across genes) is pulled out of \code{W} into
#' a per-cell offset rather than treated as a per-gene column, and the
#' per-gene intercept \code{gmean} becomes an explicit, unpenalised
#' \code{"(gmean)"} column of \code{X}.
#' @param fit a SpaNormFit.
#' @param cells which cells/spots to include: \code{"all"} (default), or
#'   \code{"fit"} for the cells the model was fit on (\code{sampling != "all"},
#'   i.e. \code{"glm"} and \code{"dispersion"}).
#' @return a list with \code{X}, \code{pen}, \code{a1} (the shared
#'   library-size coefficient, a scalar), \code{offset}, \code{A0} (starting
#'   coefficients, genes x ncol(X)), \code{psi}, \code{cells_idx} (logical,
#'   length \code{fit$ncells}) and \code{w1} (the selected cells' log library
#'   size, i.e. \code{W[cells_idx, 1]}).
#' @noRd
.spaNormPolishProblem <- function(fit, cells = c("all", "fit")) {
  cells <- match.arg(cells)
  a1s <- fit$alpha[, 1]
  if (diff(range(a1s)) > 1e-8 * max(1, abs(a1s[1]))) {
    stop("column 1 of the fit's alpha (the library-size coefficient) is not ",
         "shared across genes, so this is not a SpaNorm-model fit; use ",
         "polishNB() for a generic fitNB() fit", call. = FALSE)
  }
  idx <- if (cells == "all") rep(TRUE, nrow(fit$W)) else fit$sampling != "all"
  W <- fit$W[idx, , drop = FALSE]
  a1 <- a1s[1]
  list(X = cbind(`(gmean)` = 1, W[, -1, drop = FALSE]),
       pen = .spaNormPenalty(fit), a1 = a1, offset = a1 * W[, 1], w1 = W[, 1],
       A0 = cbind(fit$gmean, fit$alpha[, -1, drop = FALSE]),
       psi = fit$psi, cells_idx = idx)
}
