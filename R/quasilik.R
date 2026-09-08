# Quasi-likelihood dispersion for per-gene NB GLMs, edgeR v4 style, on the CPU
# or on the torch backend.
#
# edgeR computes the adjusted deviance and effective df in C (src/ql_glm.c
# compute_adjust_vec, src/ql_weights.c compute_weight), which a Bioconductor
# package may not reach into. The definition is not the Chebyshev tables --
# those are a fast path -- it is the phi >= 4.001 branch: obtain the first two
# moments of the unit deviance under the fitted NB by direct summation over the
# pmf, then match them to a scaled chi-square. If d ~ c * chisq_nu then
# E[d] = c*nu and Var[d] = 2 c^2 nu, so
#
#     nu = 2 E[d]^2 / Var[d]        the observation's effective df   (w1)
#     1/c = 2 E[d] / Var[d]         rescales its unit deviance       (w0)
#
# and s2 = sum(d_i * w0_i) / sum((1 - h_i) * w1_i). The moments depend only on
# (mu, phi), so qlDispersion() evaluates them once on a shared (log mu, log phi)
# table and interpolates -- the per-gene work is then elementwise over
# genes x cells, which is what runs on the accelerator.

#' NB unit deviance
#'
#' The negative binomial unit deviance
#' \eqn{2\{y \log(y/\mu) - (y + 1/\phi) \log((y + 1/\phi)/(\mu + 1/\phi))\}},
#' elementwise over a genes x cells matrix, on the CPU or -- when \code{y}
#' or \code{mu} is a torch tensor -- on that tensor's device.
#'
#' @param y counts, genes x cells (matrix or torch tensor).
#' @param mu fitted means, same shape (matrix or torch tensor).
#' @param phi NB dispersion: a scalar, or one value per gene (row).
#' @return a genes x cells matrix, or a torch tensor when the input was one.
#' @examples
#' y <- matrix(rpois(20, 3), 4, 5)
#' nbUnitDeviance(y, mu = y + 0.5, phi = 0.1)
#' @export
nbUnitDeviance <- function(y, mu, phi) {
  if (is_torch_tensor(y) || is_torch_tensor(mu)) {
    ref <- if (is_torch_tensor(mu)) mu else y
    y.t <- .nb_to_tensor(y, ref)
    mu.t <- .nb_to_tensor(mu, ref)
    sz <- .nb_size_tensor(1 / as.numeric(phi), ref)
    pos <- y.t > 0
    ysafe <- torch::torch_where(pos, y.t, torch::torch_ones_like(y.t))
    a <- torch::torch_where(pos, y.t * (torch::torch_log(ysafe) - torch::torch_log(mu.t)),
                            torch::torch_zeros_like(y.t))
    return(2 * (a - (y.t + sz) * (torch::torch_log(y.t + sz) - torch::torch_log(mu.t + sz))))
  }
  y <- as.matrix(y); mu <- as.matrix(mu)
  size <- rep_len(1 / as.numeric(phi), nrow(y))     # recycles down the rows
  a <- ifelse(y > 0, y * log(pmax(y, .Machine$double.xmin) / mu), 0)
  2 * (a - (y + size) * log((y + size) / (mu + size)))
}

#' First two moments of the NB unit deviance, by direct summation
#'
#' Reproduces the definition edgeR's \code{compute_weight()} approximates with
#' Chebyshev tables. The summation window is centred on \code{mu} and sized
#' from the pmf's own quantiles rather than fixed at edgeR's 50 terms from
#' zero: edgeR can truncate there because that branch is only reached when
#' \code{phi >= 4}, where the pmf piles up near zero; this is reached with
#' any dispersion.
#'
#' Elements are grouped by the width they need so a few high-mu cells do not
#' impose their window on everything, and each group accumulates in two passes
#' so memory stays O(length(mu)) rather than O(width * length(mu)).
#'
#' @param mu vector of means.
#' @param phi NB dispersion: scalar, or the same length as \code{mu}.
#' @param eps tail probability left outside the summation window at each end.
#' @param maxterms hard cap on the number of pmf terms summed per element.
#' @return list with \code{w0} (deviance rescaling) and \code{w1} (effective
#'   df), each the length of \code{mu}.
#' @examples
#' nbDevianceMoments(mu = c(0.1, 1, 10), phi = 5)
#' @export
nbDevianceMoments <- function(mu, phi, eps = 1e-10, maxterms = 20000L) {
  mu <- as.numeric(mu)
  phi <- rep_len(as.numeric(phi), length(mu))
  w0 <- w1 <- numeric(length(mu))
  ok <- which(is.finite(mu) & mu > 1e-32 & is.finite(phi) & phi > 0)
  if (!length(ok)) return(list(w0 = w0, w1 = w1))

  # Window from the pmf's own quantiles, not from mu +/- k*sd: at a sparse
  # operating point (mu ~ 0.1, phi ~ 5, so size ~ 0.2) the NB is far too
  # heavy-tailed for an sd-based window -- at mu = 1, phi = 20 a 12-sd window
  # gets w1 wrong by 41%.
  sz <- 1 / phi[ok]
  lo <- stats::qnbinom(eps, size = sz, mu = mu[ok])
  hi <- stats::qnbinom(1 - eps, size = sz, mu = mu[ok])
  wid <- pmin(maxterms, pmax(10, hi - lo))
  grp <- split(seq_along(ok), pmin(30L, ceiling(log2(wid))))

  for (ix in grp) {
    m <- mu[ok][ix]; p <- phi[ok][ix]; size <- 1 / p
    l <- lo[ix]; K <- max(wid[ix])
    unit <- function(i) {
      a <- ifelse(i > 0, i * log(pmax(i, 1) / m), 0)
      2 * (a - (i + size) * log((i + size) / (m + size)))
    }
    ed <- vd <- numeric(length(ix))
    for (j in seq_len(K + 1L) - 1L) {           # pass 1: E[d]
      i <- l + j
      ed <- ed + stats::dnbinom(i, size = size, mu = m) * unit(i)
    }
    for (j in seq_len(K + 1L) - 1L) {           # pass 2: Var[d]
      i <- l + j
      vd <- vd + stats::dnbinom(i, size = size, mu = m) * (unit(i) - ed)^2
    }
    good <- is.finite(ed) & is.finite(vd) & vd > 0
    w0[ok[ix][good]] <- (2 * ed / vd)[good]
    w1[ok[ix][good]] <- (2 * ed * ed / vd)[good]
  }
  list(w0 = w0, w1 = w1)
}

#' A shared table of NB deviance moments
#'
#' The deviance moments \code{w0} and \code{w1} (see
#' \code{\link{nbDevianceMoments}}) evaluated once on a uniform grid over
#' log mean and log dispersion, for \code{\link{qlDispersion}} to interpolate
#' onto every gene and cell. A caller that scores genes in blocks builds one
#' table over the whole range and passes it to every block, so the result is
#' invariant to how the genes are split; lookups outside the range are
#' clamped to the edge, where the moments are asymptotically flat.
#'
#' @param lmu_range range of log means the table must cover, \code{c(lo, hi)}.
#' @param lphi_range range of log dispersions, \code{c(lo, hi)}; a single
#'   value gives a one-row table.
#' @param step_mu,step_phi grid spacing along log mean and log dispersion;
#'   bilinear interpolation error scales with the square of the spacing.
#' @param ngrid,nphi explicit grid sizes, overriding the steps.
#' @return a list holding the two moment matrices (log-phi rows, log-mu
#'   columns) and the grid geometry, for \code{qlDispersion(table = )}.
#' @examples
#' tab <- qlMomentTable(log(c(1e-6, 1e3)), log(c(0.1, 10)))
#' dim(tab$w1)
#' @export
qlMomentTable <- function(lmu_range, lphi_range, step_mu = 0.08, step_phi = 0.12,
                          ngrid = NULL, nphi = NULL) {
  lo <- min(lmu_range); hi <- max(max(lmu_range), lo + 1e-8)
  K <- if (is.null(ngrid)) max(2L, ceiling((hi - lo) / step_mu) + 1L) else as.integer(ngrid)
  kn <- seq(lo, hi, length.out = K)
  pr <- range(lphi_range)
  J <- if (diff(pr) < 1e-12) 1L
       else if (is.null(nphi)) max(2L, ceiling(diff(pr) / step_phi) + 1L) else as.integer(nphi)
  lp <- if (J == 1L) pr[1] else seq(pr[1], pr[2], length.out = J)
  mk <- nbDevianceMoments(rep(exp(kn), J), rep(exp(lp), each = K))
  list(w0 = matrix(mk$w0, J, K, byrow = TRUE), w1 = matrix(mk$w1, J, K, byrow = TRUE),
       lo = lo, h = (hi - lo) / (K - 1), K = K, plo = lp[1],
       ph = if (J > 1L) (lp[J] - lp[1]) / (J - 1) else 1, J = J)
}

# the per-call table qlDispersion() builds when the caller supplies none:
# 256 points over this block's own log-mu range and up to 64 over its phi
.qlMomentTable <- function(lmu_range, lphi, ngrid, nphi, step_phi = 0.1) {
  pr <- range(lphi)
  J <- if (diff(pr) < 1e-12) 1L else min(nphi, max(2L, ceiling(diff(pr) / step_phi) + 1L))
  qlMomentTable(lmu_range, pr, ngrid = ngrid, nphi = J)
}

# per-element (row-wise phi) bilinear lookup; `lmu` a matrix or tensor of log
# means, `lphi` one value per row. Rule-2 extrapolation (clamped) at the edges.
.qlInterp <- function(tab, lmu, lphi) {
  K <- tab$K; J <- tab$J
  xp <- if (J > 1L) pmin(pmax((lphi - tab$plo) / tab$ph, 0), J - 1) else rep(0, length(lphi))
  jp <- pmin(floor(xp), max(J - 2L, 0L)); fp <- if (J > 1L) xp - jp else rep(0, length(lphi))
  if (is_torch_tensor(lmu)) {
    ref <- lmu; nr <- dim(lmu)[[1]]
    x <- torch::torch_clamp((lmu - tab$lo) / tab$h, min = 0, max = K - 1)
    i <- torch::torch_clamp(torch::torch_floor(x), min = 0, max = K - 2)
    fx <- x - i
    jt <- torch::torch_tensor(jp, dtype = ref$dtype, device = ref$device)$view(c(nr, 1L))
    fpt <- torch::torch_tensor(fp, dtype = ref$dtype, device = ref$device)$view(c(nr, 1L))
    gather <- function(M, jj, ii) {              # M is J x K; linear index j + i*J (1-based)
      flat <- torch::torch_tensor(as.numeric(M), dtype = ref$dtype, device = ref$device)
      idx <- (jj + ii * J + 1)$to(dtype = torch::torch_long())
      torch::torch_index_select(flat, 1, idx$flatten())$view(dim(lmu))
    }
    lerp <- function(M) {
      (1 - fpt) * ((1 - fx) * gather(M, jt, i) + fx * gather(M, jt, i + 1)) +
        fpt * ((1 - fx) * gather(M, pmin_t(jt + 1, J - 1), i) + fx * gather(M, pmin_t(jt + 1, J - 1), i + 1))
    }
    pmin_t <- function(a, b) torch::torch_clamp(a, max = b)
    return(list(w0 = lerp(tab$w0), w1 = lerp(tab$w1)))
  }
  nr <- nrow(lmu)
  x <- pmin(pmax((lmu - tab$lo) / tab$h, 0), K - 1)
  i <- pmin(floor(x), K - 2); fx <- x - i
  jm <- matrix(jp, nr, ncol(lmu)); fpm <- matrix(fp, nr, ncol(lmu))
  gather <- function(M, jj, ii) matrix(M[jj + ii * J + 1], nr, ncol(lmu))
  lerp <- function(M) {
    (1 - fpm) * ((1 - fx) * gather(M, jm, i) + fx * gather(M, jm, i + 1)) +
      fpm * ((1 - fx) * gather(M, pmin(jm + 1, J - 1), i) + fx * gather(M, pmin(jm + 1, J - 1), i + 1))
  }
  list(w0 = lerp(tab$w0), w1 = lerp(tab$w1))
}

#' Per-gene quasi-likelihood dispersion
#'
#' The edgeR v4 quasi-likelihood dispersion of each gene's NB GLM: the sum of
#' its unit deviances, each rescaled by the observation's deviance moments,
#' over the effective residual degrees of freedom. With \code{moments =
#' "table"} (the default) the moments are evaluated once on a shared
#' (log mu, log phi) table and interpolated, so the per-gene work is
#' elementwise and runs on the accelerator when \code{y} and \code{mu} are
#' torch tensors.
#'
#' @param y counts, genes x cells (matrix or torch tensor).
#' @param mu fitted means, same shape (matrix or torch tensor).
#' @param phi NB dispersion: scalar or one value per gene.
#' @param design the design matrix (cells x p); needed for
#'   \code{leverage = "exact"}, otherwise only its column count is used.
#' @param p the number of design columns, in place of \code{design}.
#' @param prior the average quasi-dispersion edgeR divides through by; 1
#'   leaves the parameterisation alone.
#' @param leverage \code{"trace"} spreads the p degrees of freedom evenly,
#'   which is exact to O(p/n) and is what n >> p designs want;
#'   \code{"exact"} forms per-observation hat values (O(n p^2) per gene, CPU
#'   only) and is what reproduces edgeR on its own small-n designs.
#' @param moments \code{"table"} evaluates the moments on a shared
#'   (log mu, log phi) table and interpolates onto every gene and cell (both
#'   backends); \code{"grid"} uses a per-gene log-mu grid; \code{"cell"}
#'   evaluates at every cell. The last two are CPU only.
#' @param ngrid grid points along log mu.
#' @param nphi maximum grid points along log phi for \code{moments = "table"}.
#' @param table a moments table from \code{\link{qlMomentTable}}, built once
#'   by a caller that scores genes in blocks; when \code{NULL} the table is
#'   built from this call's own range of means and dispersions.
#' @return a list of per-gene \code{deviance} (adjusted), \code{df}
#'   (effective) and \code{s2 = deviance / df}.
#' @examples
#' set.seed(1)
#' mu <- matrix(exp(rnorm(400, 0, 1)), 8, 50)
#' y <- matrix(rnbinom(400, mu = mu, size = 2), 8, 50)
#' qlDispersion(y, mu, phi = 0.5, p = 3)
#' @export
qlDispersion <- function(y, mu, phi, design = NULL, p = NULL, prior = 1,
                         leverage = c("trace", "exact"),
                         moments = c("table", "grid", "cell"),
                         ngrid = 256L, nphi = 64L, table = NULL) {
  leverage <- match.arg(leverage); moments <- match.arg(moments)
  on_device <- is_torch_tensor(y) || is_torch_tensor(mu)
  if (is.null(p)) {
    if (is.null(design)) stop("supply `design` or its column count `p`", call. = FALSE)
    p <- ncol(design)
  }
  if (on_device && (leverage != "trace" || moments != "table"))
    stop("on torch tensors qlDispersion() supports leverage = \"trace\" and moments = \"table\" only",
         call. = FALSE)
  ng <- if (on_device) dim(if (is_torch_tensor(y)) y else mu)[[1]] else nrow(as.matrix(y))
  n <- if (on_device) dim(if (is_torch_tensor(y)) y else mu)[[2]] else ncol(as.matrix(y))
  phi <- rep_len(as.numeric(phi), ng)
  phi_eff <- phi / prior                          # the deviance's dispersion

  if (moments == "table") {
    hdp <- 1 - p / n
    if (!is.finite(hdp) || hdp < 1e-4) return(list(deviance = numeric(ng), df = numeric(ng), s2 = numeric(ng)))
    if (on_device) {
      ref <- if (is_torch_tensor(mu)) mu else y
      mu.t <- .nb_to_tensor(mu, ref)
      lmu <- torch::torch_log(torch::torch_clamp(mu.t / prior, min = 1e-32))
      tab <- if (is.null(table)) {
        rng <- c(as.numeric(lmu$min()$cpu()), as.numeric(lmu$max()$cpu()))
        .qlMomentTable(rng, log(phi), ngrid, nphi)
      } else table
      w <- .qlInterp(tab, lmu, log(phi))
      udp <- nbUnitDeviance(y, mu.t, phi_eff)
      udp <- torch::torch_where(torch::torch_isfinite(udp), udp, torch::torch_zeros_like(udp))
      dev <- as.numeric(torch::torch_sum(udp * w$w0, dim = 2)$cpu())
      df <- hdp * as.numeric(torch::torch_sum(w$w1, dim = 2)$cpu())
    } else {
      mu <- as.matrix(mu)
      lmu <- log(pmax(mu / prior, 1e-32))
      tab <- if (is.null(table)) .qlMomentTable(range(lmu), log(phi), ngrid, nphi) else table
      w <- .qlInterp(tab, lmu, log(phi))
      udp <- nbUnitDeviance(y, mu, phi_eff)
      udp[!is.finite(udp)] <- 0
      dev <- rowSums(udp * w$w0)
      df <- hdp * rowSums(w$w1)
    }
    return(list(deviance = dev, df = df, s2 = ifelse(df < 1e-4, 0, dev / df)))
  }

  # per-gene paths: a log-mu grid per gene, or the moments at every cell
  y <- as.matrix(y); mu <- as.matrix(mu)
  dev <- df <- numeric(ng)
  for (g in seq_len(ng)) {
    mg <- mu[g, ]; yg <- y[g, ]
    h <- if (leverage == "exact") {
      zw <- sqrt(mg / (1 + mg * phi_eff[g]))
      rowSums(qr.Q(qr(design * zw))^2)
    } else rep(p / n, n)
    # moments are taken at mu/prior with the UNSCALED phi, matching edgeR's
    # compute_weight(u, phi, prior) call convention
    ms <- mg / prior
    m <- if (moments == "cell" || length(unique(ms)) <= ngrid) {
      nbDevianceMoments(ms, phi[g])
    } else {
      lg <- log(pmax(ms, 1e-32))
      kn <- seq(min(lg), max(lg), length.out = ngrid)
      mk <- nbDevianceMoments(exp(kn), phi[g])
      list(w0 = stats::approx(kn, mk$w0, lg, rule = 2)$y, w1 = stats::approx(kn, mk$w1, lg, rule = 2)$y)
    }
    udp <- nbUnitDeviance(yg, mg, phi_eff[g])
    hdp <- 1 - h
    drop <- !is.finite(hdp) | hdp < 1e-4
    udp[drop] <- 0; hdp[drop] <- 0
    dev[g] <- sum(udp * m$w0, na.rm = TRUE)
    df[g] <- sum(hdp * m$w1, na.rm = TRUE)
  }
  list(deviance = dev, df = df, s2 = ifelse(df < 1e-4, 0, dev / df))
}
