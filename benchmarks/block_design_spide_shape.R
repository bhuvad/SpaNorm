# Timing benchmark (not a test): spiDE's depth = "spatial_spline" design for
# one index cell type, fitted by polishNB() on its DENSE grouped absorption
# (what spiDE does: W = [patient intercepts | depth blocks | covariates |
# niche | condition x niche], absorb = the patient of each block column) and
# on the same model as a compact nbBlockDesign.
#
# spiDE's block for patient s is [1 | l_1 | l_1 B_1 | l_2 | l_2 B_2 | ...]:
# one intercept, then per section j of the patient the log library size l_j
# (centred within the section) and l_j times the section's centred 3 x 3
# (2 x 2 below 200 cells, none below 100) natural-spline tensor basis B_j
# (spiDE's .sectionBasis(), equal to tpsBasis(df = c(d, d))). Sections are
# disjoint in cells. The compact form takes group = patient and Z = each
# cell's values of its patient's block columns, every patient's block padded
# with zero columns to the widest one; spiDE's 1e-3 ridge on the block columns
# keeps a padded column identified at 0, so the fit is the dense fit (checked
# below: coefficients, log-likelihoods, dispersions).
#
# Usage (from the package root, one BLAS thread):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript benchmarks/block_design_spide_shape.R \
#     n=30000 S=50 genes=1000 out=spide_shape.csv
# key=value arguments: n index cells (30000), S patients (50), maxsec sections
# per patient, drawn uniformly from 1..maxsec (3), niche columns nL (5; the
# condition model has 2 * nL dense columns), genes (1000), paths
# (intercept,dense,compact; "intercept" is the same model without the depth
# blocks -- patient intercepts absorbed as 1x1 blocks, spiDE's design without
# depth = "spatial_spline" -- as the reference cost), seed (1), lib (a library to load SpaNorm from; default the
# source tree via pkgload), spide (a spiDE source tree: if given, the dense
# design is also built with spiDE's own .indexDesign()/.depthBlocks() and
# checked equal to the one built here), out (a CSV to append rows to).
args <- commandArgs(trailingOnly = TRUE)
opt <- list(n = 30000, S = 50, maxsec = 3, nL = 5, genes = 1000,
            paths = "intercept,dense,compact", seed = 1, lib = "", spide = "", out = "")
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  opt[[kv[1]]] <- kv[2]
}
n <- as.integer(opt$n); S <- as.integer(opt$S); maxsec <- as.integer(opt$maxsec)
nL <- as.integer(opt$nL); ng <- as.integer(opt$genes)
paths <- strsplit(opt$paths, ",", fixed = TRUE)[[1]]
if (nzchar(opt$lib)) {
  .libPaths(c(opt$lib, .libPaths()))
  suppressPackageStartupMessages(library(SpaNorm))
} else {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
}
ns <- asNamespace("SpaNorm")
cat(sprintf("SpaNorm %s, BLAS threads %s\n", utils::packageVersion("SpaNorm"),
            if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_get_num_procs() else NA))

# ---- simulate one index type ---------------------------------------------
set.seed(as.integer(opt$seed))
sizes <- as.vector(rmultinom(1, n, rlnorm(S, 0, 0.6)))
patient <- factor(rep(sprintf("P%02d", seq_len(S)), sizes), levels = sprintf("P%02d", seq_len(S)))
nsec <- sample.int(maxsec, S, replace = TRUE)
# each patient's cells split over its sections, unequally
sec <- character(n)
for (s in seq_len(S)) {
  i <- which(as.integer(patient) == s)
  pr <- rgamma(nsec[s], 2); pr <- pr / sum(pr)
  sec[i] <- sprintf("P%02d_s%d", s, sample.int(nsec[s], length(i), replace = TRUE, prob = pr))
}
perm <- sample.int(n)                                # cells not sorted by patient
patient <- patient[perm]; sec <- sec[perm]
xy <- cbind(runif(n), runif(n) * 1.3)
# the section's reference coordinates: its index cells plus twice as many
# cells of other types (spiDE builds the basis on the whole section)
cells_of <- split(seq_len(n), sec)
ext_of <- lapply(cells_of, function(i) cbind(runif(2 * length(i)), runif(2 * length(i)) * 1.3))
ref_of <- lapply(names(cells_of), function(sg) rbind(xy[cells_of[[sg]], , drop = FALSE], ext_of[[sg]]))
names(ref_of) <- names(cells_of)
llib <- rnorm(n, 7.5, 0.6)                          # log library size
L <- matrix(log1p(rexp(n * nL, 0.5)), n, nL, dimnames = list(NULL, paste0("N", seq_len(nL))))
trt_p <- rep(0:1, length.out = S)
trt <- trt_p[as.integer(patient)]

# ---- spiDE's dense design, as .indexDesign() + .depthBlocks() build it ----
pid <- as.integer(patient)
blk_cols <- list(); grp <- integer()
for (s in seq_len(S)) {
  Zs <- list()
  for (sg in unique(sec[pid == s])) {
    i <- which(pid == s & sec == sg)
    l <- llib[i] - mean(llib[i])
    df <- if (length(i) >= 200L) 3L else if (length(i) >= 100L) 2L else 0L
    z <- matrix(0, n, 1L + if (df) df^2 else 0L)
    z[i, 1] <- l
    if (df) {
      R <- ref_of[[sg]]
      z[i, -1] <- l * tpsBasis(xy[i, 1], xy[i, 2], df = c(df, df), ref.x = R[, 1], ref.y = R[, 2])[, ]
    }
    Zs[[length(Zs) + 1L]] <- z
  }
  blk_cols[[s]] <- do.call(cbind, Zs)
  grp <- c(grp, rep(s, ncol(blk_cols[[s]])))
}
P <- stats::model.matrix(~ 0 + patient)
Zdep <- do.call(cbind, blk_cols)
TN <- L * trt
W <- cbind(P, Zdep, L, TN)
nb <- S + ncol(Zdep)
absorb <- c(seq_len(S), grp, rep(NA_integer_, ncol(W) - nb))
pen <- c(rep(1e-3, nb), rep(0, ncol(W) - nb))
start <- c(rep(TRUE, S), rep(FALSE, ncol(W) - S))
block_cols <- split(seq_len(nb), factor(c(seq_len(S), grp), levels = seq_len(S)))
dense <- (nb + 1L):ncol(W)
cat(sprintf("design: %d cells, %d patients, %d sections (%s per patient), %d block columns (%d-%d per patient), %d dense columns\n",
            n, S, length(unique(sec)), paste(range(nsec), collapse = "-"), nb,
            min(lengths(block_cols)), max(lengths(block_cols)), length(dense)))

if (nzchar(opt$spide)) {
  # the same design from spiDE's own code (sourced, read only)
  env <- new.env()
  sys.source(file.path(opt$spide, "R", "fit-design.R"), envir = env)
  # all cells: the index cells first, then every section's other cells, so
  # the reference of section sg is its index cells plus its other cells
  xy_all <- rbind(xy, do.call(rbind, ext_of))
  sec_all <- c(sec, rep(names(ext_of), vapply(ext_of, nrow, 1L)))
  ik <- seq_len(n)
  blk <- env$.depthBlocks(llib, xy_all, sec_all, ik, patient, L)
  des <- env$.indexDesign(L, NULL, patient, trt = trt, tested = colnames(L)[-nL], blocks = blk)
  cat(sprintf("spiDE's own design: same W %s, same absorb %s, same pen %s\n",
              isTRUE(all.equal(unname(des$W), unname(W), tolerance = 1e-12)),
              identical(des$absorb, absorb), identical(des$pen, pen)))
}

# ---- the same model, compact ----------------------------------------------
q <- max(lengths(block_cols))
Zc <- matrix(0, n, q)
for (s in seq_len(S)) {
  i <- which(pid == s)
  Zc[i, seq_along(block_cols[[s]])] <- W[i, block_cols[[s]]]
}
t_design <- system.time(D <- nbBlockDesign(W[, dense, drop = FALSE], Zc, patient))[["elapsed"]]
# dense column j of spiDE's W -> compact column
map <- integer(ncol(W))
map[dense] <- seq_along(dense)
for (s in seq_len(S)) map[block_cols[[s]]] <- length(dense) + (s - 1L) * q + seq_along(block_cols[[s]])
pen_c <- c(rep(0, length(dense)), rep(1e-3, q))
start_c <- c(rep(FALSE, length(dense)), TRUE, rep(FALSE, q - 1L))

# ---- counts and the start spiDE uses (patient log means) -----------------
a_pat <- rnorm(S, 0, 0.5)
eta0 <- a_pat[pid] + as.numeric(L %*% rnorm(nL, 0, 0.1)) + 0.9 * (llib - ave(llib, sec))
gmean <- runif(ng, -3, 1.5)                          # sparse to moderately expressed
Y <- t(vapply(seq_len(ng), function(g) as.numeric(rnbinom(n, mu = exp(eta0 + gmean[g]), size = 2)),
              numeric(n)))
rownames(Y) <- paste0("g", seq_len(ng))
lm <- t(rowsum(t(Y), pid)) / rep(tabulate(pid, S), each = ng)
A0 <- matrix(0, ng, ncol(W), dimnames = list(rownames(Y), NULL))
A0[, seq_len(S)] <- log(lm + 1e-3)
A0c <- matrix(0, ng, ncol(D), dimnames = list(rownames(Y), NULL))
A0c[, map] <- A0

rows <- list()
record <- function(path, what, sec_) {
  rows[[length(rows) + 1L]] <<- data.frame(n = n, S = S, genes = ng, p_dense = ncol(W),
    p_x = length(dense), q = q, path = path, what = what, seconds = sec_)
  cat(sprintf("%-8s %-30s %10.3f\n", path, what, sec_))
}
tm <- function(expr) system.time(expr, gcFirst = TRUE)[["elapsed"]]
w <- runif(n, 0.2, 3)
res <- list()
for (path in paths) {
  if (path == "intercept") {
    Wi <- cbind(P, L, TN)
    ab_i <- c(rep(TRUE, S), rep(FALSE, 2L * nL))
    pen_i <- c(rep(1e-3, S), rep(0, 2L * nL))
    sol <- ns$.newtonSolver(Wi, pen_i, ab_i)
  } else if (path == "dense") {
    sol <- ns$.newtonSolver(W, pen, absorb)
    Wuse <- W
  } else {
    record(path, "build design", t_design)
    sol <- ns$.newtonSolverCompact(D, ns$.blockExpand(D, pen_c, "pen"))
    Wuse <- D
  }
  record(path, "factor(w) [gram + Schur]", median(replicate(3, tm(sol$factor(w)))))
  if (path == "intercept") {
    t_pol <- tm(fit <- polishNB(Y, Wi, A0[, c(seq_len(S), dense), drop = FALSE], rep(1, ng),
                                lambda.a = pen_i, absorb = ab_i, start.cols = ab_i,
                                psi.method = "profile"))
  } else if (path == "dense") {
    fit <- NULL
    t_pol <- tm(fit <- polishNB(Y, W, A0, rep(1, ng), lambda.a = pen, absorb = absorb,
                                start.cols = start, psi.method = "profile"))
  } else {
    t_pol <- tm(fit <- polishNB(Y, D, A0c, rep(1, ng), lambda.a = pen_c,
                                start.cols = start_c, psi.method = "profile"))
    fit$alpha <- fit$alpha[, map, drop = FALSE]     # back to spiDE's layout
  }
  record(path, sprintf("polishNB, %d genes", ng), t_pol)
  record(path, "polishNB per gene", t_pol / ng)
  record(path, "Newton iterations per gene", mean(fit$polish$iterations))
  record(path, "genes polished", sum(fit$polish$polished))
  res[[path]] <- fit
  rm(sol); invisible(gc())
}
if (all(c("dense", "compact") %in% names(res))) {
  d <- res$dense; cc <- res$compact
  ok <- d$polish$polished & cc$polish$polished
  cat(sprintf("agreement over %d genes polished on both paths: max|d alpha| %.2e (dense columns %.2e), max|d loglik|/|loglik| %.2e, max|d log psi| %.2e, same polished flags %s\n",
              sum(ok), max(abs(d$alpha[ok, ] - cc$alpha[ok, ])),
              max(abs(d$alpha[ok, dense] - cc$alpha[ok, dense])),
              max(abs(d$loglik[ok] - cc$loglik[ok]) / abs(d$loglik[ok])),
              max(abs(log(d$psi[ok]) - log(cc$psi[ok]))),
              identical(d$polish$polished, cc$polish$polished)))
}
out <- do.call(rbind, rows)
if (nzchar(opt$out)) {
  write.table(out, opt$out, sep = ",", row.names = FALSE,
              col.names = !file.exists(opt$out), append = file.exists(opt$out))
}
