# Timing benchmark (not a test): polishNB() on a per-patient block design,
# dense grouped (as.matrix(D) with absorb = <grouping>) against the compact
# nbBlockDesign, per gene, with the intercept-only design (the per-patient
# intercepts alone, absorbed as 1x1 blocks: what spiDE's engines fit today)
# as the reference cost (and that same intercept-only model passed compactly,
# nbBlockDesign(X, matrix(1, n, 1), patient), as a fourth arm), at the shape of spiDE's "spatial_spline" depth
# handling: per patient s a block [1 | l * B_s(x, y)], B_s a centred natural
# spline tensor of df 3 x 3 on the patient's section, so q = 10 block columns
# per patient, plus p_x dense columns.
#
# Usage (from the package root, one BLAS thread):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript benchmarks/block_design_timing.R \
#     n=60000 paths=intercept,dense,compact genes=3 out=timing.csv
#
# Arguments (key=value): n cells (default 60000), S patients (45), px dense
# columns (40), genes (3), paths (intercept,intercept-compact,dense,compact), reps of the solver-level
# timings (3), seed (1), lib (a library to load SpaNorm from; default: the
# source tree via pkgload), out (a CSV to append rows to).
args <- commandArgs(trailingOnly = TRUE)
opt <- list(n = 60000, S = 45, px = 40, genes = 3,
            paths = "intercept,intercept-compact,dense,compact",
            reps = 3, seed = 1, lib = "", out = "")
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  opt[[kv[1]]] <- kv[2]
}
n <- as.integer(opt$n); S <- as.integer(opt$S); px <- as.integer(opt$px)
ng <- as.integer(opt$genes); reps <- as.integer(opt$reps)
paths <- strsplit(opt$paths, ",", fixed = TRUE)[[1]]
if (nzchar(opt$lib)) {
  .libPaths(c(opt$lib, .libPaths()))
  suppressPackageStartupMessages(library(SpaNorm))
  ns <- asNamespace("SpaNorm")
} else {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  ns <- asNamespace("SpaNorm")
}
blas <- if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_get_num_procs() else NA

# ---- simulate the design --------------------------------------------------
set.seed(as.integer(opt$seed))
w_pat <- rlnorm(S, 0, 0.6)
sizes <- as.vector(rmultinom(1, n, w_pat / sum(w_pat)))
patient <- factor(rep(sprintf("P%02d", seq_len(S)), sizes), levels = sprintf("P%02d", seq_len(S)))
patient <- patient[sample.int(n)]                   # cells not sorted by patient
xy <- cbind(runif(n), runif(n) * 1.3)               # each patient its own section
llib <- rnorm(n, 0, 0.6)                            # log library size ...
l <- llib - ave(llib, patient)                      # ... centred within patient
Bsec <- matrix(0, n, 9)
for (s in levels(patient)) {
  i <- which(patient == s)
  Bsec[i, ] <- tpsBasis(xy[i, 1], xy[i, 2], df = c(3, 3))[, ]
}
Z <- cbind(int = 1, l * Bsec)                        # q = 10
colnames(Z) <- c("int", paste0("ls", 1:9))
X <- matrix(abs(rnorm(n * px)) * 0.5, n, px)         # niche-like, non-negative
colnames(X) <- paste0("x", seq_len(px))
t_design <- system.time(D <- nbBlockDesign(X, Z, patient))[["elapsed"]]
q <- ncol(Z)
a_pat <- rnorm(S, 0, 0.5)
beta <- rnorm(px, 0, 0.05)
eta0 <- as.numeric(X %*% beta) + a_pat[as.integer(patient)] + 0.9 * l
Y <- t(vapply(seq_len(ng), function(g) as.numeric(rnbinom(n, mu = exp(eta0 + log(0.5 * g)), size = 2)),
              numeric(n)))
rownames(Y) <- paste0("g", seq_len(ng))
# the per-patient sane start, as spiDE's engines build it
pid <- as.integer(patient)
lm <- t(rowsum(t(Y), pid)) / rep(tabulate(pid, S), each = ng)
A0 <- matrix(0, ng, ncol(D), dimnames = list(rownames(Y), NULL))
A0[, px + (seq_len(S) - 1L) * q + 1L] <- log(lm + 1e-3)
pen <- c(rep(0, px), 1e-3, rep(1e-3, q - 1L))        # [X | Z]
start <- c(rep(FALSE, px), TRUE, rep(FALSE, q - 1L))
pen_full <- ns$.blockExpand(D, pen, "pen")
start_full <- ns$.blockExpand(D, start, "start")
grp <- ns$.blockDenseAbsorb(D)
w <- runif(n, 0.2, 3)

rows <- list()
record <- function(path, what, sec, extra = list()) {
  rows[[length(rows) + 1L]] <<- data.frame(
    n = n, S = S, px = px, q = q, p = pcur, path = path, what = what,
    seconds = sec, blas = blas, stringsAsFactors = FALSE)
  cat(sprintf("%-8s %-28s %9.3f s\n", path, what, sec))
}
tm <- function(expr) system.time(expr, gcFirst = TRUE)[["elapsed"]]
res <- list()
pcur <- NA_integer_
for (path in paths) {
  pcur <- if (path %in% c("intercept", "intercept-compact")) S + px else ncol(D)
  if (path == "intercept") {
    # the reference: [patient indicators | X], the indicators absorbed as 1x1
    # blocks -- no library-size spline at all, so a different model
    t_build <- tm(Wi <- cbind(stats::model.matrix(~ 0 + patient), X))
    record(path, "build design", t_build)
    nest <- c(rep(TRUE, S), rep(FALSE, px))
    Wuse <- Wi; absorb <- nest; lam <- c(rep(1e-3, S), rep(0, px)); st <- nest
    A0use <- cbind(A0[, px + (seq_len(S) - 1L) * q + 1L, drop = FALSE],
                   matrix(0, ng, px))
    solver_mk <- function() ns$.newtonSolver(Wi, lam, nest)
  } else if (path == "intercept-compact") {
    t_build <- tm(Di <- nbBlockDesign(X, matrix(1, n, 1, dimnames = list(NULL, "int")), patient))
    record(path, "build design", t_build)
    Wuse <- Di; absorb <- NULL; lam <- c(rep(0, px), 1e-3)
    st <- c(rep(FALSE, px), TRUE)
    A0use <- cbind(matrix(0, ng, px), A0[, px + (seq_len(S) - 1L) * q + 1L, drop = FALSE])
    pen_i <- ns$.blockExpand(Di, lam, "pen")
    solver_mk <- function() ns$.newtonSolverCompact(Di, pen_i)
  } else if (path == "dense") {
    t_build <- tm(Wd <- as.matrix(D))
    record(path, "build design", t_build)
    Wuse <- Wd; absorb <- grp; lam <- pen_full; st <- start_full; A0use <- A0
    solver_mk <- function() ns$.newtonSolver(Wd, pen_full, grp)
  } else {
    record(path, "build design", t_design)
    Wuse <- D; absorb <- NULL; lam <- pen; st <- start; A0use <- A0
    solver_mk <- function() ns$.newtonSolverCompact(D, pen_full)
  }
  t_solver <- tm(sol <- solver_mk())
  record(path, "solver setup", t_solver)
  t_fac <- median(replicate(reps, tm(sol$factor(w))))
  record(path, "factor(w) [gram + Schur]", t_fac)
  fac <- sol$factor(w)
  sc <- rnorm(ncol(Wuse))
  t_solve <- median(replicate(reps, tm(sol$solve(fac, sc))))
  record(path, "solve(state, s)", t_solve)
  a <- A0use[1, ]
  r <- rnorm(n)
  t_eta <- median(replicate(reps, tm(ns$.designEta(Wuse, a))))
  record(path, "eta (one gene)", t_eta)
  t_score <- median(replicate(reps, tm(ns$.designCross(Wuse, r))))
  record(path, "t(W) r (one gene)", t_score)
  t_pol <- tm(fit <- polishNB(Y, Wuse, A0use, rep(1, ng), lambda.a = lam,
                              absorb = absorb, start.cols = st,
                              psi.method = "profile"))
  record(path, sprintf("polishNB, %d genes", ng), t_pol)
  record(path, "polishNB per gene", t_pol / ng)
  record(path, "Newton iterations per gene", mean(fit$polish$iterations))
  res[[path]] <- fit
  rm(Wuse, sol, fac)
  if (path == "dense") rm(Wd)
  invisible(gc())
}
if (all(c("intercept", "intercept-compact") %in% names(res))) {
  d <- res$intercept; cc <- res[["intercept-compact"]]
  # the same model in two layouts: [P | X] against [X | P]
  cat(sprintf("intercept agreement: max|d alpha_x| %.2e, max|d loglik|/|loglik| %.2e, max|d psi| %.2e\n",
              max(abs(d$alpha[, S + seq_len(px)] - cc$alpha[, seq_len(px)])),
              max(abs(d$loglik - cc$loglik) / abs(d$loglik)), max(abs(d$psi - cc$psi))))
}
if (all(c("dense", "compact") %in% names(res))) {
  d <- res$dense; cc <- res$compact
  cat(sprintf("agreement: max|d alpha_x| %.2e, max|d loglik|/|loglik| %.2e, max|d psi| %.2e, iterations %s vs %s\n",
              max(abs(d$alpha[, seq_len(px)] - cc$alpha[, seq_len(px)])),
              max(abs(d$loglik - cc$loglik) / abs(d$loglik)), max(abs(d$psi - cc$psi)),
              paste(d$polish$iterations, collapse = ","), paste(cc$polish$iterations, collapse = ",")))
}
out <- do.call(rbind, rows)
if (nzchar(opt$out)) {
  write.table(out, opt$out, sep = ",", row.names = FALSE,
              col.names = !file.exists(opt$out), append = file.exists(opt$out))
}
