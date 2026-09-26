test_that("the polished fit is at each gene's optimum with a1 held", {
  # tol is the relative log-likelihood gain at which Newton stops. At the
  # default 1e-8 the score is left at up to 2.1e-6 * sum(y) on this fixture
  # (2 of 40 genes above 1e-6; one more exact Newton step gains 4e-8 in
  # log-likelihood and moves the intercept by 3e-6), so a zero-score check at
  # 1e-6 needs the tighter stop, as in test-polishNB.R (it gives 5.7e-9)
  spe <- polishSpaNorm(.polish_spe(), tol = 1e-12, verbose = FALSE)
  f <- S4Vectors::metadata(spe)$SpaNorm
  f0 <- S4Vectors::metadata(spe)$SpaNormUnpolished
  expect_true(isPolished(f)); expect_false(isPolished(f0))
  expect_equal(f$alpha[, 1], f0$alpha[, 1])                  # a1 untouched
  p <- .spaNormPolishProblem(f)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  A <- cbind(f$gmean, f$alpha[, -1])
  for (g in seq_len(nrow(Y))) {
    mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
    sc <- crossprod(p$X, (Y[g, ] - mu) / (1 + f$psi[g] * mu)) - p$pen * A[g, ]
    expect_lt(max(abs(sc)) / sum(Y[g, ]), 1e-6, label = rownames(Y)[g])
  }
  expect_equal(f$psi, f0$psi)                                # psi.method = "fixed"
})

test_that("polishing raises every gene's penalised likelihood", {
  spe <- .polish_spe()
  f0 <- S4Vectors::metadata(spe)$SpaNorm
  f <- S4Vectors::metadata(polishSpaNorm(spe, verbose = FALSE))$SpaNorm
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  pl <- function(fit) {
    p <- .spaNormPolishProblem(fit)
    A <- cbind(fit$gmean, fit$alpha[, -1])
    vapply(seq_len(nrow(Y)), function(g) {
      mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
      sum(dnbinom(Y[g, ], size = 1 / fit$psi[g], mu = mu, log = TRUE)) -
        0.5 * sum(p$pen * A[g, ]^2)
    }, numeric(1))
  }
  expect_true(all(pl(f) >= pl(f0) - 1e-8 * abs(pl(f0))))
})

test_that("logcounts are rewritten from the polished fit", {
  spe <- .polish_spe()
  out <- polishSpaNorm(spe, verbose = FALSE)
  expect_false(identical(SummarizedExperiment::assay(out, "logcounts"),
                         SummarizedExperiment::assay(spe, "logcounts")))
  # the same result as SpaNorm() re-normalising with the polished fit cached
  again <- suppressWarnings(SpaNorm(out, df.tps = 2L, sample.p = 0.25, verbose = FALSE))
  expect_equal(SummarizedExperiment::assay(again, "logcounts"),
               SummarizedExperiment::assay(out, "logcounts"))
})

test_that("a non-integer assay is refused at the entry point", {
  spe <- .polish_spe()
  SummarizedExperiment::assay(spe, "counts") <- SummarizedExperiment::assay(spe, "counts") + 0.5
  expect_error(polishSpaNorm(spe, verbose = FALSE), "integer counts")
})

test_that("an all-zero gene is not polished and keeps its input fit", {
  # The counts are edited after the fit, which is valid: the polish reads the
  # counts, and the fit's coefficients for that gene are just a bad start.
  spe <- .polish_spe()
  cnt <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  cnt[1, ] <- 0
  SummarizedExperiment::assay(spe, "counts") <- cnt
  # held out quietly (a verbose message only), as before the batch-level rule
  expect_no_warning(out <- polishSpaNorm(spe, verbose = FALSE))
  f <- S4Vectors::metadata(out)$SpaNorm
  f0 <- S4Vectors::metadata(out)$SpaNormUnpolished
  pol <- .polishSlot(f)$genes
  expect_false(pol$polished[1])
  expect_true(all(pol$polished[-1]))
  expect_identical(pol$held_out[1], "all-zero")
  expect_true(all(is.na(pol$held_out[-1])))
  expect_identical(pol$iterations[1], 0L)
  expect_true(is.na(pol$loglik[1]))
  expect_identical(f$gmean[1], f0$gmean[1])
  expect_identical(f$alpha[1, ], f0$alpha[1, ])
  expect_identical(f$psi[1], f0$psi[1])
  expect_true(all(is.finite(f$gmean[-1])))
  expect_true(all(is.finite(f$alpha[-1, ])))
  expect_true(all(is.finite(f$psi[-1])))
  # the other genes are polished exactly as without the zero gene in the set
  ref <- S4Vectors::metadata(polishSpaNorm(.polish_spe(), verbose = FALSE))$SpaNorm
  expect_equal(f$gmean[-1], ref$gmean[-1], tolerance = 1e-10)
})

test_that("a fit whose library-size coefficient is not shared is refused by name", {
  # a fitNB()-style or hand-built fit with a per-gene alpha[, 1] is not a
  # SpaNorm-model fit: the polish holds that column as one shared offset
  spe <- .polish_spe()
  f <- S4Vectors::metadata(spe)$SpaNorm
  f@alpha[, 1] <- f@alpha[, 1] + seq(0, 0.1, length.out = nrow(f@alpha))
  S4Vectors::metadata(spe)$SpaNorm <- f
  before <- spe
  expect_error(polishSpaNorm(spe, verbose = FALSE), "column 1 .*not shared")
  # nothing was written to the object before the refusal
  expect_identical(spe, before)
  expect_null(S4Vectors::metadata(spe)$SpaNormUnpolished)
  expect_identical(S4Vectors::metadata(spe)$SpaNorm, f)
})

test_that("the polish records its settings and per-gene diagnostics", {
  spe <- .polish_spe()
  f0 <- S4Vectors::metadata(spe)$SpaNorm
  f <- S4Vectors::metadata(polishSpaNorm(spe, verbose = FALSE))$SpaNorm
  s <- .polishSlot(f)$settings
  expect_true(all(c("psi.method", "ls", "cells", "maxit", "tol", "pen",
                    "a1.input", "a1", "SpaNorm") %in% names(s)))
  expect_identical(s[c("psi.method", "ls", "cells")],
                   list(psi.method = "fixed", ls = "fixed", cells = "all"))
  expect_equal(s$pen, .spaNormPenalty(f0))
  expect_identical(s$a1, f0$alpha[1, 1])
  expect_identical(s$a1.input, s$a1)                         # ls = "fixed"
  genes <- .polishSlot(f)$genes
  expect_identical(nrow(genes), nrow(spe))
  expect_identical(rownames(genes), rownames(spe))
  expect_true(all(c("iterations", "psi_fitnb", "restarted", "capped",
                    "singular", "psi_bound", "polished", "loglik") %in% names(genes)))
  expect_true(all(is.finite(genes$loglik)))
  expect_true(methods::validObject(f))
  # the fit's own iteration trace is left alone
  expect_identical(f$loglik, f0$loglik)
})

test_that("show() of a polished fit prints the settings it was polished with", {
  spe <- .polish_spe()
  out <- polishSpaNorm(spe, psi.method = "profile", cells = "fit", tol = 1e-12,
                       verbose = FALSE)
  f <- S4Vectors::metadata(out)$SpaNorm
  expect_true(any(capture.output(show(f)) ==
                    "polished: psi.method = profile, ls = fixed, cells = fit"))
  # cells = "fit" converges each gene over the fit's cells, at its own psi (the
  # tighter tol: see the first test)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  expect_lt(max(.polish_scaled_score(f, Y, cells = "fit")), 1e-6)
})

test_that("an already polished fit is refused unless overwrite = TRUE", {
  spe <- .polish_spe()
  f_in <- S4Vectors::metadata(spe)$SpaNorm
  once <- polishSpaNorm(spe, psi.method = "profile", verbose = FALSE)
  expect_error(polishSpaNorm(once, verbose = FALSE), "already polished.*overwrite")

  twice <- polishSpaNorm(once, overwrite = TRUE, verbose = FALSE)
  f1 <- S4Vectors::metadata(once)$SpaNorm
  f2 <- S4Vectors::metadata(twice)$SpaNorm
  # SpaNormUnpolished still holds the original fit, never a polished one
  expect_identical(S4Vectors::metadata(twice)$SpaNormUnpolished, f_in)
  expect_true(isPolished(f2))
  expect_identical(.polishSlot(f2)$settings$psi.method, "fixed")
  # the second polish starts from the original: "fixed" keeps the ORIGINAL
  # dispersion, not the one the first (profile) polish estimated
  expect_false(isTRUE(all.equal(f1$psi, f_in$psi)))
  expect_identical(f2$psi, f_in$psi)
  expect_identical(.polishSlot(f2)$genes$psi_fitnb, f_in$psi)
  # and the result is the one-shot polish of the original
  direct <- S4Vectors::metadata(polishSpaNorm(spe, verbose = FALSE))$SpaNorm
  expect_identical(f2, direct)
})

test_that("an existing null is polished alike and stale SVG columns are dropped", {
  spe <- suppressWarnings(SpaNormSVG(.polish_spe(), backend = "cpu", verbose = FALSE))
  n0 <- S4Vectors::metadata(spe)$SpaNormNull
  expect_false(isPolished(n0))
  expect_warning(out <- polishSpaNorm(spe, tol = 1e-12, verbose = FALSE), "SVG")
  md <- S4Vectors::metadata(out)
  expect_true(isPolished(md$SpaNormNull))
  expect_identical(md$SpaNormNullUnpolished, n0)
  expect_identical(.polishSlot(md$SpaNormNull)$settings[c("psi.method", "ls", "cells")],
                   .polishSlot(md$SpaNorm)$settings[c("psi.method", "ls", "cells")])
  expect_false(any(c("svg.F", "svg.p", "svg.fdr") %in%
                     colnames(SummarizedExperiment::rowData(out))))
  # the null is at its own optimum too (the tighter tol: see the first test)
  Y <- as.matrix(SummarizedExperiment::assay(out, "counts"))
  expect_lt(max(.polish_scaled_score(md$SpaNormNull, Y)), 1e-6)

  # null = FALSE leaves a stored null exactly as it was
  expect_warning(out2 <- polishSpaNorm(spe, null = FALSE, verbose = FALSE), "SVG")
  expect_identical(S4Vectors::metadata(out2)$SpaNormNull, n0)
  expect_null(S4Vectors::metadata(out2)$SpaNormNullUnpolished)
})

test_that(".polishSpaNormFit() polishes a fit saved before the polish slot existed", {
  old <- readRDS(testthat::test_path("fixtures", "spanormfit_1713_noslot.rds"))
  expect_false(methods::.hasSlot(old, "polish"))
  # counts from the fixture's own model, so the polish problem is well posed
  set.seed(71)
  eta <- old$gmean + old$alpha %*% t(old$W)
  Y <- matrix(rnbinom(length(eta), mu = exp(eta), size = 1 / old$psi),
              old$ngenes, old$ncells)
  f <- .polishSpaNormFit(old, Y, tol = 1e-12)   # the tighter tol: see the first test
  expect_true(isPolished(f))
  expect_true(methods::validObject(f))
  expect_true(all(.polishSlot(f)$genes$polished))
  expect_lt(max(.polish_scaled_score(f, Y)), 1e-6)
})

test_that("a DelayedArray counts assay polishes and normalises as in memory", {
  skip_if_not_installed("DelayedArray")
  spe <- .polish_spe()
  da <- spe
  SummarizedExperiment::assay(da, "counts") <-
    DelayedArray::DelayedArray(as.matrix(SummarizedExperiment::assay(spe, "counts")))
  old <- DelayedArray::getAutoBlockSize()
  on.exit(suppressMessages(DelayedArray::setAutoBlockSize(old)))
  suppressMessages(DelayedArray::setAutoBlockSize(ncol(spe) * 8 * 7)) # ~7 genes/block
  a <- polishSpaNorm(da, verbose = FALSE)
  b <- polishSpaNorm(spe, verbose = FALSE)
  expect_identical(S4Vectors::metadata(a)$SpaNorm, S4Vectors::metadata(b)$SpaNorm)
  expect_equal(as.matrix(SummarizedExperiment::assay(a, "logcounts")),
               as.matrix(SummarizedExperiment::assay(b, "logcounts")), tolerance = 0)
})

test_that("the Seurat method polishes the fit in @misc and rewrites the data layer", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  spe <- .polish_spe()
  counts <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  xy <- SpatialExperiment::spatialCoords(spe)
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts, assay = "Spatial"))
  df <- data.frame(x = xy[, 1], y = xy[, 2], cell = colnames(counts))
  fov <- suppressWarnings(
    SeuratObject::CreateFOV(df, type = "centroids", key = "fov_", assay = "Spatial"))
  suppressWarnings(obj[["slice1"]] <- fov)
  # the SPE fixture's fit, so the two methods polish the same input
  f_in <- S4Vectors::metadata(spe)$SpaNorm
  obj@misc$SpaNorm <- f_in

  out <- polishSpaNorm(obj, assay = "Spatial", verbose = FALSE)
  ref <- polishSpaNorm(spe, verbose = FALSE)
  expect_true(isPolished(out@misc$SpaNorm))
  expect_identical(out@misc$SpaNormUnpolished, f_in)
  expect_identical(out@misc$SpaNorm, S4Vectors::metadata(ref)$SpaNorm)
  dl <- SeuratObject::LayerData(out, layer = "data", assay = "Spatial")
  expect_equal(as.matrix(dl),
               as.matrix(SummarizedExperiment::assay(ref, "logcounts")),
               tolerance = 0, ignore_attr = TRUE)
  expect_error(polishSpaNorm(out, assay = "Spatial", verbose = FALSE),
               "already polished.*overwrite")

  # a non-integer counts layer is refused here too
  bad <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts + 0.5,
                                                           assay = "Spatial"))
  suppressWarnings(bad[["slice1"]] <- fov)
  bad@misc$SpaNorm <- f_in
  expect_error(polishSpaNorm(bad, assay = "Spatial", verbose = FALSE), "integer counts")
})

# ---- final review I-2: a gene with no counts in a whole batch level -------
#
# Batch columns are unpenalised 0/1 indicators, so a gene with no counts in
# every cell of one level has no finite optimum: its level coefficient runs
# to -Inf, exactly the degeneracy an all-zero gene has in its intercept.
# Measured before the fix: the coefficient went -0.136 -> -20.9, reported
# polished, and the gene's logcounts in that level (every raw count zero)
# went from a mean 0.67 to 6.04.

test_that("hold-out groups come from the unpenalised 0/1 columns only", {
  X <- cbind(`(gmean)` = 1, bio = c(-1.2, 0.3, 0.8, -0.4, 1.1, -0.6),
             b2 = c(0, 0, 1, 1, 0, 0), b3 = c(0, 0, 0, 0, 1, 1))
  pen <- c(0, 1, 0, 0)
  # the unique rows of the unpenalised block: one group per batch level
  expect_identical(.unpenalisedGroups(X, pen), c(1L, 1L, 2L, 2L, 3L, 3L))
  # nothing unpenalised but the intercept: one group, the all-zero rule
  expect_identical(.unpenalisedGroups(X, c(0, 1, 1, 1)), rep(1L, 6))
  # an unpenalised column that is not 0/1 defines no groups
  expect_null(.unpenalisedGroups(X, c(0, 0, 0, 0)))

  Y <- rbind(c(1, 2, 3, 4, 5, 6),     # counts in every level
             c(0, 0, 0, 0, 0, 0),     # all-zero
             c(1, 1, 0, 0, 2, 0),     # none in level 2
             c(0, 0, 1, 0, 0, 1))     # none in level 1 (the reference)
  expect_identical(.polishHoldOut(Y, X, pen),
                   c(NA, "all-zero", "zero-batch-level", "zero-batch-level"))
  expect_identical(.polishHoldOut(Matrix::Matrix(Y, sparse = TRUE), X, pen),
                   .polishHoldOut(Y, X, pen))
  # no groups: the all-zero rule alone
  expect_identical(.polishHoldOut(Y, X, c(0, 0, 0, 0)), c(NA, "all-zero", NA, NA))
  expect_identical(.polishHoldOut(Y, X, c(0, 1, 1, 1)), c(NA, "all-zero", NA, NA))
})

test_that("a gene with no counts in a whole batch level is held out, with a warning", {
  spe <- .polish_batch_spe()
  g <- .polish_batch_gene
  lvl <- spe$batch == "s2"
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  expect_identical(sum(Y[g, lvl]), 0)                 # the premise
  expect_gt(sum(Y[g, !lvl]), 0)
  f0 <- S4Vectors::metadata(spe)$SpaNorm
  expect_true("batch" %in% f0$wtype)

  expect_warning(out <- polishSpaNorm(spe, verbose = FALSE),
                 "^1 gene with no counts in some batch level")
  f <- S4Vectors::metadata(out)$SpaNorm
  pol <- .polishSlot(f)$genes
  expect_true(isPolished(f))
  expect_false(pol$polished[g])
  expect_identical(pol$held_out[g], "zero-batch-level")
  expect_identical(pol$iterations[g], 0L)
  expect_true(is.na(pol$loglik[g]))
  # its input fit is kept bit for bit, the batch coefficient included
  expect_identical(f$gmean[g], f0$gmean[g])
  expect_identical(f$alpha[g, ], f0$alpha[g, ])
  expect_identical(f$psi[g], f0$psi[g])

  # the normalised values in the level it was never seen in are the unpolished
  # ones, not imputed at full strength (before the fix: mean 6.11 there,
  # against 6.24 in the level it is expressed in; now 1.54)
  l0 <- as.matrix(SummarizedExperiment::assay(spe, "logcounts"))[g, ]
  l1 <- as.matrix(SummarizedExperiment::assay(out, "logcounts"))[g, ]
  expect_equal(l1, l0)
  expect_lt(mean(l1[lvl]), mean(l1[!lvl]) - 2)
})

test_that("a gene with counts in every batch level is polished, at its optimum", {
  spe <- .polish_batch_spe()
  g <- .polish_batch_gene
  expect_warning(out <- polishSpaNorm(spe, tol = 1e-12, verbose = FALSE), "batch level")
  f <- S4Vectors::metadata(out)$SpaNorm
  pol <- .polishSlot(f)$genes
  expect_true(all(pol$polished[-g]))
  expect_true(all(is.na(pol$held_out[-g])))
  # the batch column is unpenalised, as SpaNorm() fits it
  s <- .polishSlot(f)$settings
  expect_identical(s$pen[c(FALSE, f$wtype[-1] == "batch")], 0)
  # the tighter tol: see the first test
  Y <- as.matrix(SummarizedExperiment::assay(out, "counts"))
  expect_lt(max(.polish_scaled_score(f, Y)[-g]), 1e-6)
})
