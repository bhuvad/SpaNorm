library(testthat)
library(SpatialExperiment)

# Create helper function to generate test data
create_test_spe <- function(n_genes = 100, n_spots = 50) {
  counts <- matrix(rpois(n_genes * n_spots, lambda = 5), 
                  nrow = n_genes, 
                  ncol = n_spots)
  rownames(counts) <- paste0("gene", 1:n_genes)
  colnames(counts) <- paste0("spot", 1:n_spots)
  
  coords <- data.frame(
    x = runif(n_spots),
    y = runif(n_spots)
  )
  
  spe <- SpatialExperiment(
    assays = list(counts = counts),
    colData = coords,
    spatialCoordsNames = c("x", "y")
  )
  return(spe)
}

# Add helper to create a valid SpaNorm model
create_test_spanorm <- function(n_genes = 100, n_spots = 50, n_covariates = 4) {
  W <- matrix(rnorm(n_spots * n_covariates), n_spots, n_covariates)
  SpaNormFit(
    ngenes = n_genes,
    ncells = n_spots,
    gene.model = "nb",
    df.tps = 6L,
    sample.p = 0.25,
    lambda.a = c(0.0001, 0.0002),
    batch = NULL,
    W = W,
    alpha = matrix(rnorm(n_genes * n_covariates), n_genes, n_covariates),
    gmean = rep(1, n_genes),
    psi = rep(1, n_genes),
    wtype = c("ls", "ls", rep("biology", n_covariates-2)),
    loglik = rep(-1, n_spots),
    sampling = factor(rep("dispersion", n_spots))
  )
}

test_that("SpaNormSVG fails appropriately without SpaNorm model", {
  spe <- create_test_spe()
  expect_error(
    SpaNormSVG(spe),
    "SpaNorm model"
  )
  
  # Test with NULL model
  metadata(spe)$SpaNorm <- NULL
  expect_error(
    SpaNormSVG(spe),
    "SpaNorm model"
  )
})

test_that("SpaNormSVG warns when overwriting existing results", {
  spe <- create_test_spe()
  rowData(spe)$svg.F <- NA
  rowData(spe)$svg.p <- NA
  rowData(spe)$svg.fdr <- NA
  
  metadata(spe)$SpaNorm <- create_test_spanorm(nrow(spe), ncol(spe))
  
  expect_warning(
    SpaNormSVG(spe),
    "SVG results exist"
  )
  
  # Verify old results are removed
  result <- suppressWarnings(SpaNormSVG(spe))
  expect_false(all(is.na(rowData(result)[, c("svg.F", "svg.p", "svg.fdr")])))
})

test_that("fitSpaNormTechnical handles inputs correctly", {
  Y <- matrix(rpois(100 * 50, lambda = 5), nrow = 100, ncol = 50)
  fit.spanorm <- create_test_spanorm()
  
  result <- fitSpaNormTechnical(Y, fit.spanorm, message)
  
  expect_true(is(result, "SpaNormFit"))
  expect_equal(ncol(result$W), 2)  # Should only include technical covariates
  expect_true(all(result$wtype != "biology"))
  expect_true(all(result$wtype == "ls"))
})

test_that("svgTest validates input dimensions", {
  Y <- matrix(rpois(100 * 50, lambda = 5), nrow = 100, ncol = 50)
  fit.spanorm <- create_test_spanorm(100, 50, 4)
  
  # Test mismatched genes
  fit.technical <- create_test_spanorm(99, 50, 3)  # Wrong number of genes
  expect_error(
    svgTest(Y, fit.spanorm, fit.technical),
    "number of genes differ between SpaNorm fits and/or data"
  )
  
  # Test mismatched cells
  fit.technical <- create_test_spanorm(100, 49, 3)  # Wrong number of cells
  expect_error(
    svgTest(Y, fit.spanorm, fit.technical),
    "number of cells differ between SpaNorm fits and/or data"
  )
  
  # Test non-nested technical model
  fit.technical <- create_test_spanorm(100, 50, 4)  # Same size as full model
  expect_error(
    svgTest(Y, fit.spanorm, fit.technical),
    "technical model is not nested in the full model"
  )
})

test_that("svgTest returns expected output format", {
  Y <- matrix(rpois(100 * 50, lambda = 5), nrow = 100, ncol = 50)
  fit.spanorm <- create_test_spanorm(100, 50, 4)
  fit.technical <- create_test_spanorm(100, 50, 2)
  
  result <- svgTest(Y, fit.spanorm, fit.technical)
  
  # Check structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(Y))
  expect_equal(colnames(result), c("svg.F", "svg.p", "svg.fdr"))
  
  # Check value ranges
  expect_true(all(result$svg.p >= 0 & result$svg.p <= 1))
  expect_true(all(result$svg.fdr >= 0 & result$svg.fdr <= 1))
  
  # Check FDR adjustment
  expect_equal(result$svg.fdr, p.adjust(result$svg.p, method = "fdr"))
})

test_that("SpaNormSVG handles edge cases", {
  # Test with single gene
  spe <- create_test_spe(n_genes = 1, n_spots = 50)
  metadata(spe)$SpaNorm <- create_test_spanorm(1, 50)
  expect_no_error(SpaNormSVG(spe))
  
  # Test with zero expression
  spe <- create_test_spe(n_genes = 10, n_spots = 10)
  counts(spe)[1,] <- 0  # Set first gene to zero expression
  metadata(spe)$SpaNorm <- create_test_spanorm(10, 10)
  expect_no_error(SpaNormSVG(spe))
})

test_that("SpaNormSVG results are stable across runs", {
  set.seed(42)
  spe <- create_test_spe(n_genes = 50, n_spots = 30)
  metadata(spe)$SpaNorm <- create_test_spanorm(50, 30)
  
  # Run twice and compare results
  result1 <- SpaNormSVG(spe)
  result2 <- SpaNormSVG(spe)
  
  expect_equal(
    rowData(result1)[, c("svg.F", "svg.p", "svg.fdr")],
    rowData(result2)[, c("svg.F", "svg.p", "svg.fdr")]
  )
})

test_that("SpaNormSVG validates input types", {
  # Test with invalid SpatialExperiment object
  expect_error(SpaNormSVG(list()))
  
  # Test with missing counts assay
  spe <- create_test_spe()
  assay(spe, "counts") <- NULL
  metadata(spe)$SpaNorm <- create_test_spanorm(nrow(spe), ncol(spe))
  expect_error(SpaNormSVG(spe))
})

test_that("SpaNormSVG checks model compatibility", {
  spe <- create_test_spe()
  
  # Test with mismatched model dimensions
  metadata(spe)$SpaNorm <- create_test_spanorm(nrow(spe) + 1, ncol(spe))
  expect_error(SpaNormSVG(spe))
  
  # Test with invalid model type
  metadata(spe)$SpaNorm <- list(type = "invalid")
  expect_error(SpaNormSVG(spe))
})

test_that("topSVGs returns expected results", {
  # Create test data
  spe <- create_test_spe(n_genes = 100, n_spots = 50)
  metadata(spe)$SpaNorm <- create_test_spanorm(100, 50)
  
  # Run SpaNormSVG
  spe <- SpaNormSVG(spe)
  
  # Basic functionality
  result <- topSVGs(spe, n = 10, fdr = 0.05)
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 3)
  expect_equal(colnames(result), c("svg.F", "svg.p", "svg.fdr"))
  expect_true(nrow(result) <= 10)
  
  # Check ordering
  expect_true(all(diff(result$svg.fdr) >= 0))
})

test_that("topSVGs handles edge cases", {
  spe <- create_test_spe(n_genes = 100, n_spots = 50)
  metadata(spe)$SpaNorm <- create_test_spanorm(100, 50)
  spe <- SpaNormSVG(spe)
  
  # Test with n larger than number of significant genes
  result <- topSVGs(spe, n = 1000, fdr = 0.001)
  expect_true(nrow(result) < 1000)
  
  # Test with n = 1
  result <- topSVGs(spe, n = 1)
  expect_equal(nrow(result), 1)
  
  # Test with strict FDR threshold
  result <- topSVGs(spe, fdr = 1e-10)
  expect_true(nrow(result) >= 0)
})

test_that("topSVGs validates inputs correctly", {
  spe <- create_test_spe()
  
  # Test without running SpaNormSVG first
  expect_error(
    topSVGs(spe),
    "SVG results not found"
  )
  
  # Test invalid n
  metadata(spe)$SpaNorm <- create_test_spanorm(nrow(spe), ncol(spe))
  spe <- SpaNormSVG(spe)
  expect_error(topSVGs(spe, n = 0))
  expect_error(topSVGs(spe, n = -1))
  
  # Test invalid FDR thresholds
  expect_error(topSVGs(spe, fdr = -0.1))
  expect_error(topSVGs(spe, fdr = 1.1))
})

test_that("fitSpaNormTechnical penalises the library-size terms as the full fit does", {
  # regression: the null model used to receive the raw lambda.a while the full
  # model received lambda.a * ncol(Y), so the nested null was under-penalised
  set.seed(1712)
  ngenes = 30
  ncells = 60
  Y = matrix(rpois(ngenes * ncells, 10), ngenes, ncells)
  coords = cbind(runif(ncells), runif(ncells))
  LS = colSums(Y) / mean(colSums(Y))
  quiet = function(...) invisible(NULL)

  # record the penalty vector handed to the NB fitter
  pen = list()
  orig = fitNBGivenPsi
  local_mocked_bindings(fitNBGivenPsi = function(Ysub, Wsub, psi, lambda.a, ...) {
    pen[[length(pen) + 1]] <<- setNames(lambda.a, colnames(Wsub)[-1])
    orig(Ysub, Wsub, psi, lambda.a, ...)
  })

  fit.full = fitSpaNorm(Y, coords, sample.p = 1, gene.model = "nb", df.tps = 2,
    lambda.a = c(1e-4, 2e-4), batch = NULL, LS = LS, msgfun = quiet,
    maxit.psi = 1, backend = "cpu")
  pen.full = pen[[1]]
  pen = list()
  fit.tech = fitSpaNormTechnical(Y, fit.full, quiet, maxit.psi = 1, backend = "cpu")
  pen.tech = pen[[1]]

  # shared library-size columns (the first column, logLS, is never penalised)
  ls.cols = colnames(fit.full$W)[-1][fit.full$wtype[-1] == "ls"]
  expect_gt(length(ls.cols), 0)
  expect_true(all(ls.cols %in% names(pen.tech)))
  expect_equal(pen.tech[ls.cols], pen.full[ls.cols])
  expect_equal(unname(pen.tech[ls.cols]), rep(2e-4 * ncells, length(ls.cols)))

  # the null must be fitted to the same cells as the full model
  expect_error(fitSpaNormTechnical(Y[, -1], fit.full, quiet), "number of cells")
})

# --- polished SVG consistency (Task 8) ---------------------------------
#
# `svgTest()` compares two SpaNormFit objects by an LRT; they must be
# estimated the same way, or the comparison is not a comparison of nested
# models any more. Task 7's "an existing null is polished alike and stale
# SVG columns are dropped" (tests/testthat/test-polishSpaNorm.R) already
# covers polishSpaNorm() polishing a stored null and clearing stale SVG
# columns; it is not duplicated here.

test_that("SpaNormSVG polishes the null when the full fit is polished", {
  spe <- polishSpaNorm(.polish_spe(), null = FALSE, verbose = FALSE)
  out <- SpaNormSVG(spe, verbose = FALSE)
  nul <- S4Vectors::metadata(out)$SpaNormNull
  expect_true(isPolished(nul))
  expect_identical(nul@polish$settings$psi.method,
                   S4Vectors::metadata(out)$SpaNorm@polish$settings$psi.method)
  # the pre-polish null is kept, as polishSpaNorm() keeps SpaNormNullUnpolished
  expect_false(isPolished(S4Vectors::metadata(out)$SpaNormNullUnpolished))
})

test_that("svgTest refuses a mixed pair", {
  spe <- SpaNormSVG(.polish_spe(), verbose = FALSE)          # unpolished pair
  full <- S4Vectors::metadata(spe)$SpaNorm
  nul <- S4Vectors::metadata(spe)$SpaNormNull
  Y <- SummarizedExperiment::assay(spe, "counts")
  fullP <- .polishSpaNormFit(full, as.matrix(Y))
  expect_error(svgTest(Y, fullP, nul), "polished alike")
})

test_that("SpaNormSVG refuses a polished null next to an unpolished full fit", {
  # the mirror-image mixed pair: a stored null that is MORE converged than
  # the (unpolished) full fit. SpaNormSVG() must not silently re-polish or
  # un-polish either side; it stops and names both states.
  spe <- .polish_spe()
  full <- S4Vectors::metadata(spe)$SpaNorm
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  nulP <- .polishSpaNormFit(fitSpaNormTechnical(Y, full, message, backend = "cpu"), Y,
                            verbose = FALSE)
  S4Vectors::metadata(spe)$SpaNormNull <- nulP
  expect_error(SpaNormSVG(spe, verbose = FALSE), "polishSpaNorm")
})

test_that("the polished full fit's own objective is at least the embedded polished null's", {
  # This is the property polishing actually guarantees (ruling 1, Task 8):
  # svgTest() scores winsorised mu/psi, under which the raw LRT statistic
  # can be mildly negative even for a correctly polished pair (psi.method =
  # "fixed" keeps each model's own fitNB dispersion, and winsorisation is
  # applied independently to each). What polishing DOES guarantee is that
  # the full fit's own (unwinsorised, penalised) objective, evaluated at its
  # own optimum, is at least as large as that SAME objective evaluated at
  # any other point in its parameter space -- including the null's optimum,
  # embedded with every biology column at 0 and the shared library-size
  # coefficient held at the full fit's own value.
  spe <- .polish_spe()
  full <- S4Vectors::metadata(spe)$SpaNorm
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))

  fullP <- .polishSpaNormFit(full, Y, verbose = FALSE)
  nullP <- .polishSpaNormFit(fitSpaNormTechnical(Y, full, message, backend = "cpu"), Y,
                             verbose = FALSE)

  probFull <- .spaNormPolishProblem(fullP, "all")
  A.full <- cbind(fullP@gmean, fullP@alpha[, -1, drop = FALSE])
  colnames(A.full) <- colnames(probFull$X)

  probNull <- .spaNormPolishProblem(nullP, "all")
  A.null <- cbind(nullP@gmean, nullP@alpha[, -1, drop = FALSE])
  colnames(A.null) <- colnames(probNull$X)

  # the null's columns are a subset of the full design's, matched by name;
  # everything else (the biology columns) is 0
  common <- intersect(colnames(A.full), colnames(A.null))
  expect_true(length(common) > 0)
  expect_lt(length(common), ncol(A.full))          # there IS a biology column
  A.embed <- matrix(0, nrow(A.full), ncol(A.full), dimnames = dimnames(A.full))
  A.embed[, common] <- A.null[, common]

  # SAME objective throughout: the full design/offset/psi/penalty
  offsetMat <- .offsetRows(probFull$offset, seq_len(nrow(Y)), probFull$X)
  mu.full <- .muBatch(A.full, probFull$X, offsetMat)
  mu.embed <- .muBatch(A.embed, probFull$X, offsetMat)
  obj.full <- .nbLoglikBatch(Y, mu.full, fullP@psi, A.full, probFull$pen)
  obj.embed <- .nbLoglikBatch(Y, mu.embed, fullP@psi, A.embed, probFull$pen)

  gap <- (obj.embed - obj.full) / pmax(abs(obj.full), 1)
  expect_true(all(gap <= 1e-8))
})

test_that("SpaNormSVG's unpolished SVG columns are master's statistic on the same pair", {
  # Ruling 5 (Task 8): the pairing logic only changes behaviour when a
  # polished/unpolished mismatch exists, so on a plain unpolished pair the
  # statistic is master's. The oracle below is master b5d3f8f's svgTest()
  # body, verbatim from its mu/psi lines on, applied to the full and null
  # fits SpaNormSVG() used. A stored baseline cannot hold this across
  # platforms: the fits are rebuilt at test time, and a different BLAS
  # kernel, or only a different BLAS thread count, moves them in the last
  # bits (final review I-1).
  spe <- .polish_spe()
  out <- suppressWarnings(SpaNormSVG(spe, backend = "cpu", verbose = FALSE))
  fit.spanorm <- S4Vectors::metadata(out)$SpaNorm
  fit.technical <- S4Vectors::metadata(out)$SpaNormNull
  expect_false(isPolished(fit.spanorm))
  expect_false(isPolished(fit.technical))

  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  mu <- calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  psi <- winsorisePsi(fit.spanorm$psi)
  loglik.spanorm <- rowSums(dnbinom(Y, mu = mu, size = 1 / psi, log = TRUE))
  mu <- calculateMu(fit.technical$gmean, fit.technical$alpha, fit.technical$W)
  psi <- winsorisePsi(fit.technical$psi)
  loglik.technical <- rowSums(dnbinom(Y, mu = mu, size = 1 / psi, log = TRUE))
  df1 <- ncol(fit.spanorm$W) - ncol(fit.technical$W)
  df2 <- ncol(Y) - ncol(fit.spanorm$W)
  F.lrt <- 2 * (loglik.spanorm - loglik.technical) / df1
  F.lrt <- pmax(F.lrt, 0)
  p.val <- pf(F.lrt, df1, df2, lower.tail = FALSE)
  fdr <- p.adjust(p.val, method = "fdr")
  oracle <- data.frame(svg.F = F.lrt, svg.p = p.val, svg.fdr = fdr)

  res <- getSVGResults(out)
  expect_identical(rownames(res), rownames(Y))
  expect_equal(res, oracle[rownames(res), ], tolerance = 1e-10)
})

# --- fix round 1: two polishSpaNorm() calls can polish full and null with
# different SETTINGS, both isPolished() == TRUE. The original ruling 3 only
# checked isPolished() equality ("consistent"), so this pair was accepted
# silently -- reachable through the public API (the review's exact
# sequence below), not just by hand-assembling a mismatched pair. -------

test_that("SpaNormSVG repairs a settings-mismatched polished pair (reviewer's sequence)", {
  spe <- SpaNormSVG(.polish_spe(), verbose = FALSE)
  # full + null, psi.method = "fixed"; polishSpaNorm() warns because 'spe'
  # already carries SVG results from the SpaNormSVG() call above -- expected
  expect_warning(spe <- polishSpaNorm(spe, verbose = FALSE), "SVG")
  spe <- polishSpaNorm(spe, overwrite = TRUE, null = FALSE,
                       psi.method = "profile", verbose = FALSE)    # re-polishes ONLY the full fit

  expect_message(out <- SpaNormSVG(spe, verbose = FALSE), "differ")

  full <- S4Vectors::metadata(out)$SpaNorm
  nul <- S4Vectors::metadata(out)$SpaNormNull
  expect_true(isPolished(nul))
  fs <- .polishSlot(full)$settings
  ns <- .polishSlot(nul)$settings
  expect_identical(ns$psi.method, fs$psi.method)
  expect_identical(ns$ls, fs$ls)
  expect_identical(ns$cells, fs$cells)
  expect_identical(ns$psi.method, "profile")

  # a pair polished consistently from scratch (psi.method = "profile" on
  # both fits from the start) must give the SAME svg.F, to tolerance
  ref <- SpaNormSVG(.polish_spe(), verbose = FALSE)
  expect_warning(ref <- polishSpaNorm(ref, psi.method = "profile", verbose = FALSE), "SVG")
  ref <- SpaNormSVG(ref, verbose = FALSE)

  expect_equal(SummarizedExperiment::rowData(out)$svg.F,
              SummarizedExperiment::rowData(ref)$svg.F, tolerance = 1e-8)
})

test_that("svgTest refuses a hand-made settings-mismatched polished pair", {
  spe <- .polish_spe()
  full <- S4Vectors::metadata(spe)$SpaNorm
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  fullP <- .polishSpaNormFit(full, Y, psi.method = "fixed", verbose = FALSE)
  nullP <- .polishSpaNormFit(fitSpaNormTechnical(Y, full, message, backend = "cpu"), Y,
                             psi.method = "profile", verbose = FALSE)
  expect_error(svgTest(SummarizedExperiment::assay(spe, "counts"), fullP, nullP),
              "polished alike")
})

test_that("SpaNormSVG forwards maxit/tol (not only psi.method/ls/cells) to the polished null", {
  spe <- polishSpaNorm(.polish_spe(), null = FALSE, maxit = 7L, tol = 1e-5, verbose = FALSE)
  out <- SpaNormSVG(spe, verbose = FALSE)
  fs <- .polishSlot(S4Vectors::metadata(out)$SpaNorm)$settings
  ns <- .polishSlot(S4Vectors::metadata(out)$SpaNormNull)$settings
  expect_identical(fs$maxit, 7L)
  expect_identical(fs$tol, 1e-5)
  expect_identical(ns$maxit, fs$maxit)
  expect_identical(ns$tol, fs$tol)
})

# --- final review I-2: a batch design, with a gene absent from one level ---

test_that("SpaNormSVG pairs and runs on a batch design, holding the gene out of the null too", {
  g <- .polish_batch_gene
  expect_warning(spe <- polishSpaNorm(.polish_batch_spe(), null = FALSE, verbose = FALSE),
                 "batch level")
  # the null keeps the batch column, so the same gene is held out of its polish
  expect_warning(out <- SpaNormSVG(spe, backend = "cpu", verbose = FALSE),
                 "^1 gene with no counts in some batch level.*'SpaNormNull'")
  full <- S4Vectors::metadata(out)$SpaNorm
  nul <- S4Vectors::metadata(out)$SpaNormNull
  expect_true("batch" %in% nul$wtype)
  expect_true(isPolished(nul))
  expect_true(.polishedAlike(full, nul))
  expect_false(isPolished(S4Vectors::metadata(out)$SpaNormNullUnpolished))
  np <- .polishSlot(nul)$genes
  expect_identical(np$held_out[g], "zero-batch-level")
  expect_false(np$polished[g])
  expect_true(all(np$polished[-g]))
  res <- getSVGResults(out)
  expect_identical(nrow(res), nrow(spe))
  expect_true(all(is.finite(res$svg.F)))
})
