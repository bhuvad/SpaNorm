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
