# One small real SpaNorm fit, built once per test run (~15 s).
# suppressWarnings(): scran::quickCluster() warns about irlba on a 600-spot
# subset; that is the library-size step, not the fit under test.
.polish_spe <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      data(HumanDLPFC, package = "SpaNorm", envir = environment())
      set.seed(20260925)
      top <- order(-Matrix::rowSums(SummarizedExperiment::assay(HumanDLPFC, "counts")))[1:40]
      spe <- HumanDLPFC[top, 1:600]
      cache <<- suppressWarnings(SpaNorm(spe, df.tps = 2L, sample.p = 0.25,
                                         verbose = FALSE, backend = "cpu"))
    }
    cache
  }
})

# A small real SpaNorm fit with a two-level batch factor (one unpenalised 0/1
# batch column), on which gene .polish_batch_gene has no counts anywhere in
# level "s2" (colData column 'batch'). Built once per test run.
.polish_batch_gene <- 5L
.polish_batch_spe <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      data(HumanDLPFC, package = "SpaNorm", envir = environment())
      set.seed(20260926)
      top <- order(-Matrix::rowSums(SummarizedExperiment::assay(HumanDLPFC, "counts")))[1:30]
      spe <- HumanDLPFC[top, 1:500]
      b <- factor(rep(c("s1", "s2"), each = ncol(spe) / 2))
      cnt <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
      cnt[.polish_batch_gene, b == "s2"] <- 0
      SummarizedExperiment::assay(spe, "counts") <- cnt
      spe$batch <- b
      cache <<- suppressWarnings(SpaNorm(spe, df.tps = 2L, sample.p = 0.5, batch = b,
                                         verbose = FALSE, backend = "cpu"))
    }
    cache
  }
})

# The largest absolute penalised score of each gene at a SpaNorm-model fit,
# relative to the gene's total count, over the cells `cells` selects: zero at
# each gene's own optimum with the shared library-size coefficient held.
.polish_scaled_score <- function(fit, Y, cells = "all") {
  p <- .spaNormPolishProblem(fit, cells)
  Y <- Y[, p$cells_idx, drop = FALSE]
  A <- cbind(fit$gmean, fit$alpha[, -1])
  vapply(seq_len(nrow(Y)), function(g) {
    mu <- pmax(exp(as.numeric(p$X %*% A[g, ] + p$offset)), nbMuFloor())
    sc <- crossprod(p$X, (Y[g, ] - mu) / (1 + fit$psi[g] * mu)) - p$pen * A[g, ]
    max(abs(sc)) / sum(Y[g, ])
  }, numeric(1))
}
