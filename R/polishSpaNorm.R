#' Converge each gene of a SpaNorm fit to its own optimum
#'
#' [SpaNorm()] fits every gene in one shared IRLS loop, with a gene-averaged
#' cell weight vector and an aggregate convergence criterion, so a gene can be
#' left short of its own optimum (bright, spatially restricted genes most of
#' all). `polishSpaNorm()` takes the stored fit and converges each gene
#' separately under SpaNorm's own model, stores the polished fit, and rewrites
#' the normalised assay from it.
#'
#' @param spe a SpatialExperiment or Seurat object that [SpaNorm()] has been
#'   run on, with the raw counts in its 'counts' assay (layer, for Seurat).
#' @param adj.method a character, specifying the method used to rewrite the
#'   normalised data from the polished fit (default 'auto'), as in
#'   [SpaNorm()].
#' @param scale.factor a numeric, specifying the scaling factor for the
#'   adjusted counts, as in [SpaNorm()].
#' @param psi.method how each gene's dispersion is set: `"fixed"` (the
#'   default) keeps the fit's dispersion and converges the mean under it;
#'   `"profile"` re-estimates it by profile maximum likelihood at the
#'   converged mean and re-polishes (see [polishNB()]).
#' @param ls the library-size coefficient shared by every gene: `"fixed"` (the
#'   default) holds it at the fit's value; `"joint"` re-estimates it jointly
#'   with the per-gene coefficients, at the optimum of the total penalised
#'   likelihood over the polished genes (see Details).
#' @param cells the cells/spots each gene is converged over: `"all"` (the
#'   default), or `"fit"` for the ones [SpaNorm()] sampled to fit the model.
#' @param null a logical, specifying whether to also polish the null model
#'   stored by [SpaNormSVG()] ('SpaNormNull'), when there is one (default
#'   TRUE). The SVG test compares the two fits, so they must be estimated the
#'   same way.
#' @param maxit,tol the Newton iteration cap and the relative log-likelihood
#'   tolerance, per gene.
#' @param engine,batch.size,backend passed to [polishNB()]: the batched or
#'   per-gene Newton engine, genes per batch, and the compute backend.
#' @param BPPARAM a BiocParallelParam object over gene blocks, for the polish
#'   (see [polishNB()]) and the normalisation (default
#'   \code{BiocParallel::SerialParam()}).
#' @param overwrite a logical. A fit that is already polished is refused
#'   unless this is TRUE (default FALSE); with TRUE it is polished again,
#'   starting from the unpolished fit kept as 'SpaNormUnpolished' (and
#'   'SpaNormNullUnpolished'), never from the polished one.
#' @param verbose a logical, specifying whether to show update messages
#'   (default TRUE).
#' @param assay a character, specifying the assay holding the raw counts for
#'   Seurat objects (default NULL uses the object's default assay). Ignored
#'   for SpatialExperiment objects.
#' @param ... further arguments passed to [polishNB()], such as `block.size`
#'   or `psi.range`.
#'
#' @details SpaNorm's mean for gene \eqn{g} in cell/spot \eqn{i} is
#'   \deqn{\log \mu_{gi} = \bar{\mu}_g + a_1 w_{i1} + \sum_{j \ge 2} W_{ij}
#'   \alpha_{gj},}{log mu_gi = gmean_g + a1 * W_i1 + sum_{j >= 2} W_ij
#'   alpha_gj,} where \eqn{w_{i1}} is the log library size and its
#'   coefficient \eqn{a_1} is shared by every gene. The polish maximises each
#'   gene's penalised negative binomial log-likelihood by damped Newton (see
#'   [polishNB()]), with \eqn{a_1} held at the fit's value (`ls = "fixed"`),
#'   \eqn{\bar{\mu}_g}{gmean_g} a per-gene unpenalised intercept, and the
#'   other coefficients ridge-penalised as [SpaNorm()] penalises them: by
#'   `lambda.a` for the biology and library-size splines, not at all for
#'   batch, times the number of cells/spots.
#'
#'   The polished coefficients are the optimum of the unwinsorised penalised
#'   likelihood. The normalisation (and [SpaNormSVG()]) still applies
#'   SpaNorm's usual winsorisation, exactly as for an unpolished fit: it forms
#'   the mean through [calculateMu()], which by default winsorises the upper
#'   tail of each gene's log mean, and caps the dispersions at
#'   \code{exp(median(log psi) + 4 MAD(log psi))}.
#'
#'   A gene with no counts in the polished cells/spots has no finite optimum
#'   (its intercept runs to minus infinity), so it is not polished: it keeps
#'   its input coefficients and dispersion, and its diagnostics row reads
#'   `polished = FALSE` with `iterations = 0`. A gene the Newton engine
#'   cannot converge from either start also keeps its input fit (see
#'   [polishNB()]).
#'
#'   With `ls = "joint"` the shared coefficient \eqn{a_1} is moved to the
#'   joint optimum of the total penalised log-likelihood over the polished
#'   genes: after the per-gene polish at the fit's \eqn{a_1}, Newton steps on
#'   the profile log-likelihood of \eqn{a_1} (each gene's own coefficients
#'   profiled out; the information is the Fisher information of the polish)
#'   alternate with a warm re-polish of every gene at the candidate value,
#'   and a step is halved until the total does not fall. It stops when the
#'   standardised score \eqn{|U|/\sqrt{I}} is below `1e-6`, or after 10
#'   steps; `maxit` and `tol` are the per-gene Newton's, in the cold pass and
#'   every re-polish alike. With `psi.method = "profile"` each re-polish also
#'   re-profiles the dispersion. A gene whose information is singular at a
#'   step is left out of that step's score, information and total. If the
#'   profiled information is not positive (\eqn{a_1} not identified: the log
#'   library size is collinear with the model's other terms), the fit is the
#'   `ls = "fixed"` polish, with a warning. The joint estimate weights genes by their information, so
#'   bright genes dominate it, where [SpaNorm()]'s \eqn{a_1} is an unweighted
#'   mean over genes. Genes that are not polished (no counts, or not
#'   convergeable) do not inform \eqn{a_1}; they keep their input
#'   \eqn{\bar{\mu}_g}{gmean_g}, other coefficients and dispersion exactly,
#'   but take the new shared \eqn{a_1}, since it is one value for every
#'   gene, so their fitted mean shifts by \eqn{(a_1^{new} - a_1^{old})
#'   w_{i1}}{(a1_new - a1_old) * W_i1}. With no polished gene at all,
#'   \eqn{a_1} is left at the fit's value with a warning.
#'
#'   The polished fit replaces 'SpaNorm' in the object's metadata (`@misc`
#'   for Seurat) and the input fit is kept as 'SpaNormUnpolished'. With
#'   `null = TRUE` a stored 'SpaNormNull' is polished with the same settings
#'   and the input kept as 'SpaNormNullUnpolished'. SVG results in `rowData`
#'   were computed from the unpolished fit, so they are removed with a
#'   warning; rerun [SpaNormSVG()]. The normalised assay ('logcounts', or the
#'   'data' layer for Seurat) is rewritten from the polished fit.
#'
#'   The polished fit's `polish` slot (see [isPolished()]) holds `settings`
#'   (`psi.method`, `ls`, `cells`, `maxit`, `tol`, the penalty vector `pen`,
#'   the library-size coefficient before (`a1.input`) and after (`a1`), and
#'   the SpaNorm version) and `genes`, one row per gene: [polishNB()]'s
#'   diagnostics plus `loglik`, the penalised log-likelihood over the polished
#'   cells at the returned fit (`NA` for a gene that was not polished). The
#'   fit's `loglik` slot is left as the shared fit's iteration trace. With
#'   `ls = "joint"`, `settings` also holds `ls.iterations` (accepted steps on
#'   \eqn{a_1}), `ls.maxit` and `ls.tol` (the cap on those steps and the
#'   standardised-score stop, 10 and `1e-6`), `ls.score` (the final pooled
#'   score \eqn{U}), `ls.se`
#'   (\eqn{1/\sqrt{I}}, the profiled standard error of \eqn{a_1}),
#'   `ls.singular` (genes left out of \eqn{U} and \eqn{I} at some step
#'   because their information was singular) and `ls.converged`, and the
#'   per-gene `iterations`, `capped` and `singular` include the warm
#'   re-polishes.
#'
#'   The negative binomial likelihood is defined on counts only, so a
#'   non-integer counts assay (for example a back-transform such as
#'   `2^logcounts - 1`) is refused.
#'
#' @return a SpatialExperiment or Seurat object holding the polished fit(s),
#'   with the normalised data rewritten in 'logcounts' or 'data',
#'   respectively.
#' @seealso [SpaNorm()], [polishNB()], [isPolished()].
#' @name polishSpaNorm
#'
#' @examples
#' data(HumanDLPFC)
#' \donttest{
#' top <- order(-Matrix::rowSums(SummarizedExperiment::assay(HumanDLPFC, "counts")))[1:50]
#' spe <- SpaNorm(HumanDLPFC[top, ], sample.p = 0.05, df.tps = 2, tol = 1e-2)
#' spe <- polishSpaNorm(spe)
#' isPolished(S4Vectors::metadata(spe)$SpaNorm)
#' }
#' @export
setGeneric("polishSpaNorm", function(
    spe,
    adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
    scale.factor = 1,
    psi.method = c("fixed", "profile"),
    ls = c("fixed", "joint"),
    cells = c("all", "fit"),
    null = TRUE,
    maxit = 50L,
    tol = 1e-8,
    engine = c("batch", "gene"),
    batch.size = NULL,
    backend = c("cpu", "auto", "gpu"),
    BPPARAM = BiocParallel::SerialParam(),
    overwrite = FALSE,
    verbose = TRUE,
    assay = NULL,
    ...) {
  standardGeneric("polishSpaNorm")
})

#' @rdname polishSpaNorm
setMethod(
  "polishSpaNorm",
  signature("SpatialExperiment"),
  function(spe, adj.method, scale.factor, psi.method, ls, cells, null, maxit,
           tol, engine, batch.size, backend, BPPARAM, overwrite, verbose,
           assay, ...) {
    .polishCountsCheck(checkSPE(spe))
    adj.method <- match.arg(adj.method)
    psi.method <- match.arg(psi.method)
    ls <- match.arg(ls)
    cells <- match.arg(cells)
    engine <- match.arg(engine)
    backend <- match.arg(backend)

    emat <- SummarizedExperiment::assay(spe, "counts")
    md <- S4Vectors::metadata(spe)
    res <- .polishSpaNormCore(
      .polishStoredFits(spe, md, null), emat, adj.method, scale.factor,
      psi.method, ls, cells, null, maxit, tol, overwrite, verbose,
      engine = engine, batch.size = batch.size, backend = backend,
      BPPARAM = BPPARAM, ...
    )

    for (nm in names(res$fits)) S4Vectors::metadata(spe)[[nm]] <- res$fits[[nm]]
    # the SVG statistics compare the fits this call has just replaced
    svg <- getSVGResults(spe, stop_if_missing = FALSE)
    if (!is.null(svg)) {
      warning("SVG results in 'spe' were computed from the unpolished fit and ",
              "have been removed; rerun 'SpaNormSVG'", call. = FALSE)
      rd <- SummarizedExperiment::rowData(spe)
      SummarizedExperiment::rowData(spe) <- rd[, setdiff(colnames(rd), colnames(svg)), drop = FALSE]
    }
    SummarizedExperiment::assay(spe, "logcounts") <- methods::as(res$normmat, "sparseMatrix")

    spe
  }
)

#' @rdname polishSpaNorm
setMethod(
  "polishSpaNorm",
  signature("Seurat"),
  function(spe, adj.method, scale.factor, psi.method, ls, cells, null, maxit,
           tol, engine, batch.size, backend, BPPARAM, overwrite, verbose,
           assay, ...) {
    assay_name <- if (is.null(assay)) SeuratObject::DefaultAssay(spe) else assay
    .polishCountsCheck(checkSeurat(spe, assay = assay_name))
    adj.method <- match.arg(adj.method)
    psi.method <- match.arg(psi.method)
    ls <- match.arg(ls)
    cells <- match.arg(cells)
    engine <- match.arg(engine)
    backend <- match.arg(backend)

    emat <- SeuratObject::GetAssayData(spe, layer = "counts", assay = assay_name)
    res <- .polishSpaNormCore(
      .polishStoredFits(spe, spe@misc, null, assay = assay_name), emat,
      adj.method, scale.factor, psi.method, ls, cells, null, maxit, tol,
      overwrite, verbose,
      engine = engine, batch.size = batch.size, backend = backend,
      BPPARAM = BPPARAM, ...
    )

    for (nm in names(res$fits)) spe@misc[[nm]] <- res$fits[[nm]]
    spe <- SeuratObject::SetAssayData(spe, assay = assay_name, layer = "data",
                                      new.data = methods::as(res$normmat, "dgCMatrix"))

    spe
  }
)

# Run a container checker (checkSPE()/checkSeurat()), restating a non-integer
# counts error in the polish's terms: the generic checker says only that the
# values are non-integer, and the polish's reason, and the usual cause (the
# counts overwritten with a back-transform AFTER SpaNorm() ran on them), are
# what a user needs to fix it. Any other error passes through unchanged.
.polishCountsCheck <- function(check) {
  tryCatch(check, error = function(e) {
    if (grepl("non-integer", conditionMessage(e), fixed = TRUE)) {
      stop("polishSpaNorm() needs integer counts: the negative binomial ",
           "likelihood is undefined otherwise.\n  The counts assay is not ",
           "integer-valued (e.g. a back-transform such as 2^logcounts - 1 ",
           "written over the raw counts).\n  Use the raw counts.",
           call. = FALSE)
    }
    stop(e)
  })
}

# The fits a container holds, by the names polishSpaNorm() reads and writes.
# `store` is the SpatialExperiment's metadata or the Seurat object's @misc. The
# full fit must exist and match the data (getSpaNormFit() validates both); a
# null is validated only when it is going to be polished.
.polishStoredFits <- function(spe, store, null, assay = NULL) {
  nul <- NULL
  if (null && !is.null(getSpaNormFit(spe, null = TRUE, validate = FALSE, assay = assay))) {
    nul <- getSpaNormFit(spe, null = TRUE, assay = assay)
  }
  list(SpaNorm = getSpaNormFit(spe, assay = assay), SpaNormNull = nul,
       SpaNormUnpolished = store$SpaNormUnpolished,
       SpaNormNullUnpolished = store$SpaNormNullUnpolished)
}

# The fit to polish from. An unpolished fit is its own input. A polished one
# is refused without `overwrite`, and with it the stored unpolished fit is the
# input, so that '<name>Unpolished' always holds the unpolished fit.
.polishSource <- function(fit, unpolished, name, overwrite) {
  if (!isPolished(fit)) return(fit)
  if (!isTRUE(overwrite)) {
    stop(sprintf(paste0("the '%s' fit is already polished; pass overwrite = TRUE ",
                        "to polish it again from the unpolished fit ('%sUnpolished')"),
                 name, name), call. = FALSE)
  }
  if (is.null(unpolished) || !methods::is(unpolished, "SpaNormFit") ||
      isPolished(unpolished) || unpolished$ngenes != fit$ngenes ||
      unpolished$ncells != fit$ncells) {
    stop(sprintf(paste0("the '%s' fit is polished but no matching unpolished fit ",
                        "is stored as '%sUnpolished', so there is nothing to ",
                        "polish from; rerun '%s'"),
                 name, name, if (name == "SpaNormNull") "SpaNormSVG" else "SpaNorm"),
         call. = FALSE)
  }
  unpolished
}

# Shared by the SpatialExperiment and Seurat methods: choose the input fits,
# polish them, and normalise with the polished full fit. Returns the fits to
# store, by name, and the normalised matrix; the caller writes both into its
# container. `...` goes to polishNB() (engine, batch.size, backend, BPPARAM,
# block.size, psi.range).
.polishSpaNormCore <- function(stored, emat, adj.method, scale.factor,
                               psi.method, ls, cells, null, maxit, tol,
                               overwrite, verbose,
                               BPPARAM = BiocParallel::SerialParam(), ...) {
  msgfun <- if (verbose) message else function(...) {}
  # every refusal happens before any polishing
  full.in <- .polishSource(stored$SpaNorm, stored$SpaNormUnpolished, "SpaNorm",
                           overwrite)
  null.in <- if (null && !is.null(stored$SpaNormNull)) {
    .polishSource(stored$SpaNormNull, stored$SpaNormNullUnpolished,
                  "SpaNormNull", overwrite)
  }
  nsteps <- 2L + !is.null(null.in)
  step <- 0L
  progress <- function(msg) {
    step <<- step + 1L
    msgfun(sprintf("(%d/%d) %s", step, nsteps, msg))
  }

  progress("Polishing SpaNorm model")
  full <- .polishSpaNormFit(full.in, emat, psi.method = psi.method, ls = ls,
                            cells = cells, maxit = maxit, tol = tol,
                            verbose = verbose, BPPARAM = BPPARAM, ...)
  fits <- list(SpaNorm = full, SpaNormUnpolished = full.in)
  if (!is.null(null.in)) {
    progress("Polishing null SpaNorm model")
    fits$SpaNormNull <- .polishSpaNormFit(null.in, emat, psi.method = psi.method,
                                          ls = ls, cells = cells, maxit = maxit,
                                          tol = tol, verbose = verbose,
                                          BPPARAM = BPPARAM, ...)
    fits$SpaNormNullUnpolished <- null.in
  }

  progress("Normalising data")
  adj.fun <- getAdjustmentFun(full$gene.model, adj.method)
  normmat <- normaliseBlocked(adj.fun, emat, scale.factor, full,
                              BPPARAM = BPPARAM,
                              block.size = .normaliseBlockSize(emat))

  list(fits = fits, normmat = normmat)
}

#' Polish one SpaNorm-model fit
#'
#' Maps the fit onto polishNB()'s generic problem (.spaNormPolishProblem():
#' the shared library-size coefficient as a per-cell offset, gmean as an
#' unpenalised intercept column, the wtype penalty times ncells), converges
#' every gene with counts in the polished cells, and writes the result back
#' into the fit. Genes with no counts there are not passed to polishNB(): the
#' engine drives their intercept toward -Inf and reports them polished, so
#' they keep their input coefficients and dispersion and are recorded with
#' polished = FALSE. With ls = "joint", .polishSharedLS()
#' (R/polishSpaNormJoint.R) then moves the shared a1 to the joint optimum over
#' the polished genes, and every gene, held out or not, takes the new a1.
#'
#' @param fit a SpaNormFit, including one saved before the polish slot
#'   existed (assigning the slot adds it).
#' @param Y the counts, genes x cells, matching the fit (dense, sparse or
#'   DelayedArray).
#' @param psi.method,ls,cells,maxit,tol as in polishSpaNorm().
#' @param verbose logical; report progress.
#' @param ... passed to polishNB().
#' @return the polished SpaNormFit.
#' @noRd
.polishSpaNormFit <- function(fit, Y, psi.method = "fixed", ls = "fixed",
                              cells = "all", maxit = 50L, tol = 1e-8,
                              verbose = FALSE, ...) {
  if (nrow(Y) != fit$ngenes || ncol(Y) != fit$ncells) {
    stop(sprintf("the counts (%d x %d) do not match the SpaNorm fit (%d genes x %d cells/spots)",
                 nrow(Y), ncol(Y), fit$ngenes, fit$ncells), call. = FALSE)
  }
  prob <- .spaNormPolishProblem(fit, cells)
  Yc <- if (all(prob$cells_idx)) Y else Y[, prob$cells_idx, drop = FALSE]
  # counts are non-negative, so a zero row sum is an all-zero gene
  keep <- which(Matrix::rowSums(Yc) > 0)
  nzero <- nrow(Yc) - length(keep)
  if (verbose && nzero > 0) {
    message(sprintf("  %d gene%s with no counts in the polished cells/spots %s left unpolished",
                    nzero, if (nzero == 1) "" else "s", if (nzero == 1) "is" else "are"))
  }

  # updated in place, so an unpolished gene keeps its input values bit for bit
  # and the fit keeps its own dimnames
  gmean <- fit@gmean
  alpha <- fit@alpha
  psi <- fit@psi
  a1 <- prob$a1
  genes <- data.frame(iterations = rep(0L, nrow(Yc)), psi_fitnb = fit@psi,
                      restarted = FALSE, capped = FALSE, singular = FALSE,
                      psi_bound = FALSE, polished = FALSE, loglik = NA_real_)
  # the joint step's record (ls = "joint" only); its sums run over the polished
  # genes, and a1 stays shared, so held-out genes take the new a1 below
  j <- NULL
  if (length(keep)) {
    Yk <- Yc[keep, , drop = FALSE]
    probk <- prob
    probk$A0 <- prob$A0[keep, , drop = FALSE]
    probk$psi <- prob$psi[keep]
    pol <- polishNB(Yk, probk$X, probk$A0, probk$psi, lambda.a = probk$pen,
                    offset = probk$offset,
                    start.cols = c(TRUE, rep(FALSE, ncol(probk$X) - 1)),
                    psi.method = psi.method, maxit = maxit, tol = tol,
                    verbose = verbose, ...)
    if (ls == "joint") {
      j <- .polishSharedLS(Yk, probk, pol, psi.method = psi.method,
                           maxit = maxit, tol = tol, verbose = verbose, ...)
      pol <- j$pol
      a1 <- j$a1
    }
    gmean[keep] <- pol$alpha[, 1]
    alpha[keep, -1] <- pol$alpha[, -1, drop = FALSE]
    psi[keep] <- pol$psi
    genes[keep, names(pol$polish)] <- pol$polish
    genes$loglik[keep] <- pol$loglik
  } else if (ls == "joint") {
    j <- .lsUninformed(NULL, a1)
  }
  # the library-size coefficient is one value shared by every gene
  alpha[, 1] <- a1
  if (!is.null(rownames(Y))) rownames(genes) <- rownames(Y)

  fit@gmean <- gmean
  fit@alpha <- alpha
  fit@psi <- psi
  fit@polish <- list(
    settings = c(
      list(psi.method = psi.method, ls = ls, cells = cells,
           maxit = maxit, tol = tol, pen = prob$pen,
           a1.input = prob$a1, a1 = a1),
      if (ls == "joint") {
        list(ls.iterations = j$iterations, ls.maxit = j$maxit.ls, ls.tol = j$tol.ls,
             ls.score = j$score, ls.se = j$se, ls.singular = j$singular,
             ls.converged = j$converged)
      },
      list(SpaNorm = as.character(utils::packageVersion("SpaNorm")))),
    genes = genes)
  methods::validObject(fit)
  fit
}
