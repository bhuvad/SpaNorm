#' @importFrom stats dnbinom qnbinom pnbinom mad median quantile model.matrix pf p.adjust
NULL

#' Model-based spatially variable gene (SVG) calling
#'
#' Spatially variable gene (SVG) calling using the SpaNorm model.
#'
#' @param spe a SpatialExperiment or Seurat object, with the count data stored in 'counts' or 'data' assays respectively, and a SpaNorm model fit.
#' @param backend a character, specifying the backend to use for computations. Options are "auto" (default), "cpu", or "gpu". If "auto", it will use GPU if available, otherwise CPU.
#' @param verbose a logical, specifying whether to show update messages (default TRUE).
#'
#' @details SpaNorm SVG calling works by using the SpaNorm model fit for data normalisation to perform a likelihood ratio test (LRT). The model used for normalisation is considered to be the full model. A second nested model is fit without the splines representing biology. These nested models are then compared using a LRT to identify genes where the splines representing biology contain strong signal.
#'
#' @return a SpatialExperiment or Seurat object with F-statistics, false discovery rates (FDRs). For SpatialExperiment objects, these are stored in the rowData.
#' @name SpaNormSVG
#'
#' @examples
#' 
#' library(SpatialExperiment)
#' library(ggplot2)
#' 
#' data(HumanDLPFC)
#' 
#' HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#' HumanDLPFC = SpaNormSVG(HumanDLPFC)
#' head(rowData(HumanDLPFC))
#' 
#' @export
#'
setGeneric("SpaNormSVG", function(
    spe,
    backend = c("auto", "cpu", "gpu"),
    verbose = TRUE) {
  standardGeneric("SpaNormSVG")
})

#' @rdname SpaNormSVG
setMethod(
  "SpaNormSVG",
  signature("SpatialExperiment"),
  function(spe, backend, verbose) {
    checkSPE(spe)

    # message function depending on verbose param
    msgfun = ifelse(verbose, message, \(...){})
    # Add progress tracking
    total_steps = 3
    current_step = 0
    
    report_progress <- function(msg) {
      if (verbose) {
        current_step <<- current_step + 1
        msgfun(sprintf("(%d/%d) %s", current_step, total_steps, msg))
      }
    }

    # Extract counts
    emat = SummarizedExperiment::assay(spe, "counts")
    
    # Check previous results
    results = getSVGResults(spe, stop_if_missing = FALSE)
    svg.cols = colnames(results)
    if (!is.null(results)) {
      warning("SVG results exist in 'spe' and will be overwritten")
      if (is(spe, "SpatialExperiment")) {
        cols = setdiff(colnames(SummarizedExperiment::rowData(spe)), svg.cols)
        SummarizedExperiment::rowData(spe) = SummarizedExperiment::rowData(spe)[, cols]
      }
    }

    # Retrieve and validate SpaNorm model
    report_progress("Retrieving SpaNorm model")
    fit.spanorm = getSpaNormFit(spe)

    # Fit, retrieve, or re-polish the nested null model, keeping it paired
    # with the full fit's polish state -- svgTest() compares the two fits by
    # a likelihood ratio, which is only a comparison of nested models when
    # both sides are estimated the same way (see .svgPairedNull()).
    fit.technical = getSpaNormFit(spe, null = TRUE, validate = FALSE)
    unpolished.null = S4Vectors::metadata(spe)$SpaNormNullUnpolished
    paired = .svgPairedNull(fit.spanorm, fit.technical, unpolished.null, emat,
                            msgfun, report_progress, backend)
    fit.technical = paired$fit.technical
    S4Vectors::metadata(spe)$SpaNormNull = fit.technical
    if (!is.null(paired$unpolished)) {
      S4Vectors::metadata(spe)$SpaNormNullUnpolished = paired$unpolished
    }

    # F-test
    report_progress("Finding SVGs")
    df.svg = svgTest(emat, fit.spanorm, fit.technical)
    SummarizedExperiment::rowData(spe) = cbind(SummarizedExperiment::rowData(spe), df.svg[rownames(spe), ])
    msgfun(sprintf("%d SVGs found (FDR < 0.05)", sum(df.svg$svg.fdr < 0.05)))

    spe
  }
)

#' @rdname SpaNormSVG
setMethod(
  "SpaNormSVG",
  signature("Seurat"),
  function(spe, backend, verbose) {
    stop(
      "'SpaNormSVG' currently supports 'SpatialExperiment' objects only. ",
      "Please convert your Seurat object to a SpatialExperiment to call SVGs."
    )
  }
)

#' Are a full and null SpaNormFit "polished alike"?
#'
#' TRUE iff both are unpolished, or both are polished and their
#' \code{\link{.polishSlot}()$settings} agree on \code{psi.method}, \code{ls}
#' and \code{cells}, and under \code{psi.method = "profile"} on
#' \code{psi.range} (a missing one, from a fit polished before it was
#' recorded, is polishNB()'s default; final review M-1) -- the settings that
#' define the fitted objective (\code{maxit}/\code{tol} only control how
#' tightly it is converged, and \code{a1} legitimately differs between the
#' full and null fits, so none of those is compared here). This is the ONLY
#' function either
#' \code{SpaNormSVG()} or \code{svgTest()} uses to decide whether a pair is
#' comparable (Task 8 fix round 1, ruling 1): two \code{polishSpaNorm()}
#' calls -- one polishing only the full fit, with different settings -- give
#' a pair that is \code{isPolished()} on both sides but not polished alike,
#' and comparing them silently changes the LRT (measured: median
#' \code{svg.F} +38%, FDR<0.05 calls 30 -> 33, on \code{.polish_spe()}).
#' @noRd
.polishedAlike <- function(fit.spanorm, fit.technical) {
  full.polished <- isPolished(fit.spanorm)
  null.polished <- isPolished(fit.technical)
  if (full.polished != null.polished) return(FALSE)
  if (!full.polished) return(TRUE)
  sf <- .polishSlot(fit.spanorm)$settings
  sn <- .polishSlot(fit.technical)$settings
  identical(sf$psi.method, sn$psi.method) && identical(sf$ls, sn$ls) &&
    identical(sf$cells, sn$cells) &&
    (!identical(sf$psi.method, "profile") ||
       identical(.polishPsiRange(sf), .polishPsiRange(sn)))
}

#' What `.polishedAlike()` found different, for an error/message
#'
#' Only meaningful when \code{.polishedAlike()} is FALSE for the same pair:
#' names the polish state when that is what differs, else the first
#' mismatched setting among \code{psi.method}/\code{ls}/\code{cells}, then
#' \code{psi.range} under \code{"profile"} (the same fields
#' \code{.polishedAlike()} compares) with both its values.
#' @noRd
.polishAlikeDiff <- function(fit.spanorm, fit.technical) {
  full.polished <- isPolished(fit.spanorm)
  null.polished <- isPolished(fit.technical)
  if (full.polished != null.polished) {
    return(sprintf("the polish state differs (full: %spolished, null: %spolished)",
                   if (full.polished) "" else "not ",
                   if (null.polished) "" else "not "))
  }
  sf <- .polishSlot(fit.spanorm)$settings
  sn <- .polishSlot(fit.technical)$settings
  for (field in c("psi.method", "ls", "cells")) {
    if (!identical(sf[[field]], sn[[field]])) {
      return(sprintf("'%s' differs (full: %s, null: %s)", field,
                     deparse(sf[[field]]), deparse(sn[[field]])))
    }
  }
  if (identical(sf$psi.method, "profile") &&
      !identical(.polishPsiRange(sf), .polishPsiRange(sn))) {
    return(sprintf("'psi.range' differs (full: %s, null: %s)",
                   deparse(.polishPsiRange(sf)), deparse(.polishPsiRange(sn))))
  }
  "no difference found"  # unreachable when .polishedAlike() is FALSE
}

#' Fit, retrieve, or re-polish the null model so it is paired with the full
#' fit's polish state (Task 8, spec section 5, "the pairing rule"):
#' \code{svgTest()} compares the full and null fits by a likelihood ratio,
#' which is only a comparison of nested models when both sides are estimated
#' the same way -- polished alike, or unpolished alike (\code{.polishedAlike()}).
#'
#' Four cases:
#' \itemize{
#'   \item no stored null: fit it (\code{fitSpaNormTechnical()}), and polish
#'     it too if the full fit is polished, with ALL of the full fit's
#'     recorded polish settings (\code{psi.method}, \code{ls}, \code{cells},
#'     \code{maxit}, \code{tol}, and \code{psi.range} under \code{"profile"});
#'   \item a stored null that is polished alike (\code{.polishedAlike()}
#'     TRUE): use it as is;
#'   \item a stored null whose POLISH STATE does not match: if the full fit
#'     is polished and the null is not, polish the null here (the common
#'     case: a null fit before \code{polishSpaNorm()} existed, or
#'     \code{SpaNormSVG()} rerun after \code{polishSpaNorm(spe, null =
#'     FALSE)}); if the full fit is unpolished and the null IS polished,
#'     stop rather than silently un-polishing the null or polishing the full
#'     fit here -- that is \code{polishSpaNorm()}'s job, not
#'     \code{SpaNormSVG()}'s;
#'   \item both polished, but their SETTINGS differ (fix round 1: reachable
#'     through the public API by polishing the full and null fits in
#'     separate \code{polishSpaNorm()} calls with different arguments, e.g.
#'     \code{polishSpaNorm(spe)} then \code{polishSpaNorm(spe, overwrite =
#'     TRUE, null = FALSE, psi.method = "profile")}): re-polish the null
#'     with the full fit's settings, starting from the stored
#'     \code{SpaNormNullUnpolished} when there is one (that is what the
#'     public-API sequence above leaves behind), else refit it with
#'     \code{fitSpaNormTechnical()}, with a \code{message()} naming what
#'     differed.
#' }
#'
#' @param unpolished.stored the stored \code{SpaNormNullUnpolished}, if any
#'   (\code{NULL} otherwise); used as the polish's starting point in the
#'   fourth case above so the null's own optimum is not lost.
#' @return a list with \code{fit.technical} (the paired null) and
#'   \code{unpolished} -- the pre-polish null when one was polished here
#'   (\code{NULL} otherwise), so the caller can store it as
#'   \code{'SpaNormNullUnpolished'}, as \code{polishSpaNorm()} does.
#' @noRd
.svgPairedNull <- function(fit.spanorm, fit.technical, unpolished.stored, emat,
                           msgfun, report_progress, backend) {
  full.polished <- isPolished(fit.spanorm)
  polish.null <- function(nul) {
    # ALL of the full fit's recorded settings that .polishSpaNormFit()
    # accepts, not only the ones .polishedAlike() compares (fix round 1,
    # ruling 4): psi.method/ls/cells (and psi.range, under "profile") define
    # the objective, maxit/tol define how tightly it is converged, and a
    # re-polished null should match on both. psi.range goes to polishNB()
    # through `...`; a full fit polished before it was recorded used the
    # default, which .polishPsiRange() returns (final review M-1).
    settings <- .polishSlot(fit.spanorm)$settings
    .polishSpaNormFit(nul, emat, psi.method = settings$psi.method,
                      ls = settings$ls, cells = settings$cells,
                      maxit = settings$maxit, tol = settings$tol,
                      verbose = FALSE, name = "SpaNormNull",
                      psi.range = .polishPsiRange(settings))
  }

  if (is.null(fit.technical)) {
    report_progress(if (full.polished) "Fitting and polishing Null SpaNorm model"
                    else "Fitting Null SpaNorm model")
    fit.technical <- fitSpaNormTechnical(emat, fit.spanorm, msgfun, backend = backend)
    if (!full.polished) {
      return(list(fit.technical = fit.technical, unpolished = NULL))
    }
    return(list(fit.technical = polish.null(fit.technical), unpolished = fit.technical))
  }

  if (.polishedAlike(fit.spanorm, fit.technical)) {
    report_progress("Retrieving Null SpaNorm model")
    return(list(fit.technical = fit.technical, unpolished = NULL))
  }

  null.polished <- isPolished(fit.technical)
  if (full.polished && !null.polished) {
    message("the stored null SpaNorm model is not polished; polishing it ",
            "to match the polished full model")
    report_progress("Polishing Null SpaNorm model")
    return(list(fit.technical = polish.null(fit.technical), unpolished = fit.technical))
  }
  if (!full.polished) {
    # full unpolished, null polished: .polishedAlike() would also be FALSE
    # here on a settings difference, but that cannot arise (an unpolished
    # fit has no settings), so this is always the polish-state mismatch
    stop("the full SpaNorm model is not polished but the stored null SpaNorm ",
         "model is; call 'polishSpaNorm()' to polish both (or re-fit an ",
         "unpolished null) before rerunning 'SpaNormSVG()'", call. = FALSE)
  }

  # both polished, but their settings differ (fix round 1, ruling 3): start
  # from the stored unpolished null when there is one (the public-API
  # sequence that creates this mismatch leaves it behind), else refit
  message("the stored null SpaNorm model is polished with different settings ",
          "than the full model (", .polishAlikeDiff(fit.spanorm, fit.technical),
          "); re-polishing it to match")
  report_progress("Polishing Null SpaNorm model")
  source <- if (!is.null(unpolished.stored) &&
                methods::is(unpolished.stored, "SpaNormFit") &&
                !isPolished(unpolished.stored)) {
    unpolished.stored
  } else {
    fitSpaNormTechnical(emat, fit.spanorm, msgfun, backend = backend)
  }
  list(fit.technical = polish.null(source), unpolished = source)
}

fitSpaNormTechnical <- function(Y, fit.spanorm, msgfun, ...) {
  # the null must be fitted to the same cells as the full model: the penalty
  # below is scaled by ncol(Y), as in fitSpaNorm()
  if (ncol(Y) != fit.spanorm$ncells) {
    stop(sprintf("number of cells in data (%d) differs from the SpaNorm fit (%d)", ncol(Y), fit.spanorm$ncells))
  }
  msgfun(sprintf("%d cells/spots sampled to fit model", sum(fit.spanorm$sampling != "all")))

  # select technical covariates
  wtype = fit.spanorm$wtype
  W = fit.spanorm$W[, wtype != "biology", drop = FALSE]
  wtype = wtype[wtype != "biology"]

  # create lambda.a vector
  lambda.a.vec = rep(0, length(wtype))
  lambda.a.vec[wtype == "biology"] = fit.spanorm$lambda.a[1]
  lambda.a.vec[wtype == "ls"] = fit.spanorm$lambda.a[2]
  lambda.a.vec = lambda.a.vec[-1]

  # fit model
  fit.technical = fitSpaNormNB(
    Y,
    W,
    fit.spanorm$sampling != "all",
    maxn.psi = sum(fit.spanorm$sampling == "dispersion"),
    lambda.a = lambda.a.vec * ncol(Y),
    msgfun = msgfun,
    ...,
    is.spanorm = TRUE
  )

  # create object
  fit.technical = SpaNormFit(
    ngenes = nrow(Y),
    ncells = ncol(Y),
    gene.model = fit.spanorm$gene.model,
    df.tps = fit.spanorm$df.tps,
    sample.p = fit.spanorm$sample.p,
    lambda.a =fit.spanorm$lambda.a,
    batch = fit.spanorm$batch,
    W = W,
    alpha = fit.technical$alpha,
    gmean = fit.technical$gmean,
    psi = fit.technical$psi,
    wtype = wtype,
    loglik = fit.technical$loglik,
    sampling = fit.technical$sampling
  )

  return(fit.technical)
}

svgTest <- function(Y, fit.spanorm, fit.technical) {
  # Add validation for input types
  if (!is.matrix(Y) && !methods::is(Y, "Matrix")) {
    stop("Y must be a matrix or Matrix object")
  }
  if (!methods::is(fit.spanorm, "SpaNormFit")) {
    stop("fit.spanorm must be a SpaNormFit object") 
  }
  if (!methods::is(fit.technical, "SpaNormFit")) {
    stop("fit.technical must be a SpaNormFit object")
  }
  # svgTest() is internal, but reachable via SpaNorm:::svgTest(), so this
  # guard lives here as well as in SpaNormSVG()'s pairing logic
  # (.svgPairedNull()): the LRT is only a comparison of nested models when
  # both fits were estimated the same way. .polishedAlike() is the sole
  # decision (fix round 1, ruling 1): isPolished() equality alone is not
  # enough -- two polishSpaNorm() calls can leave both fits "polished" with
  # different psi.method/ls/cells, which changes the objective just as much.
  if (!.polishedAlike(fit.spanorm, fit.technical)) {
    stop(sprintf("the full and null SpaNorm fits must be polished alike: %s",
                 .polishAlikeDiff(fit.spanorm, fit.technical)), call. = FALSE)
  }

  # Existing dimension checks
  if (length(unique(c(nrow(Y), fit.spanorm$ngenes, fit.technical$ngenes))) != 1) {
    stop("number of genes differ between SpaNorm fits and/or data")
  }
  if (length(unique(c(ncol(Y), fit.spanorm$ncells, fit.technical$ncells))) != 1) {
    stop("number of cells differ between SpaNorm fits and/or data")
  }
  if (!all(colnames(fit.technical$W) %in% colnames(fit.spanorm$W)) | ncol(fit.technical$W) == ncol(fit.spanorm$W)) {
    stop("technical model is not nested in the full model")
  }

  # convert to regular matrix
  Y = as.matrix(Y)

  # nested model residual deviance
  # calculate mu
  mu = calculateMu(fit.spanorm$gmean, fit.spanorm$alpha, fit.spanorm$W)
  # winsorize dispersion parameters (large disp slow the process)
  psi = winsorisePsi(fit.spanorm$psi)
  loglik.spanorm = rowSums(dnbinom(Y, mu = mu, size = 1 / psi, log = TRUE))

  # nested model residual deviance
  # calculate mu
  mu = calculateMu(fit.technical$gmean, fit.technical$alpha, fit.technical$W)
  # winsorize dispersion parameters (large disp slow the process)
  psi = winsorisePsi(fit.technical$psi)
  loglik.technical = rowSums(dnbinom(Y, mu = mu, size = 1 / psi, log = TRUE))

  # F-test
  df1 = ncol(fit.spanorm$W) - ncol(fit.technical$W)
  df2 = ncol(Y) - ncol(fit.spanorm$W)
  F.raw = 2 * (loglik.spanorm - loglik.technical) / df1
  # Threshold to 0 due to convergence issues
  F.lrt = pmax(F.raw, 0)
  p.val = pf(F.lrt, df1, df2, lower.tail = FALSE)
  fdr = p.adjust(p.val, method = "fdr")

  # create results table
  df.svg = data.frame(
    svg.F = F.lrt,
    svg.p = p.val,
    svg.fdr = fdr
  )
  # the unclamped statistic, for measuring how far below 0 the clamp bites
  # (additive only: svg.F above is unchanged); see Task 8 ruling 2
  attr(df.svg, "F.raw") = F.raw

  return(df.svg)
}

#' Export top SVG results to a data frame
#' 
#' @param spe a SpatialExperiment object with SVG results from SpaNormSVG.
#' @param n a numeric, specifying the number of top SVGs to call.
#' @param fdr a numeric, specifying the false discovery rate (FDR) threshold for calling SVGs.
#' 
#' @return A data frame containing the top SVGs from F-test results including F-statistics, p-values and FDR.
#' @examples
#' 
#' library(SpatialExperiment)
#' library(ggplot2)
#'
#' data(HumanDLPFC)
#'
#' HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#' HumanDLPFC = SpaNormSVG(HumanDLPFC)
#' topSVGs = topSVGs(HumanDLPFC, n = 10)
#' @export
topSVGs <- function(spe, n = 10, fdr = 1) {
  stopifnot(n > 0)
  stopifnot(fdr >= 0 && fdr <= 1)
  checkSPE(spe)
  
  # Get SVG results
  results = getSVGResults(spe)
  
  # Filter and sort results
  results = results[results$svg.fdr <= fdr, , drop = FALSE]
  results = results[order(results$svg.fdr), , drop = FALSE]
  n = min(n, nrow(results))
  results = results[seq_len(n), , drop = FALSE]
  
  return(results)
}

#' Get SVG results from a SpatialExperiment object
#'
#' @param spe a SpatialExperiment object
#' @param stop_if_missing logical indicating whether to stop if SVG results are missing (default TRUE)
#' @return A data frame containing SVG results if they exist
#' @keywords internal
getSVGResults <- function(spe, stop_if_missing = TRUE) {
  # Check input type
  if (is(spe, "SpatialExperiment")) {
    rowData = SummarizedExperiment::rowData(spe)
  } else {
    stop("'spe' must be a SpatialExperiment object")
  }
  
  # Define expected columns
  svg.cols = c("svg.F", "svg.p", "svg.fdr")
  
  # Check if results exist
  has_results = all(svg.cols %in% colnames(rowData))
  
  if (!has_results) {
    if (stop_if_missing) {
      stop("SVG results not found in 'spe'. Please run 'SpaNormSVG' first.")
    } else {
      return(NULL)
    }
  }
  
  # Get results
  results = as.data.frame(rowData[, svg.cols, drop = FALSE])
  
  return(results)
}
