#' @name SpaNormFit
#'
#' @title An S4 class to store a SpaNorm model fit
#' @slot ngenes a numeric, specifying the number of genes in the dataset.
#' @slot ncells a numeric, specifying the number of cells/spots in the dataset.
#' @slot gene.model a character, specifying the gene-specific model to used (see `getGeneModels()`).
#' @slot df.tps an integer, specifying the degrees of freedom to used for the thin plate spline.
#' @slot sample.p a numeric, specifying the proportion of samples used to approximated the model.
#' @slot lambda.a a numeric, specifying the shinkage parameter used.
#' @slot batch a vector or matrix, specifying the batch design used (if any). 
#' @slot W a matrix, specifying the covariate matrix of the linear model.
#' @slot alpha a matrix, specifying the coefficients of the linear model.
#' @slot gmean a numeric, specifying the mean estimate for each gene in the linear model.
#' @slot psi a numeric, specifying the over-dispersion parameter for each gene if a negative binomial model was used (or a vector of NAs if another gene model is used).
#' @slot wtype a factor, specifying the covariate types of columns in the covariate matrix, W. These could be "biology", "ls", or "batch".
#' @slot loglik a numeric, specifying the log-likelihood of the model at each external iteration.
#' @slot sampling a factor, specifying the cells/spots used for dispersion estimation ('dispersion'), GLM fitting ('glm' and 'dispersion'), all other cells/spots ('all').
#' 
#' @param x an object of class SpaNormFit.
#' @param name a character, specifying the name of the slot to retrieve.
#' @return Return value varies depending on method.
#' @examples
#' example(SpaNorm)
NULL

#' @rdname SpaNormFit
#' @export
#' @import methods
setClass(
  Class = "SpaNormFit",
  slots = c(
    ngenes = "numeric",
    ncells = "numeric",
    gene.model = "character",
    df.tps = "integer",
    sample.p = "numeric",
    lambda.a = "numeric",
    batch = "ANY",
    W = "matrix",
    alpha = "matrix",
    gmean = "numeric",
    psi = "numeric",
    wtype = "factor",
    loglik = "numeric",
    sampling = "factor"
  )
)

#' @rdname SpaNormFit
#' @export
setMethod(
  f = "$",
  signature = "SpaNormFit",
  definition = function(x, name) {
    return(slot(x, name))
  }
)

setMethod(
  f = "show",
  signature = "SpaNormFit",
  definition = function(object) {
    # present sampling
    sampling = object@sampling
    sampling = sprintf("all (%d), glm (%d), dispersion (%d)", length(sampling), sum(sampling != "all"), sum(sampling == "dispersion"))

    cat(
      is(object)[[1]],
      sprintf("Data: %d genes, %d cells/spots", object@ngenes, object@ncells),
      sprintf("Gene model: %s", object@gene.model),
      sprintf("Degrees of freedom for TPS (y,x): Biology (%d,%d), LS (%d,%d)", object@df.tps[2], object@df.tps[1], object@df.tps[2], object@df.tps[1]),
      sprintf("Spots/cells sampled: %s%%", signif(object@sample.p * 100, 3)),
      sprintf("Regularisation parameter: %s", signif(object@lambda.a, 3)),
      sprintf("Batch: %s", utils::capture.output(utils::str(object@batch))),
      sprintf("log-likelihood (per-iteration): %s", utils::capture.output(utils::str(object@loglik))),
      sprintf("W: %s", utils::capture.output(utils::str(object@W))),
      sprintf("alpha: %s", utils::capture.output(utils::str(object@alpha))),
      sprintf("gmean: %s", utils::capture.output(utils::str(object@gmean))),
      sprintf("psi: %s", utils::capture.output(utils::str(object@psi))),
      sprintf("wtype: %s", utils::capture.output(utils::str(object@wtype))),
      sprintf("sampling: %s", sampling),
      sep = "\n"
    )
  }
)

validSpaNormFit <- function(object) {
  if (any(object@df.tps <= 0)) {
    stop("'df.tps' should be greater than 0")
  }
  if (!object@gene.model %in% getGeneModels()) {
    stop(sprintf("'gene.model' should be one of: %s", paste(getGeneModels(), collapse = ", ")))
  }
  if (object@sample.p <= 0 | object@sample.p > 1) {
    stop("'sample.p' should be in the interval (0,1]")
  }
  if (object@lambda.a <= 0) {
    stop("'lambda.a' should be greater than 0")
  }
  if (any(object@loglik > 0)) {
    stop("'loglik' should be less than or equal to 0")
  }
  if (!all(object@sampling %in% c("all", "glm", "dispersion")) || !is.factor(object@sampling)) {
    stop("'sampling' should be a factor containing 'all', 'glm', or 'dispersion'")
  }

  # check dimensions
  if (length(unique(c(ncol(object@alpha), ncol(object@W), length(object@wtype)))) > 1) {
    stop("ncol of 'alpha', ncol of 'W', and length of 'wtype' do not match")
  }
  if (nrow(object@alpha) != object@ngenes) {
    stop("nrow of 'alpha' does not match 'ngenes")
  }
  if (length(object@gmean) != object@ngenes) {
    stop("length of 'gmean' does not match 'ngenes")
  }
  if (length(object@psi) != object@ngenes) {
    stop("length of 'psi' does not match 'ngenes")
  }
  if (nrow(object@W) != object@ncells) {
    stop("nrow of 'W' does not match 'ncells")
  }
  if (!is.null(object@batch) && is.vector(object@batch) && length(object@batch) != object@ncells) {
    stop("length of 'batch' does not match 'ncells'")
  }
  if (length(object@sampling) != object@ncells) {
    stop("length of 'sampling' does not match 'ncells'")
  }
  if (!is.null(object@batch) && is.matrix(object@batch) && nrow(object@batch) != object@ncells) {
    stop("nrow of 'batch' does not match 'ncells'")
  }

  # check NAs
  if (any(is.na(object@W))) {
    stop("'W' cannot have missing values")
  }
  if (any(is.na(object@alpha))) {
    stop("'alpha' cannot have missing values")
  }
  if (any(is.na(object@gmean))) {
    stop("'gmean' cannot have missing values")
  }
  if (object@gene.model %in% c("nb") & any(is.na(object@psi))) {
    stop("'psi' cannot have missing values")
  }
  if (any(is.na(object@wtype))) {
    stop("'wtype' cannot have missing values")
  }
  if (any(is.na(object@sampling))) {
    stop("'sampling' cannot have missing values")
  }
  if (!is.null(object@batch) && any(is.na(object@batch))) {
    stop("'batch' cannot have missing values")
  }

  # check wtype levels
  if (!all(levels(object@wtype) %in% c("biology", "ls", "batch"))) {
    stop("'wtype' should have values that are either 'biology', 'ls', or 'batch'.")
  }

  TRUE
}

setValidity("SpaNormFit", validSpaNormFit)

SpaNormFit <- function(ngenes, ncells, gene.model, ..., df.tps, sample.p, lambda.a, W, alpha, gmean, wtype, loglik, batch = NULL, psi = NULL, sampling) {
  if (!gene.model %in% getGeneModels()) {
    stop(sprintf("'gene.model' should be one of: %s", paste(getGeneModels(), collapse = ", ")))
  }

  if (is.null(psi)) {
    psi = rep(0, ngenes)
  }
  if (!gene.model %in% c("nb")) {
    psi = rep(NA, ngenes)
  }
  if (is.character(wtype)) {
    wtype = factor(wtype, levels = c("batch", "biology", "ls"))
  }

  new(
    "SpaNormFit",
    ngenes = ngenes,
    ncells = ncells,
    gene.model = gene.model,
    df.tps = df.tps,
    sample.p = sample.p,
    lambda.a = lambda.a,
    batch = batch,
    W = W,
    alpha = alpha,
    gmean = gmean,
    psi = psi,
    wtype = wtype,
    loglik = loglik,
    sampling = sampling
  )
}
