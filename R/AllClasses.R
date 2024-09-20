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
    isbio = "logical",
    loglik = "numeric"
  )
)

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
    cat(
      is(object)[[1]],
      sprintf("Data: %d genes, %d cells/spots", object@ngenes, object@ncells),
      sprintf("Gene model: %s", object@gene.model),
      sprintf("Degrees of freedom (TPS): %d", object@df.tps),
      sprintf("Spots/cells sampled: %s%%", signif(object@sample.p * 100, 3)),
      sprintf("Regularisation parameter: %s", signif(object@lambda.a, 3)),
      sprintf("Batch: %s", utils::capture.output(str(object@batch))),
      sprintf("log-likelihood (per-iteration): %s", utils::capture.output(str(object@loglik))),
      sprintf("W: %s", utils::capture.output(str(object@W))),
      sprintf("alpha: %s", utils::capture.output(str(object@alpha))),
      sprintf("gmean: %s", utils::capture.output(str(object@gmean))),
      sprintf("psi: %s", utils::capture.output(str(object@psi))),
      sep = "\n"
    )
  }
)

validSpaNormFit <- function(object) {
  if (object@df.tps <= 0) {
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
  if (!any(object@isbio)) {
    stop("'isbio' should have at least one TRUE value")
  }

  # check dimensions
  if (length(unique(c(ncol(object@alpha), ncol(object@W), length(object@isbio)))) > 1) {
    stop("ncol of 'alpha', ncol of 'W', and length of 'isbio' do not match")
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

  TRUE
}

setValidity("SpaNormFit", validSpaNormFit)

SpaNormFit <- function(ngenes, ncells, gene.model, ..., df.tps, sample.p, lambda.a, W, alpha, gmean, isbio, loglik, batch = NULL, psi = NULL) {
  if (is.null(psi)) {
    psi = rep(0, ngenes)
  }
  if (!gene.model %in% c("nb")) {
    psi = rep(NA, ngenes)
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
    isbio = isbio,
    loglik = loglik
  )
}
