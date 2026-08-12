# Fit a per-gene negative binomial GLM

Fits a per-gene negative binomial regression over an arbitrary design
matrix using SpaNorm's IRLS engine (iteratively reweighted least squares
with per-gene dispersion estimated by
[`edgeR::estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)).
This exposes the fitting machinery independently of SpaNorm's spatial
model: the fitted model is `log(mu) = W %*% t(alpha)` with no built-in
intercept, so encode one as a column of `W` when a per-gene baseline is
required. It corresponds to the internal fit with `is.spanorm = FALSE`
(no shared library-size column, ridge applied to every column of `W`).

## Usage

``` r
fitNB(
  Y,
  W,
  idx = rep(TRUE, ncol(Y)),
  lambda.a = 0,
  winsor = DEFAULT_WINSOR,
  maxit.psi = 25,
  maxit.nb = 50,
  tol = 1e-04,
  ...,
  offset = NULL,
  psi = NULL,
  backend = c("auto", "cpu", "gpu"),
  verbose = TRUE
)
```

## Arguments

- Y:

  a genes x cells matrix of counts (dense, sparse, or DelayedArray).

- W:

  a cells x covariates numeric design matrix.

- idx:

  a logical vector (length `ncol(Y)`) selecting the cells used to fit
  the model (default: all cells).

- lambda.a:

  a numeric ridge penalty on the columns of `W`: a single value or a
  per-column vector (default 0, i.e. unregularised).

- winsor:

  a numeric, the number of MADs at which per-gene coefficients are
  winsorised during fitting (default 4). Must be a single positive
  number; `Inf` disables winsorisation entirely.

- maxit.psi:

  a numeric, the maximum number of dispersion iterations.

- maxit.nb:

  a numeric, the maximum number of NB mean IRLS iterations.

- tol:

  a numeric, the convergence tolerance.

- ...:

  additional fitting parameters forwarded to the internal fitter, e.g.
  `maxn.psi` (dispersion-estimation subsample size) or `step.factor`
  (IRLS step-halving factor).

- offset:

  either `NULL` (the default, no offset) or a genes x cells numeric
  matrix with the same dimensions as `Y`, added to the linear predictor
  with a coefficient fixed at 1, so that each fitted log-mean
  `log(mu[g, c])` gains `offset[g, c]`. Use this (rather than an extra
  column of `W`) for a known, non-estimated effect – a fitted column's
  coefficient floats and can absorb signal that correlates with the
  known effect, whereas an offset cannot. The offset is subset by `idx`
  internally, alongside `Y`, and enters the dispersion estimation as
  well as the mean.

- psi:

  either `NULL` (the default: per-gene dispersions are estimated by the
  usual outer loop, via
  [`edgeR::estimateDisp`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html))
  or a numeric vector of length `nrow(Y)` of per-gene NB dispersions
  (`size = 1/psi`). Supplied dispersions are used as-is – no
  re-estimation, no winsorisation – and the outer dispersion loop is
  bypassed entirely: the coefficients come from a single IRLS fit at the
  given `psi` (so `maxit.psi` is ignored, and the returned `loglik` has
  one element). The dispersions need not come from an identical design –
  values estimated on a design nested in (or equal to) the fitting
  design are appropriate, e.g. pooled across a coarser model, which errs
  conservative. `offset` and `psi` compose.

- backend:

  a character, the compute backend ('auto', 'cpu', or 'gpu').

- verbose:

  a logical, whether to print progress messages (default TRUE).

## Value

a list with per-gene coefficients `alpha` (genes x covariates),
dispersions `psi`, a `gmean` element (always zero – the generic fit has
no intercept term), the `sampling` factor (with a supplied `psi` no
dispersion subsample is drawn, so its `"dispersion"` level is absent:
cells are `"glm"` if used for fitting, else `"all"`), and
per-outer-iteration `loglik` (length 1 when `psi` is supplied).

## Examples

``` r
set.seed(1)
Y <- matrix(rpois(20 * 50, 5), 20, 50)
W <- cbind(1, scale(seq_len(50))) # intercept + one covariate
fit <- fitNB(Y, W, verbose = FALSE)
str(fit$alpha)
#>  num [1:20, 1:2] 1.65 1.6 1.72 1.62 1.62 ...
```
