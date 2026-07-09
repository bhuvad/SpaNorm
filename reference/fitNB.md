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

- backend:

  a character, the compute backend ('auto', 'cpu', or 'gpu').

- verbose:

  a logical, whether to print progress messages (default TRUE).

## Value

a list with per-gene coefficients `alpha` (genes x covariates),
dispersions `psi`, a `gmean` element (always zero – the generic fit has
no intercept term), the `sampling` factor, and per-iteration `loglik`.

## Examples

``` r
set.seed(1)
Y <- matrix(rpois(20 * 50, 5), 20, 50)
W <- cbind(1, scale(seq_len(50))) # intercept + one covariate
fit <- fitNB(Y, W, verbose = FALSE)
str(fit$alpha)
#>  num [1:20, 1:2] 1.65 1.6 1.72 1.62 1.62 ...
```
