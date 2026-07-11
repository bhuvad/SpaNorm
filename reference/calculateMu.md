# Compute fitted means from a negative binomial GLM fit

Computes the fitted mean matrix `mu = exp(gmean + tcrossprod(alpha, W))`
from the per-gene coefficients of a negative binomial GLM, with optional
per-gene winsorisation of the log-means to `median +/- winsor * MAD` to
bound the influence of extreme fitted values. Exposed so downstream
packages (e.g. spiDE) can reconstruct fitted means from a
[`fitNB`](https://bhuvad.github.io/spaNorm/reference/fitNB.md) result;
for the generic (intercept-free) fit pass `gmean = rep(0, nrow(alpha))`.

## Usage

``` r
calculateMu(gmean, alpha, W, winsor = DEFAULT_WINSOR)
```

## Arguments

- gmean:

  a numeric vector of per-gene intercepts (length `nrow(alpha)`).

- alpha:

  a genes x covariates matrix of coefficients.

- W:

  a cells x covariates design matrix.

- winsor:

  a numeric, the number of MADs at which per-gene log-means are
  winsorised (default 4); `Inf` disables winsorisation.

## Value

a genes x cells matrix of fitted means.

## Examples

``` r
set.seed(1)
Y <- matrix(rpois(20 * 50, 5), 20, 50)
W <- cbind(1, scale(seq_len(50)))
fit <- fitNB(Y, W, verbose = FALSE)
mu <- calculateMu(fit$gmean, fit$alpha, W)
dim(mu)
#> [1] 20 50
```
