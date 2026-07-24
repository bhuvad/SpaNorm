# Negative binomial density (mean/size parameterisation), GPU-aware

As [`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html) (with `mu`
instead of `prob`), but runs the log-pmf on the accelerator when `x` or
`mu` is a torch tensor and an accelerator is active. Targets the domain
SpaNorm actually uses (non-negative integer `x`, `mu > 0`) – outside
that domain this is not a drop-in replacement for
[`stats::dnbinom`](https://rdrr.io/r/stats/NegBinomial.html) (e.g.
`mu = 0` gives `NaN` rather than a point mass at 0). Exposed so
downstream packages (e.g. spiDE) computing their own per-gene NB
log-likelihoods reuse the same GPU/CPU-dual-path NB math SpaNorm's own
fitting uses, rather than duplicating it.

## Usage

``` r
dnbinom_gpu(x, mu, size, log = FALSE)
```

## Arguments

- x:

  quantiles (matrix or torch tensor).

- mu:

  the mean (matrix or torch tensor, same shape as `x`).

- size:

  the NB size parameter (`1/dispersion`); a length-`nrow(x)` vector
  broadcasts row-wise, matching how
  [`stats::dnbinom`](https://rdrr.io/r/stats/NegBinomial.html) recycles
  it.

- log:

  a logical, return the log-density (default `FALSE`).

## Value

the (log-)density, as a matrix or torch tensor matching the input.

## Examples

``` r
dnbinom_gpu(matrix(0:3, 2, 2), mu = matrix(1, 2, 2), size = c(2, 2))
#>           [,1]       [,2]
#> [1,] 0.4444444 0.14814815
#> [2,] 0.2962963 0.06584362
```
