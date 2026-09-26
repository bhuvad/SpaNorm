# Does a SpaNormFit carry a converged per-gene polish?

Does a SpaNormFit carry a converged per-gene polish?

## Usage

``` r
isPolished(fit)
```

## Arguments

- fit:

  an object of class SpaNormFit.

## Value

`TRUE` iff `fit` has a non-empty `polish` slot, i.e. it was returned by
[`polishSpaNorm()`](https://bhuvad.github.io/spaNorm/reference/polishSpaNorm.md)
(or otherwise had its `polish` slot set). `FALSE` for a fresh
`fitSpaNorm()`/[`SpaNorm()`](https://bhuvad.github.io/spaNorm/reference/SpaNorm.md)
result, and for an object saved before this slot existed (checked via
[`methods::.hasSlot()`](https://rdrr.io/r/methods/slot.html), since
`fit$polish`/`fit@polish` would error on such an object).

## Examples

``` r
data(HumanDLPFC)
# \donttest{
spe <- SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#> (1/2) Fitting SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1185444.215197
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1185444.215197
#> iter:  1, iter:  2, log-likelihood: -839781.689791
#> iter:  1, iter:  3, log-likelihood: -746711.491368
#> iter:  1, iter:  4, log-likelihood: -730472.651002
#> iter:  1, iter:  5, log-likelihood: -728064.170265
#> iter:  1, iter:  6, log-likelihood: -727668.850714
#> iter:  1, iter:  7, log-likelihood: -727584.876497
#> iter:  1, iter:  8, log-likelihood: -727560.469215 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -727322.375691
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -727322.375691
#> iter:  2, iter:  2, log-likelihood: -727197.061571
#> iter:  2, iter:  3, log-likelihood: -727192.389828 (converged)
#> iter:  3, log-likelihood: -727192.389828 (converged)
#> (2/2) Normalising data
isPolished(S4Vectors::metadata(spe)$SpaNorm)
#> [1] FALSE
# }
```
