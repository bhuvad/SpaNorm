# The floor on a fitted negative binomial mean

The polish engine floors every fitted mean at this value, and a caller
that rebuilds the mean from polished coefficients (for inference, say)
should floor at the same value, so both stages work on one mean
function. It is a numerical guard against
[`exp()`](https://rdrr.io/r/base/Log.html) underflowing to exactly zero,
not a statistical clamp: `exp(-30)` is far below any mean the model can
meaningfully estimate.

## Usage

``` r
nbMuFloor()
```

## Value

`exp(-30)`.

## Examples

``` r
nbMuFloor()
#> [1] 9.357623e-14
```
