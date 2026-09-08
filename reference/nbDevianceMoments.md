# First two moments of the NB unit deviance, by direct summation

Reproduces the definition edgeR's `compute_weight()` approximates with
Chebyshev tables. The summation window is centred on `mu` and sized from
the pmf's own quantiles rather than fixed at edgeR's 50 terms from zero:
edgeR can truncate there because that branch is only reached when
`phi >= 4`, where the pmf piles up near zero; this is reached with any
dispersion.

## Usage

``` r
nbDevianceMoments(mu, phi, eps = 1e-10, maxterms = 20000L)
```

## Arguments

- mu:

  vector of means.

- phi:

  NB dispersion: scalar, or the same length as `mu`.

- eps:

  tail probability left outside the summation window at each end.

- maxterms:

  hard cap on the number of pmf terms summed per element.

## Value

list with `w0` (deviance rescaling) and `w1` (effective df), each the
length of `mu`.

## Details

Elements are grouped by the width they need so a few high-mu cells do
not impose their window on everything, and each group accumulates in two
passes so memory stays O(length(mu)) rather than O(width \* length(mu)).

## Examples

``` r
nbDevianceMoments(mu = c(0.1, 1, 10), phi = 5)
#> $w0
#> [1] 1.788589 4.251873 3.001213
#> 
#> $w1
#> [1] 0.5233816 2.8091676 2.8650045
#> 
```
