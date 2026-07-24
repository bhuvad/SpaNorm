# Convert a torch tensor (or R object) to a base R matrix

Converts a torch tensor to a base R matrix (moving it to the CPU first
if needed); non-tensor input is coerced via
[`as.matrix`](https://rdrr.io/r/base/matrix.html). Exposed so downstream
packages (e.g. spiDE) doing their own GPU-blocked computation can
convert results back to plain R without depending on torch internals.

## Usage

``` r
toRMatrix(x)
```

## Arguments

- x:

  a torch tensor or an R object coercible via `as.matrix`.

## Value

a base R matrix.

## Examples

``` r
toRMatrix(matrix(1:4, 2, 2))
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
```
