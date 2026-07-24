# Build a diagonal matrix from a vector, GPU-aware

As [`diag`](https://rdrr.io/r/base/diag.html), but returns a torch
tensor (dense, to match the CPU fallback's dense output) when an
accelerator is active.

## Usage

``` r
diag_mat(vec, backend = c("auto", "cpu", "gpu"))
```

## Arguments

- vec:

  a numeric vector (or torch tensor).

- backend:

  one of `"auto"` (default), `"cpu"` or `"gpu"`.

## Value

a diagonal matrix (base matrix or torch tensor).

## Examples

``` r
diag_mat(1:3, backend = "cpu")
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    2    0
#> [3,]    0    0    3
```
