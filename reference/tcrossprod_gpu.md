# Tensor-aware `tcrossprod`

Computes `x %*% t(y)` (or `x %*% t(x)` if `y` is `NULL`) using the
accelerator if either argument is a torch tensor and an accelerator is
active; falls back to
[`tcrossprod`](https://rdrr.io/r/base/crossprod.html) otherwise. Exposed
so downstream packages (e.g. spiDE) can build their own GPU-blocked
linear algebra with the same dual (tensor-or-matrix) semantics SpaNorm's
own fitting uses.

## Usage

``` r
tcrossprod_gpu(x, y = NULL)
```

## Arguments

- x:

  a matrix or torch tensor.

- y:

  a matrix or torch tensor, or `NULL` (default).

## Value

`x %*% t(y)`, as a matrix or torch tensor matching the input.

## Examples

``` r
m <- matrix(rnorm(12), nrow = 3)
tcrossprod_gpu(m)
#>            [,1]       [,2]       [,3]
#> [1,]  4.0193579 -0.9301318  0.7354387
#> [2,] -0.9301318  0.4280649 -0.2297178
#> [3,]  0.7354387 -0.2297178  1.9126805
```
