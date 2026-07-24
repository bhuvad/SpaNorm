# Does an object hold any NA/NaN/Inf?

`TRUE` if `x` (tensor or matrix) holds any `NA`/`NaN`/ `Inf`. On a
tensor this is a scalar reduction copied back, not a full GPU-\>CPU
materialisation of `x`.

## Usage

``` r
hasBadValues(x)
```

## Arguments

- x:

  a matrix or torch tensor.

## Value

a logical.

## Examples

``` r
hasBadValues(matrix(c(1, NA, 3, 4), 2, 2))
#> [1] TRUE
```
