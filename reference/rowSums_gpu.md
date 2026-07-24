# Row sums, GPU-aware

As [`rowSums2`](https://rdrr.io/pkg/matrixStats/man/rowSums2.html), but
runs on the accelerator when `mat` is a torch tensor and an accelerator
is active.

## Usage

``` r
rowSums_gpu(mat)
```

## Arguments

- mat:

  a matrix or torch tensor.

## Value

a numeric vector (or torch tensor) of row sums.

## Examples

``` r
rowSums_gpu(matrix(1:6, 2, 3))
#> [1]  9 12
```
