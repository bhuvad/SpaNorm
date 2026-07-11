# Invert a symmetric positive-definite matrix

Inverts a (typically small) symmetric positive-definite matrix using a
Cholesky factorisation, falling back to a general solve when the
Cholesky fails. When a torch backend is active and `mat` is a tensor,
the factorisation runs on the accelerator (or CPU for MPS) and the
result is returned on the input device. Exposed for downstream packages
(e.g. spiDE) that reuse SpaNorm's fitting machinery for Wald-type
inference.

## Usage

``` r
invert_mat(mat)
```

## Arguments

- mat:

  a symmetric positive-definite matrix (base matrix or torch tensor).

## Value

the matrix inverse, in the same representation as `mat`.

## Examples

``` r
m <- crossprod(matrix(rnorm(30), 10, 3))
invert_mat(m)
#>             [,1]        [,2]        [,3]
#> [1,]  0.06695051 -0.01098789  0.05112750
#> [2,] -0.01098789  0.18775345 -0.01322082
#> [3,]  0.05112750 -0.01322082  0.24485474
```
