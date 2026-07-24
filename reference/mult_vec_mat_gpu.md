# Broadcast-multiply a vector and a matrix, GPU-aware

As
[`add_vec_mat_gpu`](https://bhuvad.github.io/spaNorm/reference/add_vec_mat_gpu.md),
but multiplying rather than adding.

## Usage

``` r
mult_vec_mat_gpu(vec, mat, backend = c("auto", "cpu", "gpu"))
```

## Arguments

- vec:

  a numeric vector (or torch tensor).

- mat:

  a matrix (or torch tensor).

- backend:

  one of `"auto"` (default), `"cpu"` or `"gpu"`.

## Value

`mat * vec`, broadcast as per
[`add_vec_mat_gpu`](https://bhuvad.github.io/spaNorm/reference/add_vec_mat_gpu.md).

## Examples

``` r
mult_vec_mat_gpu(1:3, matrix(1, 3, 4), backend = "cpu")
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    1    1    1
#> [2,]    2    2    2    2
#> [3,]    3    3    3    3
```
