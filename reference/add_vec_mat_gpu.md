# Broadcast-add a vector to a matrix, GPU-aware

Adds `vec` to `mat`, broadcasting per-row if `length(vec) == nrow(mat)`
or per-column if `length(vec) == ncol(mat)` (matching the GPU path's
broadcast orientation explicitly, rather than relying on R's
column-major recycling). Runs on the accelerator when either argument is
a torch tensor and an accelerator is active.

## Usage

``` r
add_vec_mat_gpu(vec, mat, backend = c("auto", "cpu", "gpu"))
```

## Arguments

- vec:

  a numeric vector (or torch tensor).

- mat:

  a matrix (or torch tensor).

- backend:

  one of `"auto"` (default), `"cpu"` or `"gpu"`.

## Value

`mat + vec`, broadcast as described above.

## Examples

``` r
add_vec_mat_gpu(1:3, matrix(1, 3, 4), backend = "cpu")
#>      [,1] [,2] [,3] [,4]
#> [1,]    2    2    2    2
#> [2,]    3    3    3    3
#> [3,]    4    4    4    4
```
