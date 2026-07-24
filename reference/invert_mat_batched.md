# Invert a batch of symmetric positive-definite matrices

Batched counterpart of
[`invert_mat`](https://bhuvad.github.io/spaNorm/reference/invert_mat.md):
inverts a `(batch, n, n)` stack of symmetric positive-definite matrices
in a single Cholesky factorisation call rather than looping
[`invert_mat()`](https://bhuvad.github.io/spaNorm/reference/invert_mat.md)
over each slice. SpaNorm's own model never needs this (its per-gene fit
shares a single coefficient across all genes, so it only ever inverts
one matrix at a time); it is exposed for downstream packages (e.g.
spiDE) that need a genuinely per-gene/per-feature covariance – a
different `n x n` matrix for every slice of the batch, varying with a
per-gene working weight – where looping
[`invert_mat()`](https://bhuvad.github.io/spaNorm/reference/invert_mat.md)'s
R-level `tryCatch`/Cholesky call once per gene is the dominant cost.

## Usage

``` r
invert_mat_batched(mat)
```

## Arguments

- mat:

  a `(batch, n, n)` torch tensor, or a base R array of the same shape
  (`dim(mat) == c(batch, n, n)`).

## Value

the batch of matrix inverses, in the same representation as `mat`.

## Details

Unlike
[`invert_mat()`](https://bhuvad.github.io/spaNorm/reference/invert_mat.md)'s
per-matrix `tryCatch`, a failure here (e.g. one non-SPD slice) falls
back to a general solve for the *whole* batch, not just the failing
slice – correctness is unaffected (the fallback solves every slice
correctly regardless of conditioning), it is simply not the exact same
numerical routine applied per-slice in that (rare) case.

## Examples

``` r
m <- array(0, c(4, 3, 3))
for (i in 1:4) m[i, , ] <- crossprod(matrix(rnorm(9), 3, 3)) + diag(3)
invert_mat_batched(m)
#> , , 1
#> 
#>           [,1]        [,2]        [,3]
#> [1,] 0.5409032  0.02995611 -0.34225119
#> [2,] 0.3085952 -0.08673537 -0.26976422
#> [3,] 0.5043533 -0.16815179 -0.08212729
#> [4,] 0.7557993  0.20616298  0.23024131
#> 
#> , , 2
#> 
#>             [,1]      [,2]         [,3]
#> [1,]  0.02995611 0.4863725 -0.040903185
#> [2,] -0.08673537 0.3418442 -0.004202047
#> [3,] -0.16815179 0.3653400 -0.058594193
#> [4,]  0.20616298 0.2549388  0.123926094
#> 
#> , , 3
#> 
#>             [,1]         [,2]      [,3]
#> [1,] -0.34225119 -0.040903185 0.6582193
#> [2,] -0.26976422 -0.004202047 0.8759369
#> [3,] -0.08212729 -0.058594193 0.2439739
#> [4,]  0.23024131  0.123926094 0.2287382
#> 
```
