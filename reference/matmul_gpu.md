# Tensor-aware matrix multiplication

Computes `x %*% y` using the accelerator if either argument is a torch
tensor and an accelerator is active (vectors are reshaped so the result
matches base R's `%*%` shape semantics); falls back to base `%*%`
otherwise. Named explicitly (rather than overriding `%*%`) so base
matrix multiplication semantics are preserved everywhere else. Exposed
so downstream packages (e.g. spiDE) can build batched linear algebra
(e.g. a shared design matrix against many per-gene/per-block weight
vectors) without depending on torch internals.

## Usage

``` r
matmul_gpu(x, y)
```

## Arguments

- x, y:

  a matrix, numeric vector, or torch tensor.

## Value

`x %*% y`, as a matrix or torch tensor matching the input.

## Examples

``` r
matmul_gpu(matrix(1:6, 2, 3), matrix(1:6, 3, 2))
#>      [,1] [,2]
#> [1,]   22   49
#> [2,]   28   64
```
