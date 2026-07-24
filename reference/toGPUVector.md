# Convert a vector to a torch tensor on the active backend

As
[`toGPUMatrix`](https://bhuvad.github.io/spaNorm/reference/toGPUMatrix.md),
for a numeric vector; optionally broadcasts a length-1 `vec` to length
`n`.

## Usage

``` r
toGPUVector(vec, n = NULL, backend = c("auto", "cpu", "gpu"))
```

## Arguments

- vec:

  a numeric vector (or an existing torch tensor, returned as-is).

- n:

  if supplied, the target length; a length-1 `vec` is broadcast to
  length `n`.

- backend:

  one of `"auto"` (default), `"cpu"` or `"gpu"`.

## Value

a torch tensor, or `vec` unchanged.

## Examples

``` r
toGPUVector(1:4, backend = "cpu")
#> [1] 1 2 3 4
```
