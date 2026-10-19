# Convert a matrix to a torch tensor on the active backend

Converts a base R matrix to a torch tensor on the resolved backend
device and dtype
([`getBackendDevice`](https://bhuvalab.github.io/SpaNorm/reference/getBackendDevice.md)/[`getBackendDtype`](https://bhuvalab.github.io/SpaNorm/reference/getBackendDtype.md))
if an accelerator is available and requested; otherwise returns `mat`
unchanged (already a tensor, or `backend = "cpu"`, or no accelerator
present). Exposed so downstream packages (e.g. spiDE) share the exact
same device/dtype resolution SpaNorm's own fitting uses.

## Usage

``` r
toGPUMatrix(mat, ..., backend = c("auto", "cpu", "gpu"))
```

## Arguments

- mat:

  a base R matrix (or an existing torch tensor, returned as-is).

- ...:

  ignored.

- backend:

  one of `"auto"` (default), `"cpu"` or `"gpu"`.

## Value

a torch tensor, or `mat` unchanged.

## Examples

``` r
toGPUMatrix(matrix(1:4, 2, 2), backend = "cpu")
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
```
