# Bytes per element for the active backend dtype

The element width (in bytes) of
[`getBackendDtype()`](https://bhuvad.github.io/spaNorm/reference/getBackendDtype.md)'s
dtype – 4 for float32 (MPS), 8 for float64 (cuda/cpu). Used for sizing
memory budgets against genes x cells tensors; exposed so downstream
packages (e.g. spiDE) doing their own GPU memory-budget accounting don't
need to duplicate this device/dtype mapping.

## Usage

``` r
gpuDtypeBytes()
```

## Value

a numeric, 4 or 8.

## Examples

``` r
gpuDtypeBytes()
#> [1] 8
```
