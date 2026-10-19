# Is an accelerated (GPU) torch device in use?

Is an accelerated (GPU) torch device in use?

## Usage

``` r
checkGPU()
```

## Value

a logical, `TRUE` iff
[`getBackendDevice()`](https://bhuvalab.github.io/SpaNorm/reference/getBackendDevice.md)
resolves to a non-`"cpu"` device.

## Examples

``` r
checkGPU()
#> [1] FALSE
```
