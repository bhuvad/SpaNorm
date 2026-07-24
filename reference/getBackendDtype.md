# Resolve the active torch backend dtype

The torch floating-point dtype to use for accelerated tensors: float64
everywhere except MPS, which cannot represent float64 and so is capped
at float32 (its maximum precision).

## Usage

``` r
getBackendDtype()
```

## Value

a torch dtype
([`torch::torch_float64()`](https://torch.mlverse.org/docs/reference/torch_dtype.html)
or
[`torch::torch_float32()`](https://torch.mlverse.org/docs/reference/torch_dtype.html)).

## Examples

``` r
if (checkGPU()) getBackendDtype()
```
