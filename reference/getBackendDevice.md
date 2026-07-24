# Resolve the active torch backend device

Resolves (once per session, then cached – see
[`resetGPUCache`](https://bhuvad.github.io/spaNorm/reference/resetGPUCache.md))
the torch device to use for accelerated operations: `"cuda"` if a CUDA
GPU is available, else `"mps"` on Apple Silicon, else `"cpu"`. Exposed
for downstream packages (e.g. spiDE) that push their own tensors onto
the same device SpaNorm is using.

## Usage

``` r
getBackendDevice()
```

## Value

a character, one of `"cuda"`, `"mps"` or `"cpu"`.

## Examples

``` r
getBackendDevice()
#> [1] "cpu"
```
