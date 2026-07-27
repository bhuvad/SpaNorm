# Set the accelerator memory budget for blocked fitting

Explicitly sets the session-cached GPU memory budget (see
[`getGPUMemoryBudget`](https://bhuvad.github.io/spaNorm/reference/getGPUMemoryBudget.md)),
bypassing auto-detection for the rest of the session, or until
[`resetGPUCache`](https://bhuvad.github.io/spaNorm/reference/resetGPUCache.md)
clears it. Useful when auto-detection is unreliable – e.g. on a
MIG-partitioned GPU, where `nvidia-smi`-based detection may not resolve
the correct instance on every driver/scheduler combination – so the
correct budget (e.g. read from `nvidia-smi -L`'s reported MIG instance
size) can be set once rather than passed via `gpu.mem.budget` on every
call.

## Usage

``` r
setGPUMemoryBudget(bytes)
```

## Arguments

- bytes:

  a single positive number (bytes), or `Inf` to disable blocking
  outright.

## Value

`NULL`, invisibly.

## Examples

``` r
setGPUMemoryBudget(1e9)
getGPUMemoryBudget()
#> [1] 1e+09
resetGPUCache()
```
