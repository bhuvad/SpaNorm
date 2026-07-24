# Determine the accelerator memory budget for blocked fitting

Resolves (once per session, then cached – see
[`resetGPUCache`](https://bhuvad.github.io/spaNorm/reference/resetGPUCache.md))
the accelerator memory budget (bytes) to use for memory-aware
gene-blocked fitting: free CUDA memory via `nvidia-smi`, a conservative
fraction of total system memory for MPS (unified memory has no
per-process query), or a fixed fallback if detection fails. A
user-supplied `gpu.mem.budget` bypasses detection entirely – set it to
`Inf` to disable blocking outright. Returns `Inf` when no accelerator is
in use. Exposed so downstream packages (e.g. spiDE) doing their own
GPU-memory-aware blocking share the same budget/cache as SpaNorm's own
fitting.

## Usage

``` r
getGPUMemoryBudget(gpu.mem.budget = NULL)
```

## Arguments

- gpu.mem.budget:

  `NULL` (default, auto-detect) or a single positive number (bytes).

## Value

a numeric, the memory budget in bytes (`Inf` if no accelerator is in
use).

## Examples

``` r
getGPUMemoryBudget()
#> [1] Inf
```
