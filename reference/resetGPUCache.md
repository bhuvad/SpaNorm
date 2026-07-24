# Clear the cached GPU device/memory-budget state

Clears the session-level cache of the resolved torch device and memory
budget (set by
[`getBackendDevice`](https://bhuvad.github.io/spaNorm/reference/getBackendDevice.md)/[`getGPUMemoryBudget`](https://bhuvad.github.io/spaNorm/reference/getGPUMemoryBudget.md)).
Useful in tests, or if the available device/memory changes within a
session.

## Usage

``` r
resetGPUCache()
```

## Value

`NULL`, invisibly.

## Examples

``` r
resetGPUCache()
```
