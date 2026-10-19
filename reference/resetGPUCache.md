# Clear the cached GPU device/memory-budget state

Clears the session-level cache of the resolved torch device and memory
budget (set by
[`getBackendDevice`](https://bhuvalab.github.io/SpaNorm/reference/getBackendDevice.md)/[`getGPUMemoryBudget`](https://bhuvalab.github.io/SpaNorm/reference/getGPUMemoryBudget.md)).
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
