# Filter genes based on expression

This function computes the size factors using a fast but inaccurate
approach. Size factors are computed using the direct estimate of library
sizes (sum of all counts). Though fast, this approach does not cater for
compositional biases in the data and therefore is less accurate than
scran-based estimates.

## Usage

``` r
fastSizeFactors(spe)

# S4 method for class 'SpatialExperiment'
fastSizeFactors(spe)
```

## Arguments

- spe:

  a SpatialExperiment, Seurat, or SpatialFeatureExperiment object
  containing count data.

## Value

a SpatialExperiment, Seurat, or SpatialFeatureExperiment, containing
size factors in the 'sizeFactor' column of the column annotation.

## Examples

``` r
data(HumanDLPFC)
HumanDLPFC <- fastSizeFactors(HumanDLPFC)
head(HumanDLPFC$sizeFactor)
#> AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 AAACAGAGCGACTCCT-1 AAACAGCTTTCAGAAG-1 
#>          1.3854989          1.5920004          1.0583586          0.9358734 
#> AAACAGGGTCTATATT-1 AAACATTTCCCGGATT-1 
#>          1.7868073          1.0057331 
```
