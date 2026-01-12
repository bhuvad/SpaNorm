# Filter genes based on expression

This function removes genes that are very lowly expressed.

## Usage

``` r
filterGenes(spe, prop = 0.1)

# S4 method for class 'SpatialExperiment'
filterGenes(spe, prop = 0.1)

# S4 method for class 'Seurat'
filterGenes(spe, prop = 0.1)
```

## Arguments

- spe:

  a SpatialExperiment, Seurat, or SpatialFeatureExperiment object
  containing count data.

- prop:

  a numeric, indicating the proportion of loci/cells where the gene
  should be expressed (default is 0.1, i.e., genes should be expressed
  in at least 10% of the loci/cells).

## Value

a logical vector encoding which genes should be kept for further
analysis.

## Examples

``` r
data(HumanDLPFC)
keep <- filterGenes(HumanDLPFC)
table(keep)
#> keep
#> TRUE 
#> 5076 
```
