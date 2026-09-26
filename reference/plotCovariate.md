# Diagnostic plot of predicted expression for a covariate

This function can be used to spatially visualise the library size,
biology or batch specific effect modelled for each gene.

## Usage

``` r
plotCovariate(spe, covariate = c("biology", "ls", "batch"), ...)
```

## Arguments

- spe:

  a SpatialExperiment object.

- covariate:

  a character, specifying the type of covariate to be plot: "biology"
  (default), "ls" to plot the library size effect, and "batch" to plot
  the batch-specific effect.

- ...:

  additional parameters to be passed to the
  [plotSpatial](https://bhuvad.github.io/spaNorm/reference/plotSpatial.md)
  function.

## Value

a ggplot2 object

## Examples

``` r
library(SpatialExperiment)
library(ggplot2)

data(HumanDLPFC)
# \donttest{
HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#> (1/2) Fitting SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1126007.520596
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1126007.520596
#> iter:  1, iter:  2, log-likelihood: -797161.862744
#> iter:  1, iter:  3, log-likelihood: -714170.013420
#> iter:  1, iter:  4, log-likelihood: -700816.489950
#> iter:  1, iter:  5, log-likelihood: -699099.951708
#> iter:  1, iter:  6, log-likelihood: -698800.476586
#> iter:  1, iter:  7, log-likelihood: -698727.814555
#> iter:  1, iter:  8, log-likelihood: -698706.008926 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -698391.055515
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -698391.055515
#> iter:  2, iter:  2, log-likelihood: -698261.599091
#> iter:  2, iter:  3, log-likelihood: -698254.283587 (converged)
#> iter:  3, log-likelihood: -698254.283587 (converged)
#> (2/2) Normalising data
# plot spatial region annotations
p1 <- plotCovariate(HumanDLPFC, covariate = "biology", colour = ENSG00000075624) +
  scale_colour_viridis_c(option = "F")
p1


p2 <- plotCovariate(HumanDLPFC, covariate = "ls", colour = ENSG00000075624) +
  scale_colour_viridis_c(option = "F")
p2

# }
```
