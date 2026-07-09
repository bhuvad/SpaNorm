# Export top SVG results to a data frame

Export top SVG results to a data frame

## Usage

``` r
topSVGs(spe, n = 10, fdr = 1)
```

## Arguments

- spe:

  a SpatialExperiment object with SVG results from SpaNormSVG.

- n:

  a numeric, specifying the number of top SVGs to call.

- fdr:

  a numeric, specifying the false discovery rate (FDR) threshold for
  calling SVGs.

## Value

A data frame containing the top SVGs from F-test results including
F-statistics, p-values and FDR.

## Examples

``` r

library(SpatialExperiment)
library(ggplot2)

data(HumanDLPFC)

HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#> (1/2) Fitting SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1118949.014392
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1118949.014392
#> iter:  1, iter:  2, log-likelihood: -798368.531434
#> iter:  1, iter:  3, log-likelihood: -715073.474205
#> iter:  1, iter:  4, log-likelihood: -700895.334178
#> iter:  1, iter:  5, log-likelihood: -698953.831943
#> iter:  1, iter:  6, log-likelihood: -698599.754859
#> iter:  1, iter:  7, log-likelihood: -698510.763916
#> iter:  1, iter:  8, log-likelihood: -698484.373174 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -698242.639330
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -698242.639330
#> iter:  2, iter:  2, log-likelihood: -698117.656127
#> iter:  2, iter:  3, log-likelihood: -698108.879713 (converged)
#> iter:  3, log-likelihood: -698108.879713 (converged)
#> (2/2) Normalising data
HumanDLPFC = SpaNormSVG(HumanDLPFC)
#> (1/3) Retrieving SpaNorm model
#> (2/3) Fitting Null SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1118949.014392
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1118949.014392
#> iter:  1, iter:  2, log-likelihood: -796647.625266
#> iter:  1, iter:  3, log-likelihood: -719440.482083
#> iter:  1, iter:  4, log-likelihood: -707998.494771
#> iter:  1, iter:  5, log-likelihood: -706919.057060
#> iter:  1, iter:  6, log-likelihood: -706895.545484 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -706866.664769
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -706866.664769
#> iter:  2, iter:  1, log-likelihood: -706866.664769
#> iter:  2, iter:  1, log-likelihood: -706866.664769
#> iter:  2, iter:  2, log-likelihood: -706866.664769
#> iter:  2, iter:  2, log-likelihood: -706866.664769
#> iter:  2, iter:  2, log-likelihood: -706866.664769
#> iter:  2, iter:  3, log-likelihood: -706866.664769 (converged)
#> iter:  3, log-likelihood: -706866.664769 (converged)
#> (3/3) Finding SVGs
#> 1534 SVGs found (FDR < 0.05)
topSVGs = topSVGs(HumanDLPFC, n = 10)
```
