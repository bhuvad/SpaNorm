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
HumanDLPFC = SpaNormSVG(HumanDLPFC)
#> (1/3) Retrieving SpaNorm model
#> (2/3) Fitting Null SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1126007.520596
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1126007.520596
#> iter:  1, iter:  2, log-likelihood: -799092.515118
#> iter:  1, iter:  3, log-likelihood: -720678.671068
#> iter:  1, iter:  4, log-likelihood: -708864.101492
#> iter:  1, iter:  5, log-likelihood: -707693.635599
#> iter:  1, iter:  6, log-likelihood: -707677.962961 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -707642.791806
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -707642.791806
#> iter:  2, iter:  1, log-likelihood: -707642.791806
#> iter:  2, iter:  1, log-likelihood: -707642.791806
#> iter:  2, iter:  2, log-likelihood: -707642.791806
#> iter:  2, iter:  2, log-likelihood: -707642.791806
#> iter:  2, iter:  2, log-likelihood: -707642.791806
#> iter:  2, iter:  3, log-likelihood: -707642.791806 (converged)
#> iter:  3, log-likelihood: -707642.791806 (converged)
#> (3/3) Finding SVGs
#> 1348 SVGs found (FDR < 0.05)
topSVGs = topSVGs(HumanDLPFC, n = 10)
```
