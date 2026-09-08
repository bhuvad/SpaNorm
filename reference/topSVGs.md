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
#> iter:  1, log-likelihood: -1115103.280782
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1115103.280782
#> iter:  1, iter:  2, log-likelihood: -794964.478436
#> iter:  1, iter:  3, log-likelihood: -712182.411762
#> iter:  1, iter:  4, log-likelihood: -698208.242581
#> iter:  1, iter:  5, log-likelihood: -696345.622476
#> iter:  1, iter:  6, log-likelihood: -696012.591368
#> iter:  1, iter:  7, log-likelihood: -695930.346690
#> iter:  1, iter:  8, log-likelihood: -695906.099189 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -695637.573816
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -695637.573816
#> iter:  2, iter:  2, log-likelihood: -695507.572354
#> iter:  2, iter:  3, log-likelihood: -695499.852807 (converged)
#> iter:  3, log-likelihood: -695499.852807 (converged)
#> (2/2) Normalising data
HumanDLPFC = SpaNormSVG(HumanDLPFC)
#> (1/3) Retrieving SpaNorm model
#> (2/3) Fitting Null SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1115103.280782
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1115103.280782
#> iter:  1, iter:  2, log-likelihood: -794161.257817
#> iter:  1, iter:  3, log-likelihood: -717114.154534
#> iter:  1, iter:  4, log-likelihood: -705708.802500
#> iter:  1, iter:  5, log-likelihood: -704629.727873
#> iter:  1, iter:  6, log-likelihood: -704600.117394 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -704565.889641
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -704565.889641
#> iter:  2, iter:  1, log-likelihood: -704565.889641
#> iter:  2, iter:  1, log-likelihood: -704565.889641
#> iter:  2, iter:  2, log-likelihood: -704565.889641
#> iter:  2, iter:  2, log-likelihood: -704565.889641
#> iter:  2, iter:  2, log-likelihood: -704565.889641
#> iter:  2, iter:  3, log-likelihood: -704565.889641 (converged)
#> iter:  3, log-likelihood: -704565.889641 (converged)
#> (3/3) Finding SVGs
#> 1430 SVGs found (FDR < 0.05)
topSVGs = topSVGs(HumanDLPFC, n = 10)
```
