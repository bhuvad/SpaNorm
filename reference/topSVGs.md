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
#> iter:  1, log-likelihood: -1119713.505106
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1119713.505106
#> iter:  1, iter:  2, log-likelihood: -798605.721672
#> iter:  1, iter:  3, log-likelihood: -715923.179376
#> iter:  1, iter:  4, log-likelihood: -702159.496628
#> iter:  1, iter:  5, log-likelihood: -700328.159628
#> iter:  1, iter:  6, log-likelihood: -700015.031384
#> iter:  1, iter:  7, log-likelihood: -699941.719733
#> iter:  1, iter:  8, log-likelihood: -699919.137466 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -699611.667949
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -699611.667949
#> iter:  2, iter:  2, log-likelihood: -699479.739050
#> iter:  2, iter:  3, log-likelihood: -699470.613563 (converged)
#> iter:  3, log-likelihood: -699470.613563 (converged)
#> (2/2) Normalising data
HumanDLPFC = SpaNormSVG(HumanDLPFC)
#> (1/3) Retrieving SpaNorm model
#> (2/3) Fitting Null SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1119713.505106
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1119713.505106
#> iter:  1, iter:  2, log-likelihood: -799017.372419
#> iter:  1, iter:  3, log-likelihood: -721729.185444
#> iter:  1, iter:  4, log-likelihood: -710184.610200
#> iter:  1, iter:  5, log-likelihood: -709113.039844
#> iter:  1, iter:  6, log-likelihood: -709108.399010 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -709071.412745
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -709071.412745
#> iter:  2, iter:  1, log-likelihood: -709071.412745
#> iter:  2, iter:  1, log-likelihood: -709071.412745
#> iter:  2, iter:  2, log-likelihood: -709071.412745
#> iter:  2, iter:  2, log-likelihood: -709071.412745
#> iter:  2, iter:  2, log-likelihood: -709071.412745
#> iter:  2, iter:  3, log-likelihood: -709071.412745 (converged)
#> iter:  3, log-likelihood: -709071.412745 (converged)
#> (3/3) Finding SVGs
#> 1200 SVGs found (FDR < 0.05)
topSVGs = topSVGs(HumanDLPFC, n = 10)
```
