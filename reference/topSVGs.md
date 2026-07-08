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
#> iter:  1, log-likelihood: -1151095.650438
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1151095.650438
#> iter:  1, iter:  2, log-likelihood: -811517.870744
#> iter:  1, iter:  3, log-likelihood: -724866.413098
#> iter:  1, iter:  4, log-likelihood: -709875.790517
#> iter:  1, iter:  5, log-likelihood: -707690.692163
#> iter:  1, iter:  6, log-likelihood: -707348.747499
#> iter:  1, iter:  7, log-likelihood: -707280.952543 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -706914.055060
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -706914.055060
#> iter:  2, iter:  2, log-likelihood: -706767.384132
#> iter:  2, iter:  3, log-likelihood: -706755.422571 (converged)
#> iter:  3, log-likelihood: -706755.422571 (converged)
#> (2/2) Normalising data
HumanDLPFC = SpaNormSVG(HumanDLPFC)
#> (1/3) Retrieving SpaNorm model
#> (2/3) Fitting Null SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1151095.650438
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1151095.650438
#> iter:  1, iter:  2, log-likelihood: -817593.526908
#> iter:  1, iter:  3, log-likelihood: -733088.659408
#> iter:  1, iter:  4, log-likelihood: -718980.135957
#> iter:  1, iter:  5, log-likelihood: -717114.990731
#> iter:  1, iter:  6, log-likelihood: -716990.831370
#> iter:  1, iter:  6, log-likelihood: -716990.831370
#> iter:  1, iter:  6, log-likelihood: -716990.831370
#> iter:  1, iter:  7, log-likelihood: -716990.831370 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -716732.056292
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -716732.056292
#> iter:  2, iter:  1, log-likelihood: -716732.056292
#> iter:  2, iter:  1, log-likelihood: -716732.056292
#> iter:  2, iter:  2, log-likelihood: -716732.056292
#> iter:  2, iter:  2, log-likelihood: -716732.056292
#> iter:  2, iter:  2, log-likelihood: -716732.056292
#> iter:  2, iter:  3, log-likelihood: -716732.056292 (converged)
#> iter:  3, estimating gene-wise dispersion
#> iter:  3, log-likelihood: -716732.056292
#> iter:  3, fitting NB model
#> iter:  3, iter:  1, log-likelihood: -716732.056292
#> iter:  3, iter:  1, log-likelihood: -716732.056292
#> iter:  3, iter:  1, log-likelihood: -716732.056292
#> iter:  3, iter:  2, log-likelihood: -716732.056292
#> iter:  3, iter:  2, log-likelihood: -716732.056292
#> iter:  3, iter:  2, log-likelihood: -716732.056292
#> iter:  3, iter:  3, log-likelihood: -716732.056292 (converged)
#> iter:  4, log-likelihood: -716732.056292 (converged)
#> (3/3) Finding SVGs
#> 1050 SVGs found (FDR < 0.05)
topSVGs = topSVGs(HumanDLPFC, n = 10)
```
