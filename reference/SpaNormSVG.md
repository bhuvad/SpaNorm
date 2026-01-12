# Model-based spatially variable gene (SVG) calling

Spatially variable gene (SVG) calling using the SpaNorm model.

## Usage

``` r
SpaNormSVG(spe, backend = c("auto", "cpu", "gpu"), verbose = TRUE)

# S4 method for class 'SpatialExperiment'
SpaNormSVG(spe, backend = c("auto", "cpu", "gpu"), verbose = TRUE)
```

## Arguments

- spe:

  a SpatialExperiment or Seurat object, with the count data stored in
  'counts' or 'data' assays respectively, and a SpaNorm model fit.

- backend:

  a character, specifying the backend to use for computations. Options
  are "auto" (default), "cpu", or "gpu". If "auto", it will use GPU if
  available, otherwise CPU.

- verbose:

  a logical, specifying whether to show update messages (default TRUE).

## Value

a SpatialExperiment or Seurat object with F-statistics, false discovery
rates (FDRs). For SpatialExperiment objects, these are stored in the
rowData.

## Details

SpaNorm SVG calling works by using the SpaNorm model fit for data
normalisation to perform a likelihood ratio test (LRT). The model used
for normalisation is considered to be the full model. A second nested
model is fit without the splines representing biology. These nested
models are then compared using a LRT to identify genes where the splines
representing biology contain strong signal.

## Examples

``` r
library(SpatialExperiment)
library(ggplot2)

data(HumanDLPFC)

HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#> (1/2) Fitting SpaNorm model
#> 201 cells/spots sampled to fit model
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, estimating gene-wise dispersion
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, log-likelihood: -1182732.284681
#> iter:  1, fitting NB model
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  1, log-likelihood: -1182732.284681
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  2, log-likelihood: -834607.897784
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  3, log-likelihood: -743087.131350
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  4, log-likelihood: -727106.108939
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  5, log-likelihood: -724746.362495
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  6, log-likelihood: -724375.196568
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  7, log-likelihood: -724300.281954
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  8, log-likelihood: -724279.738344 (converged)
#> iter:  2, estimating gene-wise dispersion
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, log-likelihood: -723961.820395
#> iter:  2, fitting NB model
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  1, log-likelihood: -723961.820395
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  2, log-likelihood: -723841.634557
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  3, log-likelihood: -723834.091681 (converged)
#> iter:  3, log-likelihood: -723834.091681 (converged)
#> (2/2) Normalising data
HumanDLPFC = SpaNormSVG(HumanDLPFC)
#> (1/3) Retrieving SpaNorm model
#> (2/3) Fitting Null SpaNorm model
#> 201 cells/spots sampled to fit model
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, estimating gene-wise dispersion
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, log-likelihood: -1182732.284681
#> iter:  1, fitting NB model
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  1, log-likelihood: -1182732.284681
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  2, log-likelihood: -837874.274131
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  3, log-likelihood: -750511.292381
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  4, log-likelihood: -735776.020158
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  5, log-likelihood: -733842.363539
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  6, log-likelihood: -733682.885089
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  6, log-likelihood: -733682.885089
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  6, log-likelihood: -733682.885089
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  7, log-likelihood: -733682.885089 (converged)
#> iter:  2, estimating gene-wise dispersion
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, log-likelihood: -733655.276018
#> iter:  2, fitting NB model
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  1, log-likelihood: -733655.276018
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  1, log-likelihood: -733655.276018
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  1, log-likelihood: -733655.276018
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  2, log-likelihood: -733655.276018
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  2, log-likelihood: -733655.276018
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  2, log-likelihood: -733655.276018
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  3, log-likelihood: -733655.276018 (converged)
#> iter:  3, log-likelihood: -733655.276018 (converged)
#> (3/3) Finding SVGs
#> 1285 SVGs found (FDR < 0.05)
head(rowData(HumanDLPFC))
#> DataFrame with 6 rows and 5 columns
#>                   gene_name   gene_biotype     svg.F       svg.p     svg.fdr
#>                 <character>    <character> <numeric>   <numeric>   <numeric>
#> ENSG00000188976       NOC2L protein_coding   0.00000 1.00000e+00 1.00000e+00
#> ENSG00000188290        HES4 protein_coding  24.58751 3.90378e-20 1.04844e-18
#> ENSG00000187608       ISG15 protein_coding   0.00000 1.00000e+00 1.00000e+00
#> ENSG00000188157        AGRN protein_coding   4.22614 2.04571e-03 9.23845e-03
#> ENSG00000078808        SDF4 protein_coding   0.00000 1.00000e+00 1.00000e+00
#> ENSG00000176022     B3GALT6 protein_coding   0.00000 1.00000e+00 1.00000e+00
```
