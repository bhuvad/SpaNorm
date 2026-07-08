# Model-based spatially variable gene (SVG) calling

Spatially variable gene (SVG) calling using the SpaNorm model.

## Usage

``` r
SpaNormSVG(spe, backend = c("auto", "cpu", "gpu"), verbose = TRUE)

# S4 method for class 'SpatialExperiment'
SpaNormSVG(spe, backend = c("auto", "cpu", "gpu"), verbose = TRUE)

# S4 method for class 'Seurat'
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
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1188651.839051
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1188651.839051
#> iter:  1, iter:  2, log-likelihood: -838534.577238
#> iter:  1, iter:  3, log-likelihood: -746271.504655
#> iter:  1, iter:  4, log-likelihood: -729897.735928
#> iter:  1, iter:  5, log-likelihood: -727408.241156
#> iter:  1, iter:  6, log-likelihood: -727021.233884
#> iter:  1, iter:  7, log-likelihood: -726948.410401
#> iter:  1, iter:  8, log-likelihood: -726929.357994 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -726617.561752
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -726617.561752
#> iter:  2, iter:  2, log-likelihood: -726494.458843
#> iter:  2, iter:  3, log-likelihood: -726488.535634 (converged)
#> iter:  3, log-likelihood: -726488.535634 (converged)
#> (2/2) Normalising data
HumanDLPFC = SpaNormSVG(HumanDLPFC)
#> (1/3) Retrieving SpaNorm model
#> (2/3) Fitting Null SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1188651.839051
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1188651.839051
#> iter:  1, iter:  2, log-likelihood: -842378.925238
#> iter:  1, iter:  3, log-likelihood: -753995.514810
#> iter:  1, iter:  4, log-likelihood: -738739.077835
#> iter:  1, iter:  5, log-likelihood: -736640.376089
#> iter:  1, iter:  6, log-likelihood: -736458.772822
#> iter:  1, iter:  6, log-likelihood: -736458.772822
#> iter:  1, iter:  6, log-likelihood: -736458.772822
#> iter:  1, iter:  7, log-likelihood: -736458.772822 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -736440.814582
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -736440.814582
#> iter:  2, iter:  1, log-likelihood: -736440.814582
#> iter:  2, iter:  1, log-likelihood: -736440.814582
#> iter:  2, iter:  2, log-likelihood: -736440.814582
#> iter:  2, iter:  2, log-likelihood: -736440.814582
#> iter:  2, iter:  2, log-likelihood: -736440.814582
#> iter:  2, iter:  3, log-likelihood: -736440.814582 (converged)
#> iter:  3, log-likelihood: -736440.814582 (converged)
#> (3/3) Finding SVGs
#> 1237 SVGs found (FDR < 0.05)
head(rowData(HumanDLPFC))
#> DataFrame with 6 rows and 5 columns
#>                   gene_name   gene_biotype     svg.F       svg.p    svg.fdr
#>                 <character>    <character> <numeric>   <numeric>  <numeric>
#> ENSG00000188976       NOC2L protein_coding   0.00000 1.00000e+00 1.0000e+00
#> ENSG00000188290        HES4 protein_coding  17.99709 1.16743e-14 1.9493e-13
#> ENSG00000187608       ISG15 protein_coding   0.00000 1.00000e+00 1.0000e+00
#> ENSG00000188157        AGRN protein_coding   2.96211 1.86299e-02 7.3764e-02
#> ENSG00000078808        SDF4 protein_coding   0.00000 1.00000e+00 1.0000e+00
#> ENSG00000176022     B3GALT6 protein_coding   0.00000 1.00000e+00 1.0000e+00
```
