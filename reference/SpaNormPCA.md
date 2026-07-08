# GLM-based (SpaNorm) PCA

GLM-based PCA using the SpaNorm model. The null model is considered to
consist of the library size effects, batch effects, and the gene mean.
GLM-PCA is approximated by regressing the null model from the data, and
performing PCA on the residuals (Pearson or deviance).

## Usage

``` r
SpaNormPCA(
  spe,
  nsvgs = 3000,
  ncomponents = 50,
  svg.fdr = 1,
  BSPARAM = bsparam(),
  BPPARAM = SerialParam(),
  residuals = c("deviance", "pearson"),
  name = "PCA"
)

# S4 method for class 'SpatialExperiment'
SpaNormPCA(
  spe,
  nsvgs = 3000,
  ncomponents = 50,
  svg.fdr = 1,
  BSPARAM = bsparam(),
  BPPARAM = SerialParam(),
  residuals = c("deviance", "pearson"),
  name = "PCA"
)

# S4 method for class 'Seurat'
SpaNormPCA(
  spe,
  nsvgs = 3000,
  ncomponents = 50,
  svg.fdr = 1,
  BSPARAM = bsparam(),
  BPPARAM = SerialParam(),
  residuals = c("deviance", "pearson"),
  name = "PCA"
)
```

## Arguments

- spe:

  a SpatialExperiment or Seurat object, with the count data stored in
  'counts' or 'data' assays respectively, and a SpaNorm model fit.

- nsvgs:

  the number of SVGs to use for PCA.

- ncomponents:

  the number of components to compute.

- svg.fdr:

  the FDR threshold for SVG calling.

- BSPARAM:

  a BiocSingularParam object specifying which algorithm should be used
  to perform the PCA.

- BPPARAM:

  a BiocParallelParam object specifying whether the PCA should be
  parallelized.

- residuals:

  the type of residuals to use for PCA. Either "deviance" (default) or
  "pearson".

- name:

  the name of the reducedDim to store the PCA results.

## Value

a SpatialExperiment or Seurat object with PCA results. For
SpatialExperiment objects, these are stored in the reducedDims.

## Details

SpaNorm PCA works by using the SpaNorm model fit for data normalisation
to approximate a GLM-based PCA as described in Townes et al. (Genome
Biology, 2019). The model used for normalisation represents the library
size effects and the gene mean. Regressing these covariates, we remain
with the deviance or Pearson residuals, upon which PCA can be performed
to approximate the GLM-PCA.

## Examples

``` r

library(SpatialExperiment)
#> Loading required package: SingleCellExperiment
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
library(ggplot2)

data(HumanDLPFC)

HumanDLPFC = SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#> (1/2) Fitting SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1150766.971816
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1150766.971816
#> iter:  1, iter:  2, log-likelihood: -818097.360502
#> iter:  1, iter:  3, log-likelihood: -730318.790270
#> iter:  1, iter:  4, log-likelihood: -715481.292625
#> iter:  1, iter:  5, log-likelihood: -713421.842805
#> iter:  1, iter:  6, log-likelihood: -713066.396466
#> iter:  1, iter:  7, log-likelihood: -712983.525523
#> iter:  1, iter:  8, log-likelihood: -712958.763665 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -712615.885064
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -712615.885064
#> iter:  2, iter:  2, log-likelihood: -712458.388066
#> iter:  2, iter:  3, log-likelihood: -712449.412685 (converged)
#> iter:  3, log-likelihood: -712449.412685 (converged)
#> (2/2) Normalising data
HumanDLPFC = SpaNormSVG(HumanDLPFC)
#> (1/3) Retrieving SpaNorm model
#> (2/3) Fitting Null SpaNorm model
#> 201 cells/spots sampled to fit model
#> iter:  1, estimating gene-wise dispersion
#> iter:  1, log-likelihood: -1150766.971816
#> iter:  1, fitting NB model
#> iter:  1, iter:  1, log-likelihood: -1150766.971816
#> iter:  1, iter:  2, log-likelihood: -818404.573175
#> iter:  1, iter:  3, log-likelihood: -736635.867673
#> iter:  1, iter:  4, log-likelihood: -723772.272788
#> iter:  1, iter:  5, log-likelihood: -722451.327139
#> iter:  1, iter:  6, log-likelihood: -722434.308242 (converged)
#> iter:  2, estimating gene-wise dispersion
#> iter:  2, log-likelihood: -722420.428878
#> iter:  2, fitting NB model
#> iter:  2, iter:  1, log-likelihood: -722420.428878
#> iter:  2, iter:  1, log-likelihood: -722420.428878
#> iter:  2, iter:  1, log-likelihood: -722420.428878
#> iter:  2, iter:  2, log-likelihood: -722420.428878
#> iter:  2, iter:  2, log-likelihood: -722420.428878
#> iter:  2, iter:  2, log-likelihood: -722420.428878
#> iter:  2, iter:  3, log-likelihood: -722420.428878 (converged)
#> iter:  3, log-likelihood: -722420.428878 (converged)
#> (3/3) Finding SVGs
#> 1245 SVGs found (FDR < 0.05)
HumanDLPFC = SpaNormPCA(HumanDLPFC)
reducedDims(HumanDLPFC)
#> List of length 1
#> names(1): PCA
```
