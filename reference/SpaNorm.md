# Spatially-dependent normalisation for spatial transcriptomics data

Performs normalisation of spatial transcriptomics data using
spatially-dependent spot- and gene- specific size factors.

## Usage

``` r
SpaNorm(
  spe,
  sample.p = 0.25,
  gene.model = c("nb"),
  adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
  scale.factor = 1,
  df.tps = 6,
  lambda.a = 1e-04,
  batch = NULL,
  tol = 1e-04,
  step.factor = 0.5,
  maxit.nb = 50,
  maxit.psi = 25,
  maxn.psi = 500,
  overwrite = FALSE,
  backend = c("auto", "cpu", "gpu"),
  verbose = TRUE,
  ...
)

# S4 method for class 'SpatialExperiment'
SpaNorm(
  spe,
  sample.p = 0.25,
  gene.model = c("nb"),
  adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
  scale.factor = 1,
  df.tps = 6,
  lambda.a = 1e-04,
  batch = NULL,
  tol = 1e-04,
  step.factor = 0.5,
  maxit.nb = 50,
  maxit.psi = 25,
  maxn.psi = 500,
  overwrite = FALSE,
  backend = c("auto", "cpu", "gpu"),
  verbose = TRUE,
  ...
)

# S4 method for class 'Seurat'
SpaNorm(
  spe,
  sample.p = 0.25,
  gene.model = c("nb"),
  adj.method = c("auto", "logpac", "pearson", "medbio", "meanbio"),
  scale.factor = 1,
  df.tps = 6,
  lambda.a = 1e-04,
  batch = NULL,
  tol = 1e-04,
  step.factor = 0.5,
  maxit.nb = 50,
  maxit.psi = 25,
  maxn.psi = 500,
  overwrite = FALSE,
  backend = c("auto", "cpu", "gpu"),
  verbose = TRUE,
  ...
)
```

## Arguments

- spe:

  a SpatialExperiment or Seurat object, with the count data stored in
  'counts' or 'data' assays respectively.

- sample.p:

  a numeric, specifying the (maximum) proportion of cells/spots to
  sample for model fitting (default is 0.25).

- gene.model:

  a character, specifying the model to use for gene/protein abundances
  (default 'nb'). This should be 'nb' for count based datasets.

- adj.method:

  a character, specifying the method to use to adjust the data (default
  'auto', see details)

- scale.factor:

  a numeric, specifying the sample-specific scaling factor to scale the
  adjusted count.

- df.tps:

  a numeric, of length 1 or 2, specifying the maximum degrees of freedom
  for the thin-plate spline for the biology and library size effects,
  respectively (default is 6, see details).

- lambda.a:

  a numeric, of length 1 or 2, specifying the smoothing parameter for
  regularizing regression coefficients (default is 0.0001, see details).
  Actual lambda.a used is lambda.a \* ncol(spe).

- batch:

  a vector or numeric matrix, specifying the batch design to regress out
  (default NULL, representing no batch effects). See details for more
  information on how to define this variable.

- tol:

  a numeric, specifying the tolerance for convergence (default is 1e-4).

- step.factor:

  a numeric, specifying the multiplicative factor to decrease IRLS step
  by when log-likelihood diverges (default is 0.5).

- maxit.nb:

  a numeric, specifying the maximum number of IRLS iteration for
  estimating NB mean parameters for a given dispersion parameter
  (default is 50).

- maxit.psi:

  a numeric, specifying the maximum number of IRLS iterations to
  estimate the dispersion parameter (default is 25).

- maxn.psi:

  a numeric, specifying the maximum number of cells/spots to sample for
  dispersion estimation (default is 500).

- overwrite:

  a logical, specifying whether to force recomputation and overwrite an
  existing fit (default FALSE). Note that if df.tps, batch, lambda.a, or
  gene.model are changed, the model is recomputed and overwritten.

- backend:

  a character, specifying the backend to use for computations (default
  'auto', see details). If 'gpu', GPU-based computations are used if
  available, otherwise CPU-based computations are used.

- verbose:

  a logical, specifying whether to show update messages (default TRUE).

- ...:

  other parameters fitting parameters.

## Value

a SpatialExperiment or Seurat object with the adjusted data stored in
'logcounts' or 'data', respectively.

## Details

SpaNorm works by first fitting a spatial regression model for library
size to the data. Normalised data can then be computed using various
adjustment approaches. When a negative binomial gene-model is used, the
data can be adjusted using the following approaches: 'logpac',
'pearson', 'medbio', and 'meanbio'.

The `df.tps` parameter specifies the degrees of freedom for the
thin-plate spline. If only 1 value is provided, it specifies the degrees
of freedom of the biology with the degrees of freedom of the library
size being half of that. If 2 values are provided, the first value
specifies the degrees of freedom of the biology and the second value
specifies the degrees of freedom of the library size. For rectangular
tissues, df.tps specifies the degrees of freedom along the length, with
the degrees of freedom along the width calculated ceiling(width / length
\* df.tps).

Similarly, the `lambda.a` parameter specifies the smoothing parameter
for regularizing regression coefficients. If only 1 value is provided,
it specifies the lambda.a for both the biology and library size
functions. If 2 values are provided, the first value specifies the
lambda.a for the biology and the second value specifies the lambda.a for
the library size. Batch effects are not regularised.

Batch effects can be specified using the `batch` parameter. If this
parameter is a vector, a design matrix will be created within the
function using `model.matrix`. If a custom design is provided in the
form of a numeric matrix, this should ideally be created using
`model.matrix`. The batch matrix should be created with an intercept
term. The SpaNorm function will automatically detect the intercept term
and remove the relevant column. Alternatively, users can subset the
model matrix to remove this column manually. Please note that the model
formula should include the intercept term and that the intercept column
should be subset out after.

## Examples

``` r
data(HumanDLPFC)
# \donttest{
SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#> Loading required namespace: SpatialExperiment
#> (1/2) Fitting SpaNorm model
#> 201 cells/spots sampled to fit model
#> Downloading uv...
#> Done!
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, estimating gene-wise dispersion
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, log-likelihood: -1137587.588303
#> iter:  1, fitting NB model
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  1, iter:  1, log-likelihood: -1137587.588303
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
#> iter:  1, iter:  2, log-likelihood: -812406.176544
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
#> iter:  1, iter:  3, log-likelihood: -727378.805684
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
#> iter:  1, iter:  4, log-likelihood: -713068.723490
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
#> iter:  1, iter:  5, log-likelihood: -711121.782393
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
#> iter:  1, iter:  6, log-likelihood: -710793.210862
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
#> iter:  1, iter:  7, log-likelihood: -710718.091129
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
#> iter:  1, iter:  8, log-likelihood: -710697.120927 (converged)
#> iter:  2, estimating gene-wise dispersion
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, log-likelihood: -710376.479496
#> iter:  2, fitting NB model
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> Hint: To use tensorflow with `py_require()`, call `py_require("tensorflow")` at the start of the R session
#> iter:  2, iter:  1, log-likelihood: -710376.479496
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
#> iter:  2, iter:  2, log-likelihood: -710235.662645
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
#> iter:  2, iter:  3, log-likelihood: -710227.877473 (converged)
#> iter:  3, log-likelihood: -710227.877473 (converged)
#> (2/2) Normalising data
#> class: SpatialExperiment 
#> dim: 5076 4015 
#> metadata(1): SpaNorm
#> assays(2): counts logcounts
#> rownames(5076): ENSG00000188976 ENSG00000188290 ... ENSG00000198727
#>   ENSG00000278817
#> rowData names(2): gene_name gene_biotype
#> colnames(4015): AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 ...
#>   TTGTTTCCATACAACT-1 TTGTTTGTGTAAATTC-1
#> colData names(3): cell_count sample_id AnnotatedCluster
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
# }
```
