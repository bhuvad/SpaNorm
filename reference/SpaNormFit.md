# An S4 class to store a SpaNorm model fit

An S4 class to store a SpaNorm model fit

## Usage

``` r
# S4 method for class 'SpaNormFit'
x$name
```

## Arguments

- x:

  an object of class SpaNormFit.

- name:

  a character, specifying the name of the slot to retrieve.

## Value

Return value varies depending on method.

## Slots

- `ngenes`:

  a numeric, specifying the number of genes in the dataset.

- `ncells`:

  a numeric, specifying the number of cells/spots in the dataset.

- `gene.model`:

  a character, specifying the gene-specific model to used (see
  `getGeneModels()`).

- `df.tps`:

  an integer, specifying the degrees of freedom to used for the thin
  plate spline.

- `sample.p`:

  a numeric, specifying the proportion of samples used to approximated
  the model.

- `lambda.a`:

  a numeric, specifying the shinkage parameter used.

- `batch`:

  a vector or matrix, specifying the batch design used (if any).

- `W`:

  a matrix, specifying the covariate matrix of the linear model.

- `alpha`:

  a matrix, specifying the coefficients of the linear model.

- `gmean`:

  a numeric, specifying the mean estimate for each gene in the linear
  model.

- `psi`:

  a numeric, specifying the over-dispersion parameter for each gene if a
  negative binomial model was used (or a vector of NAs if another gene
  model is used).

- `wtype`:

  a factor, specifying the covariate types of columns in the covariate
  matrix, W. These could be "biology", "ls", or "batch".

- `loglik`:

  a numeric, specifying the log-likelihood of the model at each external
  iteration.

- `sampling`:

  a factor, specifying the cells/spots used for dispersion estimation
  ('dispersion'), GLM fitting ('glm' and 'dispersion'), all other
  cells/spots ('all').

## Examples

``` r
example(SpaNorm)
#> 
#> SpaNrm> data(HumanDLPFC)
#> 
#> SpaNrm> ## No test: 
#> SpaNrm> ##D SpaNorm(HumanDLPFC, sample.p = 0.05, df.tps = 2, tol = 1e-2)
#> SpaNrm> ## End(No test)
#> SpaNrm> 
#> SpaNrm> 
#> SpaNrm> 
```
