# Per-gene quasi-likelihood dispersion

The edgeR v4 quasi-likelihood dispersion of each gene's NB GLM: the sum
of its unit deviances, each rescaled by the observation's deviance
moments, over the effective residual degrees of freedom. With
`moments = "table"` (the default) the moments are evaluated once on a
shared (log mu, log phi) table and interpolated, so the per-gene work is
elementwise and runs on the accelerator when `y` and `mu` are torch
tensors.

## Usage

``` r
qlDispersion(
  y,
  mu,
  phi,
  design = NULL,
  p = NULL,
  prior = 1,
  leverage = c("trace", "exact"),
  moments = c("table", "grid", "cell"),
  ngrid = 256L,
  nphi = 64L,
  table = NULL
)
```

## Arguments

- y:

  counts, genes x cells (matrix or torch tensor).

- mu:

  fitted means, same shape (matrix or torch tensor).

- phi:

  NB dispersion: scalar or one value per gene.

- design:

  the design matrix (cells x p); needed for `leverage = "exact"`,
  otherwise only its column count is used.

- p:

  the number of design columns, in place of `design`.

- prior:

  the average quasi-dispersion edgeR divides through by; 1 leaves the
  parameterisation alone.

- leverage:

  `"trace"` spreads the p degrees of freedom evenly, which is exact to
  O(p/n) and is what n \>\> p designs want; `"exact"` forms
  per-observation hat values (O(n p^2) per gene, CPU only) and is what
  reproduces edgeR on its own small-n designs.

- moments:

  `"table"` evaluates the moments on a shared (log mu, log phi) table
  and interpolates onto every gene and cell (both backends); `"grid"`
  uses a per-gene log-mu grid; `"cell"` evaluates at every cell. The
  last two are CPU only.

- ngrid:

  grid points along log mu.

- nphi:

  maximum grid points along log phi for `moments = "table"`.

- table:

  a moments table from
  [`qlMomentTable`](https://bhuvad.github.io/spaNorm/reference/qlMomentTable.md),
  built once by a caller that scores genes in blocks; when `NULL` the
  table is built from this call's own range of means and dispersions.

## Value

a list of per-gene `deviance` (adjusted), `df` (effective) and
`s2 = deviance / df`.

## Examples

``` r
set.seed(1)
mu <- matrix(exp(rnorm(400, 0, 1)), 8, 50)
y <- matrix(rnbinom(400, mu = mu, size = 2), 8, 50)
qlDispersion(y, mu, phi = 0.5, p = 3)
#> $deviance
#> [1] 117.30917  97.58064 101.99442  96.83781 109.72176  77.87750 104.79576
#> [8]  77.78427
#> 
#> $df
#> [1] 94.67550 91.83868 91.35217 95.66567 91.40969 90.64134 89.48627 91.02706
#> 
#> $s2
#> [1] 1.2390658 1.0625223 1.1164969 1.0122525 1.2003296 0.8591831 1.1710819
#> [8] 0.8545181
#> 
```
