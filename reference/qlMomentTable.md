# A shared table of NB deviance moments

The deviance moments `w0` and `w1` (see
[`nbDevianceMoments`](https://bhuvad.github.io/spaNorm/reference/nbDevianceMoments.md))
evaluated once on a uniform grid over log mean and log dispersion, for
[`qlDispersion`](https://bhuvad.github.io/spaNorm/reference/qlDispersion.md)
to interpolate onto every gene and cell. A caller that scores genes in
blocks builds one table over the whole range and passes it to every
block, so the result is invariant to how the genes are split; lookups
outside the range are clamped to the edge, where the moments are
asymptotically flat.

## Usage

``` r
qlMomentTable(
  lmu_range,
  lphi_range,
  step_mu = 0.08,
  step_phi = 0.12,
  ngrid = NULL,
  nphi = NULL
)
```

## Arguments

- lmu_range:

  range of log means the table must cover, `c(lo, hi)`.

- lphi_range:

  range of log dispersions, `c(lo, hi)`; a single value gives a one-row
  table.

- step_mu, step_phi:

  grid spacing along log mean and log dispersion; bilinear interpolation
  error scales with the square of the spacing.

- ngrid, nphi:

  explicit grid sizes, overriding the steps.

## Value

a list holding the two moment matrices (log-phi rows, log-mu columns)
and the grid geometry, for `qlDispersion(table = )`.

## Examples

``` r
tab <- qlMomentTable(log(c(1e-6, 1e3)), log(c(0.1, 10)))
dim(tab$w1)
#> [1]  40 261
```
